# Can frmtmb fit the models in eLife 56938?

Shinn, Lam & Murray (2020), *A flexible framework for simulating and fitting generalized
drift-diffusion models*, eLife 9:e56938 - the PyDDM paper. Model definitions and the
likelihood are pinned in `model-notes.md`; every claim below is backed by a probe script
in this directory.

**Short answer: yes, all five of the paper's models are reachable, but only two of them
are reachable by spelling alone.** The analytically tractable DDMs fit today through
`frmtmb.ddm::wiener()` and cost one line each. The paper's headline generalized DDM has
no closed-form density and needs a Fokker-Planck solve per likelihood evaluation - and
that solve turns out to be within reach of RTMB, which was the open question. A working
6-parameter GDDM fit, converging in under 3 seconds with all six parameters recovered
inside 1.4 standard errors, is in `probe-05-gddm-accuracy.R`.

## Mapping table

| Paper model | Verdict | Spelling / cost |
|---|---|---|
| **M1** 3-par DDM, Eq. 14 (without the lapse overlay) | FITS TODAY | `frm(bf(rt \| vint(upper) ~ 0 + coh, bias = 0.5), family = wiener(), data = d)` |
| **M1** as published, *with* the `0.95 p + 0.025` overlay | FITS WITH WORKAROUND | a 12-line `custom_family()`; see below |
| **M2** 8-par DDM, Eq. 16 (no overlay in the paper) | FITS TODAY | as M1 but `~ 0 + cohf` |
| **M4** 18-par DDM, Eq. 17 | FITS TODAY | `bf(rt \| vint(upper) ~ 0 + cohf, bs ~ 0 + cohf, ndt ~ 0 + cohf, bias = 0.5)` |
| **M3** 11-par full DDM, Eq. 15 | BLOCKED | needs `s_v`, `s_z`, `s_T`; a quadrature family, medium size |
| **M5** 6-par GDDM, Eq. 13 - the headline | BLOCKED, but proven feasible | needs a `gddm()` structured family in `frmtmb.ddm`; the solver is prototyped and works |

M1, M2 and M4 are the models the maintainer's anchor (1) predicted: constant drift,
constant bounds, point start. They reduce to the analytic Navarro-Fuss density, and
frmtmb's per-dpar regressions do the rest for free - M4 needs no new machinery at all,
just three formulas instead of one, which is a nicer construction than the paper's own
EZ-diffusion moment matching since it is full ML.

## What the probes established

**The density is exact.** `frmtmb.ddm`'s Navarro-Fuss implementation matches `RWiener`
to a maximum absolute log-density error of 1.4e-14 across drift, boundary, response and
`u = t/a^2` down to 0.0015 (`probe-03`). Nothing downstream is limited by the density.

**M1 recovers.** Over 25 replicates of 2640 trials at Roitman dimensions, bias is within
1.5 Monte Carlo standard errors on all three parameters, at 0.69 s per fit
(`probe-02`). On the real Roitman monkey 1 data (2611 trials after the paper's own
trimming), BIC orders the models M4 < M2 < M1, the saturated model beating the
3-parameter DDM by about 620 BIC - the same qualitative conclusion the paper draws, that
the standard DDM is not adequate for these data.

**The lapse overlay is not cosmetic.** With 6 percent contamination and no lapse
component, M1's drift estimate falls from a true 10.5 to 5.37, the boundary inflates from
1.6 to 2.08, and non-decision time collapses from 0.28 to 0.045 (`probe-06`). The
overlay is load-bearing, which is why treating "M1 without the overlay" as the paper's M1
would be wrong.

**The GDDM fits by gradient-based ML.** `probe-04` and `probe-05` build the 6-parameter
GDDM likelihood on an RTMB tape and fit it. AD gradients match central differences to
1e-9 on all six parameters. Recovery against the model's own density, n = 24000:

| | truth | estimate | se | z |
|---|---|---|---|---|
| mu0 | 8.00 | 8.018 | 0.072 | 0.25 |
| alpha | 0.80 | 0.800 | 0.006 | -0.05 |
| leak | 1.20 | 1.106 | 0.137 | -0.69 |
| B0 | 1.50 | 1.549 | 0.041 | 1.20 |
| tau | 1.00 | 0.984 | 0.016 | -1.05 |
| t_nd | 0.25 | 0.245 | 0.004 | -1.36 |

47 iterations, 2.6 s of optimization, max absolute gradient 2.3e-4 at the optimum, and
the estimates are stable to 3-4 decimals between a 200x201 and a 200x361 grid.

## The three things that had to change to get there

PyDDM's solver cannot be transcribed onto an AD tape, because three of its steps branch
on values. Each needed a replacement, and each replacement is the interesting part.

**1. The moving domain.** This is the hard one the maintainer flagged. The collapsing
bound shrinks the spatial grid, and PyDDM sandwiches the bound between two integer grids
with `weight_inner` / `weight_outer` (`pyddm/model.py`, `solve_numerical`). The weights
are continuous in `B0` and `tau` but kinked, and the grid indices themselves depend on
parameters. That is precisely why the paper reaches for differential evolution: its
objective is not differentiable, so no gradient method was ever an option for them.

The fix is to change variables rather than to chase the bound. With `y = x / B(t)` the
walls sit at +/-1 for all `t`, the grid is fixed, and no index depends on a parameter:

    dy = [ mu0 (C/Cmax)^alpha / B(t) - l y + y/tau ] dt + (1/B(t)) dW

Every coefficient is now a smooth function of every parameter, the two-grid sandwich
disappears, and a bonus falls out: the paper's Table 1 rules out Crank-Nicolson whenever
bounds are time-varying, so PyDDM's own GDDM fits run on backward Euler at `O(dt)`. After
rescaling, the bound is constant and Crank-Nicolson at `O(dt^2)` is available.

**2. The early exit.** PyDDM breaks out of the time loop once surviving mass drops below
1e-4. A tape cannot record that. Dropping it costs time, not correctness.

**3. The non-decision-time lookup.** Indexing the density at `round((rt - t_nd)/dt)`
makes an integer index depend on a parameter. Replaced by convolution with a cubic
B-spline of fixed length spanning the whole admissible `t_nd` range, written in
truncated-power form so it needs no branch: `max(x,0) = (x + |x|)/2` and `abs()` is
overloaded. C2 in `t_nd`, and a partition of unity so the shift conserves mass exactly.

## Two findings that transfer beyond this paper

**Renormalization is not optional.** The discretized solve loses a little probability
mass - 5.3e-4 at a 200x201 grid - and the loss *depends on the parameters*, because
configurations that absorb faster lose less. Maximizing an un-renormalized `sum(n log p)`
therefore pays a hidden bonus for fast absorption. Measured on data simulated from the
model itself, at the same grid:

| | truth | renorm TRUE | renorm FALSE |
|---|---|---|---|
| leak | 1.20 | 1.106 | 1.572 |
| B0 | 1.50 | 1.549 | 1.430 |

The bias is 3-4 standard errors and it lands squarely on the two parameters the paper's
scientific claim is about, leak and bound height. Dividing each condition's defective
density by its own total mass removes it entirely. Anyone fitting a GDDM by ML on a
discretized grid needs this, PyDDM included.

**The paper's coherence nonlinearity has an undefined gradient at zero coherence.**
`d/d(alpha)` of `C^alpha` is `C^alpha log(C)`, which is NaN at `C = 0` - and the Roitman
design contains a 0 percent coherence condition. PyDDM never trips over it because
differential evolution takes no derivatives. Any gradient-based implementation must
special-case `C = 0`, where the drift is identically zero and carries no `alpha`
dependence at all. `COH[j]` is data, so branching on it at tape-build time is legal.

## Which route for the GDDM: `frm_ode` or a structured family?

The maintainer asked whether method-of-lines on the Fokker-Planck could ride the existing
`frm_ode` / RTMBode seam. **It should not.** Three reasons, in descending order of how
fatal they are.

*Wrong seam.* `frm_ode()` produces a value for the linear predictor - it computes `mu`.
The GDDM needs the PDE solution to *be* the likelihood, evaluated at each trial's RT on
the correct boundary, as a defective bivariate density. There is no predictor to put the
answer in. This alone settles it.

*Stiffness.* Discretizing the Fokker-Planck operator in space gives eigenvalues spanning
roughly `[-4D/h^2, 0]`, and the rescaled diffusion coefficient `D = 1/(2 B(t)^2)` grows
as the bound collapses. At `B0 = 1.5`, `tau = 1`, `ny = 201`:

| t | B(t) | D | \|lambda\|max | explicit step limit |
|---|---|---|---|---|
| 0 | 1.500 | 0.222 | 9.1e3 | 2.2e-4 s |
| 1 | 0.552 | 1.642 | 6.7e4 | 3.0e-5 s |
| 2 | 0.203 | 12.13 | 5.0e5 | 4.0e-6 s |

A non-stiff integrator would need on the order of 495,000 steps over the 2 s window.
Crank-Nicolson is unconditionally stable and does it in 200.

*Adaptive stepping is not tape-safe.* An adaptive solver chooses step sizes by comparing
error estimates against tolerances - a value-dependent branch. Taping it freezes whatever
step sequence happened at tape-build time, so the recorded derivative is of a different
function than the one being evaluated. A fixed-step scheme has no such problem.

**The structured-family route is the right one**, with one nuance. The GDDM likelihood
*does* factorize over trials - it is a plain sum of per-row log densities, so on the face
of it this is `lpdf` territory, not `frmtmb_structure()` territory. The reason to reach
for a structure anyway is economy: the expensive object is shared by every row in a
condition. `lpdf` receives per-row dpar vectors and cannot discover that they are
constant within a condition, because that would mean comparing AD values. So a plain
`custom_family` would do 2611 PDE solves where 6 suffice. `frmtmb_structure()` with a
`frame_block` carrying the condition index does 6 and gathers - which is exactly the
economy PyDDM gets from `condition_combinations`.

The design wart to be aware of: the family would have to declare which dpars are allowed
to vary and group on a user-supplied grouping, because it cannot verify that, say,
`leak ~ subject` is constant within a condition. That is a documented constraint, not a
blocker.

## Cost, honestly

Measured at a 200x201 grid over 6 conditions (`probe-05`):

- objective: **0.4 ms** per call
- gradient: **28 ms** per call
- optimization: **2.6 s**, 47 iterations
- **tape build: 37-77 s**, and 65 s at PyDDM's default resolution (200x361)

Evaluation is cheap; the tape build is the bottleneck, and it scales with
`nt * ny * n_conditions` because RTMB records the Thomas sweep one scalar at a time
through an R-level loop. For the Roitman design this is a one-off minute per fit, which
is tolerable. For a design with many more conditions it would not be, and the escape
hatch is a compiled TMB template rather than R-level taping - a much larger undertaking
that should not be attempted before someone actually needs it.

## Features this queues

1. **`gddm()` structured family in `frmtmb.ddm`.** Solver, B-spline shift,
   renormalization, structure wrapper, family plumbing, tests, vignette. By analogy with
   the existing `wiener()` family, roughly 600-900 lines in the extension. **Needs no
   core grammar change** - that is the main result of this investigation. Everything the
   prototype needed already exists.

2. **Lift `wiener()`'s `max_ndt` guard when a lapse component is present.** Today
   `wiener(max_ndt = 0.4)` is refused outright whenever `min(rt)` is smaller, so
   `mixture(wiener(), contaminant)` cannot even be constructed. The refusal is correct for
   a bare Wiener model and wrong in a mixture, where the contaminant is precisely what
   covers the unreachable rows. PyDDM's tutorial fit puts `nondectime` at 0.211 s against
   a fastest RT of 0.203 s, so this is not a hypothetical corner.

3. **Make the Wiener density return `-Inf`, not `NaN`, below the non-decision time.**
   Measured in `probe-06`: it currently returns `NaN`, which would poison a log-sum-exp
   even if gap 2 were closed. Small fix, and it makes gap 2 safe to close.

4. **Ratcliff across-trial variability (`s_v`, `s_z`, `s_T`) for M3.** Gauss-Legendre
   quadrature inside an `lpdf`. `s_v` has a closed form; `s_z` and `s_T` need nodes. This
   is a self-contained family, not a grammar change.

5. **`init_dpars` validation.** Passing scalars instead of `function(y, aterms)` is
   accepted silently by both `custom_family()` and `check_custom_family()`, then fails
   deep inside `frm()` with `could not find function "init_fn"`. `check_custom_family()`
   is the natural place to catch it. Cost: a few lines.

6. Minor, observed in passing: `fitted(fit, dpar = "bs")` returns the same vector for
   every `dpar` - the argument appears to be silently ignored. Not chased further.

## Deliberate omissions

- **The Evans & Hawkins human dataset** (OSF 2vnam) was not fetched. It is pooled across
  subjects and the paper's own conclusion there is that the GDDM does no better than the
  standard DDM, so it adds nothing to a feasibility question that the monkey data has not
  already answered.
- **Monkey 2 (`monkey == 2`)** was not fit. Same models, more runtime, no new information.
- **The paper's own fitted parameter values were not reproduced, because they do not
  exist.** There is no parameter table, and no numeric BIC, log-likelihood or MSE anywhere
  in the article, its legends, its tables, its appendix or its back matter; the model
  comparison is unlabeled bar plots and the JATS carries no `source-data`. The only
  fitted numbers in the ecosystem are PyDDM's tutorial values for a *similar but not
  identical* model. This is why the prototype leans on simulation-based recovery.
- **M3 was not prototyped.** Its gap is a quadrature, which is well understood and
  uninteresting next to the GDDM question.
- **No PyDDM was installed or run.** Its source was read as the ground truth for the
  likelihood, which is what the brief asked for; running it would have added a Python
  toolchain for no extra evidence.
- **No timing comparison against PyDDM.** Different machine, different language, different
  optimizer. The numbers above are absolute costs, not a benchmark.

## Data provenance

`roitman_rts.csv` is redistributed in PyDDM's docs (MIT licensed) and derives from the
Roitman & Shadlen (2002) public release at the Shadlen lab. It was used only from the
session scratchpad and is **not** committed to this repository. Probe scripts locate it
via the `EL_SCRATCH` environment variable and skip the real-data sections when absent.

## Reproducing

    R_LIBS="<scratch>/el-lib;<user lib>" R CMD INSTALL -l <scratch>/el-lib extensions/frmtmb.ddm
    EL_SCRATCH=<scratch> Rscript dev/feasibility-elife56938/probe-0N-....R

from the worktree root. `probe-01` needs nothing but base R. `probe-04` and `probe-05`
need RTMB; `probe-02`, `probe-03`, `probe-06` and `probe-07` need `frmtmb.ddm` and
`RWiener`.
