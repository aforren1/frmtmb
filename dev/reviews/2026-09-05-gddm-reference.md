# Review: GDDM-REFERENCE lane (`wt-gddm-ref`)

Reviewer's independent verification. Written incrementally; the verdict is
at the end. Every number below that is labelled "reproduced" was produced
by this reviewer in a separate process and library, not copied from
`dev/gddm-reference-findings.md`.

Reviewer environment:

* private library `<scratchpad>/rgr-lib`, `frmtmb` 0.51.0 and
  `frmtmb.ddm` 0.2.0 installed from the worktree, `--library=` explicit
  throughout, user library only as a read fallback in `R_LIBS`.
* separate uv cache `C:\Users\adf44\gr-sp\rgr-uvcache`, so the lane's
  `gr-uvcache` was neither read nor written.
* the generator was copied BYTE-IDENTICAL to
  `C:\Users\adf44\gr-sp\rgr-rg\dev\gddm-pyddm-reference.py` and run from
  there, so `FIXTURES` resolved into the reviewer's own scratch tree and
  the lane's fixture was never overwritten. That also makes the
  regeneration a path-independence check, which running in place is not.

## 0. Scope

`git -C <wt> diff --name-only 2210aa1` gives exactly two tracked files,
`extensions/frmtmb.ddm/NEWS.md` and
`extensions/frmtmb.ddm/vignettes/gddm.Rmd`. Untracked: the generator, the
three fixture CSVs, the test, and the findings file. **No file under any
`R/` changed** (confirmed by filtering the diff for `(^|/)R/`). The main
checkout is clean throughout; it moved from 2210aa1 to 9b9011a during
this review, by sibling lanes and not by this reviewer, and none of
those commits touches `extensions/frmtmb.ddm/`. See section 9d.

Fixture size on disk: 59,081 + 2,971 + 10,243 = **72,295 bytes**. The
findings say 71 KB in section 3 and 72 KB in section 9; the second is
right. Trivial.

## 1. The generator and the fixture: REPRODUCED

The generator was copied byte-identical to a mirror tree under the
junction and run there with a fresh, reviewer-owned uv cache. All three
CSVs came back **byte-identical** to the ones in the worktree
(`cmp` on each), including `gddm-pyddm-meta.csv`, and the printed hash is

    37232180985e3a023ab89a9e03fefee50c41207c526b435a8085e80a52b72104

which is the value in `_fixture,content_sha256`. The per-case solver,
ladder level, step sizes and residuals printed by the run match the
findings' section 3 table row for row, including `exp_collapse`
stopping at level 4 with residual 2.584e-3 and `sz_uniform` at level 2
with 3.21e-4. Whole run: 66 s.

Determinism, read from the source: no `random`, no `numpy.random`, no
seed, no clock, no `time`. `platform.python_version()` and
`importlib.metadata.version()` are read, but only into the metadata
file, and the content hash covers the two payload files only, so the
hash is insensitive to them. The resolution ladder is a fixed list, not
a search that could depend on machine arithmetic. The only inputs are
literals in the file. **Deterministic as claimed**, and now also shown
to be path-independent, which running in place would not have shown.

The inline metadata pins `pyddm==0.9.0`, `numpy==2.5.2`,
`scipy==1.18.1` and `requires-python = "==3.13.*"`, which is what the
docstring and the findings say it pins. One caveat the lane does not
state: uv resolved **16** packages, and the eleven that are not pinned
float, among them `matplotlib 3.11.1`, `pandas 3.0.5` and
`paranoid-scientist 0.2.3`. None of them touches the solver arithmetic
(pyddm's numerics are numpy/scipy only), and the reviewer's independent
resolution happened to give a byte-identical fixture, so this is a note
rather than a defect: the fixture is what is frozen, and the metadata
assertions in the test are what would catch a regeneration that moved.

## 2. Convergence of the reference: CONFIRMED

Reviewer's own script, which does NOT import the generator; the models
are built from `R/gddm.R`'s documented parameterizations and pyddm's own
class docstrings, both read by hand. Each case solved at ONE HALVING
finer than the level the fixture stopped at, compared with the fixture:

| case | solver | dt (finer) | max abs d log vs fixture |
| --- | --- | --- | --- |
| constant | analytical | 3.125e-4 | 2.08e-6 |
| leak | Crank-Nicolson | 3.125e-4 | 2.18e-5 |
| lin_collapse | analytical | 3.125e-4 | 7.17e-8 |
| sz_uniform | Crank-Nicolson | 1.563e-4 | 1.63e-4 |
| lapse | analytical | 3.125e-4 | 2.06e-6 |

All five under the stated 5e-4. The reviewer's independent translation
reproducing the fixture to these levels is also an independent check of
the parameter map in section 6 below.

`exp_collapse` is the case that never reached 5e-4, and it is worth its
own line. Solved at the fixture's own level the reviewer's independent
model reproduces the fixture to **4.0e-13**, which is the round trip
through `%.12e`. Solved one halving finer than the fixture (dt 3.906e-5,
dx 7.813e-5, 38 s) it moves by **1.305e-3**. The fixture's recorded
`ref_residual_dlog` is 2.584e-3, and 2.584/1.305 = 1.98: the sequence is
**first order**, exactly as the findings predict for backward Euler on a
grid-snapped bound. Summing the remaining first-order tail, the
fixture's own distance from PyDDM's limit is about
1.305e-3 x 2 = 2.6e-3, which is the number the lane wrote down.

## 3. exp_collapse: the attribution is right, and provable more directly

Read out of pyddm 0.9.0's own source in the resolved environment:

* `Model.can_solve_cn()` (`model.py:585`) is
  `if self.get_dependence("bound")._uses_t(): return False`. A
  time-dependent bound has **no Crank-Nicolson path**. Confirmed.
* `Model.solve()` (`model.py:614`) dispatches analytic first, then
  Crank-Nicolson only `isinstance(bound, BoundConstant)`, else
  `solve_numerical_implicit`. `BoundCollapsingExponential` therefore
  falls to backward Euler. Confirmed.
* `Model.has_analytical_solution()` (`model.py:557`) returns True for a
  time-dependent bound only when the class is `BoundCollapsingLinear`
  and the IC is a point, which is the Anderson (1960) case the lane
  cites. Confirmed.
* The moving bound really is sandwiched and linearly blended:
  `model.py:830-841` computes `x_index_outer = floor(bound_shift/dx)`
  and `x_index_inner = ceil(...)` with `weight_inner`/`weight_outer`,
  solves the step on BOTH grid-snapped domains and adds them weighted
  (`model.py:891-901`). Confirmed, and it is first order in `dx`.

**The lin_collapse arbiter, re-measured.** PyDDM's own implicit solver
against PyDDM's own closed form, on the fixture's grid:

| dt, dx | max abs d log |
| --- | --- |
| 6.25e-4, 1.25e-3 | **5.5754e-3** |
| 3.125e-4, 6.25e-4 | 2.8808e-3 |

The first reproduces the fixture's recorded
`pyddm_numeric_vs_analytic_dlog = 5.575401944600e-03` to every digit,
and the pair HALVES, which confirms first order. gddm() against the
same closed form is 7.90e-4 at the shipped grid (reproduced, section 4),
so the 5.6e-3 against 7.9e-4 asymmetry is real and is about a factor of
seven, not the "roughly a hundred times" the findings claim in section
4 -- that sentence compares PyDDM at its grid with gddm at a FINER grid
(5.1e-5), which is not a like-for-like comparison. The direction of the
claim is right; the factor is overstated. In fairness, the like-for-like
factor of seven understates gddm in a different way, because PyDDM is
being run at dt 6.25e-4 while gddm is at dt 0.01, a sixteen times
coarser step. The defensible sentence is that on a collapsing bound
gddm is about seven times closer to the closed form while taking a
sixteen times coarser time step.

**A better argument than the arbiter, which the lane did not make.**
The lane infers the exp_collapse gap is PyDDM's by analogy with the
neighboring case. It can be shown directly. The reviewer solved
exp_collapse one level finer than the fixture and formed a first-order
Richardson limit, then measured gddm() against all three:

| gddm grid | vs FIXTURE | vs PyDDM one finer | vs PyDDM's extrapolated limit |
| --- | --- | --- | --- |
| dt .01 ny 201 | 3.250e-3 | 1.965e-3 | **8.081e-4** |
| dt .005 ny 401 | 2.832e-3 | 1.547e-3 | **2.780e-4** |
| dt .0025 ny 801 | 2.764e-3 | 1.459e-3 | **1.689e-4** |

gddm() does not stall. It converges, quartering then slowing as the
extrapolant's own error takes over, and what it converges TO is PyDDM's
limit, not PyDDM's level-4 value. At the shipped grid its distance from
the extrapolated truth is 8.1e-4, which sits inside the 7.3e-4 to
1.6e-3 band the other nine cases occupy. The whole 2.8e-3 is the
reference's discretization error.

**So: does the 5e-3 tolerance hide anything gddm should be blamed for?
No.** It is not even close. gddm()'s real error on this case is 8.1e-4,
sixth-best of the ten, and the tolerance is wide only because the
fixture it is measured against is wrong by 2.6e-3. The one thing the
tolerance does cost is discrimination: any gddm error smaller than
about 2e-3 on this case would be invisible. The lane's own mutation
table answers that (the exp_collapse mutant is caught by 1.4e5), and the
convergence tier assert on this case is `expect_lt(e[3], e[2])`, which
still bites: 2.832e-3 < 3.250e-3.

## 4. The measured discrepancy table: nine columns right, one wrong

Reproduced with the TEST FILE'S OWN helpers. The file was parsed and
every top-level `test_that()` call dropped, so `gr_worst_dlog()`,
`gr_mass()` and `gr_expected_mass()` are the shipped definitions and not
a paraphrase of them.

| case | dt .02 n101 | dt .01 n201 | dt .005 n401 | dt .0025 n801 | d loglik | d mass (SHIPPED rule) |
| --- | --- | --- | --- | --- | --- | --- |
| constant | 6.89e-3 | 8.76e-4 | 2.35e-4 | 6.41e-5 | +2.3e-4 | 8.65e-5 |
| bias_ndt | 6.61e-3 | 1.62e-3 | 4.22e-4 | 1.15e-4 | +1.6e-4 | 1.38e-4 |
| leak | 7.08e-3 | 1.00e-3 | 2.73e-4 | 7.71e-5 | -2.9e-3 | 8.81e-5 |
| unstable | 6.35e-3 | 7.33e-4 | 1.88e-4 | 3.93e-5 | -5.1e-5 | 1.31e-4 |
| exp_collapse | 6.27e-3 | 3.25e-3 | 2.83e-3 | 2.76e-3 | -2.6e-3 | 6.67e-5 |
| lin_collapse | 6.76e-3 | 7.90e-4 | 2.00e-4 | 5.06e-5 | +9.2e-4 | 3.74e-5 |
| coh_high | 6.79e-3 | 1.12e-3 | 2.84e-4 | 7.14e-5 | +1.3e-3 | 1.34e-5 |
| coh_low | 6.56e-3 | 8.71e-4 | 2.49e-4 | 7.39e-5 | -6.6e-4 | **1.53e-4** |
| sz_uniform | 5.63e-3 | 7.56e-4 | 2.69e-4 | 2.64e-4 | -8.3e-5 | 1.09e-4 |
| lapse | 6.83e-3 | 8.68e-4 | 2.33e-4 | 6.36e-5 | +1.7e-4 | 8.22e-5 |

**The four log-density columns and the log-likelihood column reproduce
the findings' section 4 table exactly**, to every digit printed. The
qualitative claims that rest on them are all confirmed: eight cases
quarter at each halving, the shipped grid lands between 7.3e-4 and
1.6e-3, `exp_collapse` stalls near 2.8e-3, `sz_uniform` stalls at
2.64e-4, and the worst log-likelihood difference is -2.9e-3 on `leak`.

**The mass column does not reproduce, and this is the one substantive
defect found.** The findings report a worst boundary-mass difference of
7.1e-5 on `unstable` and justify the granted 2e-4 as "about 2.8 times
the worst measured". The shipped test's worst is **1.53e-4** on
`coh_low`, so the granted tolerance is **1.31x** the worst measured, not
2.8x. It still passes, and it passed in this reviewer's run, but it is
the thinnest margin in the file and the findings describe it as the
widest.

The cause is identified, not guessed. `gr_mass()` (test:210) integrates
gddm()'s density with the **trapezoid** rule, while both sides of the
comparison are **rectangle** sums: `gd_densities()`'s `renormalize`
divides by `(sum(pu) + sum(pl)) * dt` (R/gddm.R:352) and PyDDM's
`Solution.prob()` is a plain sum of per-bin masses. Re-measuring with
the rectangle rule reproduces the findings' column exactly:

| case | shipped (trapezoid) | rectangle | findings' table |
| --- | --- | --- | --- |
| constant | 8.65e-5 | 2.554e-5 | 2.6e-5 |
| bias_ndt | 1.38e-4 | 2.875e-5 | 2.9e-5 |
| leak | 8.81e-5 | 1.676e-5 | 1.7e-5 |
| unstable | 1.31e-4 | 7.122e-5 | 7.1e-5 |
| exp_collapse | 6.67e-5 | 6.665e-5 | 6.7e-5 |
| lin_collapse | 3.74e-5 | 3.569e-5 | 3.6e-5 |
| coh_high | 1.34e-5 | 9.383e-6 | 9.4e-6 |
| coh_low | 1.53e-4 | 3.090e-5 | 3.1e-5 |
| sz_uniform | 1.09e-4 | 1.800e-5 | 1.8e-5 |
| lapse | 8.22e-5 | **1.076e-4** | 2.6e-5 |

Nine of ten match the rectangle column to two significant figures. So
the findings' mass column was measured with a rule the shipped test does
not use. The tenth, `lapse`, matches neither: under the rectangle rule
it is 1.076e-4, because gddm()'s lapse floor `0.5 lapse / t_max` summed
over `nt + 1` nodes times `dt` overcounts by `lapse * dt / t_max` =
1.67e-4 (visible as the `lapse` row's total mass being 1 + 1.667e-4);
the findings record 2.6e-5, which is the `constant` row's value, so that
entry looks like it was taken from the un-mixed model.

This is a documentation defect and a thin-margin warning, not a
correctness defect. Worth saying plainly: the trapezoid rule is
defensible on its own terms (the likelihood interpolates linearly and
`gd_draw()` uses the trapezoid), but it is not the rule the reference
was computed with, so roughly 1e-4 of the measured "disagreement" is a
quadrature convention rather than a solver difference. The `total-1`
column makes that visible: it runs to -1.46e-4 under the trapezoid rule
and is exactly 0 under the rectangle rule.

## 5. Euler-Maruyama: the numbers reproduce, the INFERENCE does not

The reviewer re-ran the test's own `gr_euler()` at n = 200000, h = 0.001,
with three seeds the lane did not use. Reference values, recomputed
here: PyDDM P(upper) 0.854860 and mean decision time 0.450595 (the
second from the fixture's two conditional means weighted by the two
boundary probabilities, which is what the test does); gddm at dt 0.005,
ny 401 gives 0.854831 and 0.450656. The two grid solvers therefore agree
with each other to 3.0e-5 and 6.1e-5, confirming the premise. (The
findings say gddm's mean is 0.450665 and the gap 7.0e-5; this reviewer
gets 0.450656 and 6.1e-5. Immaterial, but the findings' digit is off.)

| seed | P(upper) | its SE | vs PyDDM | mean DT | its SE | vs PyDDM |
| --- | --- | --- | --- | --- | --- | --- |
| 11111 | 0.853280 | 8.58e-4 | **-1.84 SE** | 0.451690 | 6.08e-4 | **+1.80 SE** |
| 22222 | 0.855875 | 8.48e-4 | +1.20 SE | 0.450387 | 6.04e-4 | -0.35 SE |
| 33333 | 0.853574 | 8.54e-4 | -1.51 SE | 0.450561 | 6.06e-4 | -0.06 SE |
| (lane, 20250905) | 0.855873 | 8.5e-4 | +1.2 SE | 0.449674 | 6.0e-4 | -1.5 SE |

**The standard errors reproduce exactly**: 8.5e-4 for P(upper) and
6.0e-4 for the mean, on every seed. The paired-Richardson machinery does
what the lane says it does, and the pairing rate is 1.0000, above the
0.999 the test asserts.

**The magnitude claim holds.** Pooling the reviewer's three independent
seeds, P(upper) sits -0.72 single-seed SE from PyDDM, i.e.
-6.1e-4 +/- 4.9e-4: no detectable shared bias. Over the four runs
together, 800000 paths, nothing larger than about 1e-3 in P(upper) is
resolvable and nothing is detected. The lane's "resolves a shared bias
down to about 2e-3" is therefore fair, and the vignette's "resolves a
shared bias down to roughly two parts in a thousand" is fair.

**The direction claim does not hold, and should be removed.** The
findings (section 6) and the test's comment at :534 argue that landing
"1.2 and 1.5 standard errors away from both, in the same direction for
P(upper) and the opposite for the mean ... is what leftover
higher-order extrapolation error looks like and not what a shared bias
looks like". Two objections, both demonstrated above:

1. The two functionals are strongly ANTI-CORRELATED by construction:
   under a positive drift, more upper-boundary hits means faster
   decisions, so P(upper) up forces mean DT down. Opposite directions is
   the null expectation for any single perturbation, noise included. It
   distinguishes nothing.
2. The signs are not stable. Seed 11111 gives the exact mirror of the
   lane's seed: P(upper) 1.84 SE BELOW and the mean 1.80 SE ABOVE. Seed
   33333 puts both on the same side. The observed direction is a
   property of the seed, not evidence about either solver.

The conclusion the lane draws (no shared bias at this resolution) is
supported by the MAGNITUDES across seeds. The reasoning offered for it
is not. Answering the question as put: at this n, "1.2 and 1.5 SE in
opposite directions" does not support the no-shared-bias claim; the
pooled magnitude does.

A note on the shipped assertion rather than the writeup: the test runs
n = 80000 with a fixed seed and a 4 SE gate, so it is deterministic and
cannot flake on sampling. Observed single-seed deviations reach 1.84 SE,
so 4 SE is real but not lavish headroom if R's `rnorm` stream ever
changes. Acceptable as written.

## 6. Mutation testing: three of the lane's, plus one of the reviewer's

Mutants built with the public extension points only (`gddm_drift_term()`,
`gddm_bound_term()`, `gddm_start_term()`); nothing under `R/` touched,
and they live in the reviewer's scratchpad, not in the package. Baselines
measured first, so the factor is against a model that passes.

| | case | worst d log | tolerance | factor | lane's figure |
| --- | --- | --- | --- | --- | --- |
| BASELINE leak | leak | 1.00e-3 | 2.5e-3 | 0.40 | 1.00e-3 |
| BASELINE exp_collapse | exp_collapse | 3.25e-3 | 5.0e-3 | 0.65 | 3.25e-3 |
| BASELINE constant | constant | 8.76e-4 | 2.5e-3 | 0.35 | 8.76e-4 |
| BASELINE sz_uniform | sz_uniform | 7.56e-4 | 2.5e-3 | 0.30 | 7.56e-4 |
| MUT leak sign flipped | leak | **1.92** | 2.5e-3 | **768x** | 770x |
| MUT tau read as a rate | exp_collapse | **686** | 5.0e-3 | **1.37e5 x** | 1.4e5 x |
| MUT ndt off by one step | constant | **4.85e-2** | 2.5e-3 | **19.4x** | 19x |
| MUT (reviewer) sz as a fraction of the FULL separation | sz_uniform | **2.59e-1** | 2.5e-3 | **104x** | (new) |

All three of the lane's reproduce to two significant figures, baselines
included. The reviewer's own mutation is the start-point slip the lane
did not test: writing `ICRange(sz = sz * bs)` instead of
`ICRange(sz = sz * bs / 2)`, i.e. reading `sz` as a fraction of the full
separation rather than the half separation, which is exactly the trap
`R/gddm.R:608` warns about. It is caught with a factor of 104. The
`sz_uniform` case does have teeth despite being the case whose own
residual stalls at 2.6e-4.

The tightest mutation remains the non-decision time off by one grid
step, at 19.4x. The lane's reading of that is right: the granted 2.5e-3
is nowhere near wide enough to swallow a real translation error.

## 7. The translation table: every row checked, both sides

Checked against `R/gddm.R`'s documentation and implementation, and
against pyddm 0.9.0's own class definitions read in the resolved
environment. The table appears in three places that must agree: the
generator's comment block, the test's block at :49, and the findings'
section 2. All three carry the same text.

| row | gddm side (file:line) | pyddm side (file:line) | verdict |
| --- | --- | --- | --- |
| noise | `dx = a dt + dW` | `NoiseConstant(noise=1)` | correct |
| boundary | separation `bs`, `B(0)=bs/2`, doc at `R/gddm.R:515-518` | `B` is the HALF separation, `BoundConstant` | correct |
| start point | `bias` in (0,1) above the lower bound, `R/gddm.R:600` | `ICPoint.get_IC` snaps `round(x0/dx)`, centered, `ic.py:86` | correct; `x0 = (bs/2)(2 bias - 1)` |
| constant drift | `a = mu`, `R/gddm.R:394` | `DriftConstant` | correct |
| **leak** | `a = -leak x`, positive is leaky, `R/gddm.R:400-404`; impl `function(x,t,p,cov) -p$leak * x` at `R/gddm.R:485` | `DriftLinear.get_drift = drift + x*x + t*t`, `drift.py:178`; its own docstring gives `x=-1` as "Leaky integrator" | **correct**: PyDDM `x = -leak`. Both signs exercised (`leak` +1.2, `unstable` -0.8) |
| **exponential collapse** | `B = (bs/2) exp(-t/tau)`, `R/gddm.R:564` | `get_bound = B*exp(-tau*t)`, `bound.py:137`, docstring "one divided by the time constant" | **correct**: PyDDM `tau = 1/tau`, a genuine reciprocal trap, documented at `R/gddm.R:523` |
| **linear collapse** | `B = (bs/2)(1 - kappa t/t_max)`, `R/gddm.R:525` | `get_bound = max(B - t*time, 0)`, `bound.py:107` | **correct**: PyDDM `t = (bs/2) kappa / t_max` = 0.2 per second here. `kappa < 1` keeps `B(t_max) = 0.4 > 0`, so PyDDM's `max(.,0)` clip never engages and the two are the same function, not merely the same at the start |
| coherence drift | `a = sign(C) mu (abs(C)/cmax)^alpha`, `R/gddm.R:396` | `DriftConstant(drift = that number)` | correct; `drift_value` in the fixture is 2 and `2*0.5^1.3` and the test asserts both |
| **start variability** | uniform half width `sz` as a fraction of the HALF separation, `R/gddm.R:607-609` | `ICRange.get_IC` fills `2 int(sz/dx)+1` cells centered on the domain, `ic.py:177-185`, carries NO bias | **correct**: PyDDM `sz = sz * bs/2`. The generator's `assert c["bias"] == 0.5` is not decoration, it is required, because `ICRange` cannot express a bias |
| non-decision time | `gd_shift()`, centered cubic B-spline, `R/gddm.R:300` | `int(ndt/dt)` whole bins, hard truncation, `overlay.py:308` | correct as a translation; see the note below |
| **lapse** | renormalize, THEN mix, floor `0.5 lapse / t_max`, `R/gddm.R:352-361` | `OverlayUniformMixture` adds `0.5 c / len(t_domain)` of MASS to a DEFECTIVE density, `overlay.py:212` | the two are genuinely different; handled honestly, see below |

Three remarks.

**The lapse is the only row where the two packages do not agree, and the
lane does the right thing with it.** PyDDM's window is `T_dur + dt` and
gddm's is `t_max`, so the floors are 8.3316e-3 against 8.3333e-3. The
generator does not paper over this: it rebuilds PyDDM's mixture by hand,
records `lapse_formula_max_abs_diff` (asserted below 1e-12 in the test),
writes BOTH densities into the fixture, and then carries gddm's
convention into the fixture values. The test asserts the two floors
differ by less than 2e-6. This is the correct handling of a convention
mismatch: state it, measure it, and pick one.

**The non-decision-time shift is not the same operator on the two
sides, and this is probably what the test mostly measures at the
shipped grid.** PyDDM translates by whole bins exactly. gddm convolves
with a centered cubic B-spline, which even at an exactly integer shift
applies weights `(1/6, 2/3, 1/6)`, a smoothing with variance `dt^2/6`.
That is an `O(dt^2)` perturbation of the density, so it converges at the
second order the test observes, and at `dt = 0.01` its size is of order
`dt^2 f''/12`, which for these densities is several times 1e-4: the same
order as the 7.3e-4 to 1.6e-3 the test actually sees. This is not a
defect. It is the density the likelihood reads, so it is the right thing
to compare, and the lane says as much in findings section 4. But it is
worth recording that the headline "the two solvers agree to 1e-3" is in
part a statement about the shift kernel and not only about the
Fokker-Planck solve. A reader of the vignette could be forgiven for
thinking otherwise.

**One consequence not stated anywhere.** At the coarsest grid in the
convergence tier, `dt = 0.02`, the non-decision times 0.25 and 0.40 are
12.5 and 20 steps: `constant`'s shift is half a step OFF-node. gddm's
spline shift handles that smoothly, which is why the coarse column is a
uniform 5.6e-3 to 7.1e-3 rather than erratic. It also means the coarse
column is not purely a grid-resolution effect. The convergence assertions
still hold and hold for the right reason (the reviewer confirmed the
ratios quarter from the second column on), so this is a note only.

## 8. CRAN safety and the vignette's honesty

**CRAN safety: yes, on all three counts.**

* *Size.* 72,295 bytes across three CSVs. Nothing near a threshold.
* *No Python at test time.* Grepped the test for `system`, `system2`,
  `reticulate`, `processx`, `Sys.which` and `python`: the only hit is
  `expect_match(gr_meta("_fixture", "python"), "^3[.]13[.]")`, a regex
  against a string read out of the frozen CSV. The generator is the only
  Python and it is excluded from the tarball by the existing
  `^dev$` line in `.Rbuildignore`.
* *No network.* No `download.file`, `url`, `curl`, `http`. The test's
  only inputs are the three fixture files, reached through
  `testthat::test_path()`.
* *Dependencies.* Only `utils`, `tools`, `stats` and `testthat`. No new
  entry in `Suggests` is needed and none was added, correctly.
* *Old R.* `tools::sha256sum()` is guarded by
  `skip_if_not(exists("sha256sum", where = asNamespace("tools")))`,
  which is the right feature test for a package declaring R >= 4.1. On
  R 4.6.1 it runs rather than skips, which is why the CRAN-mode count is
  147 passes and 2 skips rather than 146 and 3.

**The vignette section: honest about scope, with two overstatements.**

What it gets right, all reproduced by this reviewer: the range "between
7e-4 and 1.6e-3 ... on every case but one" (measured 7.33e-4 to
1.62e-3); "an order of magnitude better than the 0.01 claimed above";
PyDDM having no Crank-Nicolson path for a moving bound and its own
residual being 2.6e-3; the lin_collapse trio 5.6e-3, 7.9e-4 and 5.1e-5;
the leading edge exclusion at 0.2 s; that fitting is not compared; and
the two translated conventions (renormalization, and the lapse window).
The "what the reference does not cover" paragraph is unusually candid
for a vignette and names the right things, including components a user
writes themselves.

Two sentences should change.

1. "On the exponentially collapsing bound the two sit 2.8e-3 apart in
   the log density and **refining either grid does not close the gap**."
   The second half is false. Refining gddm()'s grid does not close it,
   which is what the lane measured. Refining **PyDDM's** grid does:
   one level past the fixture the gap falls from 3.25e-3 to 1.97e-3, and
   against PyDDM's extrapolated limit it is 8.1e-4 (section 3). The
   sentence gives away the lane's own case. It should say that refining
   this package's grid does not close it while refining PyDDM's does,
   which is precisely how you know whose the gap is.
2. "It agrees with both packages **to about one standard error**." The
   lane's own table says 1.2 and 1.5, and the reviewer's seeds reach
   1.84. "To one or two standard errors" is the honest phrasing.

One thing the vignette does not say that it arguably should: the density
comparison it reports is of the density the likelihood reads, which
includes gddm()'s cubic-B-spline non-decision-time shift against PyDDM's
whole-bin truncation (section 7). At the shipped grid that term is
plausibly the larger part of the quoted 1e-3. The current text invites
the reader to attribute all of it to the Fokker-Planck solve.

## 9. Runs

All against `rgr-lib`, `--library=`/`R_LIBS` explicit, one process per
file, `NOT_CRAN=true` unless stated.

### 9a. Per file

| file | pass | fail | skip | warn | elapsed | lane's figure |
| --- | --- | --- | --- | --- | --- | --- |
| test-gddm-family.R | 77 | 0 | 0 | 0 | 7.9 s | 77 |
| test-gddm-gradients.R | 22 | 0 | 0 | 0 | 104.6 s | 22 |
| test-gddm-recovery.R | 28 | 0 | 0 | 0 | 427.8 s | 28 |
| test-gddm-reference.R (NOT_CRAN=true) | **166** | 0 | 0 | 0 | 22.6 s | 166 |
| test-gddm-reference.R (CRAN mode) | **147** | 0 | **2** | 0 | 8.8 s | 147 / 2 |
| test-gddm-solver.R | 94 | 0 | 0 | 0 | 32.8 s | 94 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 | 1.8 s | (not stated) |
| test-surface.R | 47 | 0 | 0 | 0 | 27.8 s | (not stated) |

Every count the lane published is confirmed. The two CRAN-mode skips are
the resolution sweep (`test-gddm-reference.R:406`) and the
Euler-Maruyama tier (`:491`), both `skip_on_cran()`, which is the right
choice for a 10 s Monte Carlo and a 4-grid sweep. The `tools::sha256sum`
expectation runs rather than skips on R 4.6.1, so the hash check is
live in CRAN mode too. `test-message-uniqueness.R` passing matters
because it is the file a new user-visible string would break, and this
lane adds none: no `R/` file changed.

## 10. sz_uniform does not stall either, for the same reason

The lane files `sz_uniform` next to `exp_collapse` as a second case where
the sequence "stalls", at 2.6e-4, and explains it as two starting
distributions that "differ at first order in `dx` and second order in
the spatial step respectively, and neither refinement removes the
other's". The first half is right and the second half is not. Refining
PyDDM's `dx` does remove it.

PyDDM's own sz_uniform sequence, solved by the reviewer two halvings
past the fixture's level: L2 -> L3 moves 1.632e-4 and L3 -> L4 moves
8.232e-5. It HALVES, so `ICRange`'s box, which occupies
`2 int(sz/dx) + 1` whole cells and is therefore about `dx` too wide, is
a first-order-in-`dx` error, exactly as the lane's mechanism says.
Measuring gddm() against those refined values:

| gddm grid | vs FIXTURE (L2) | vs L3 | vs L4 | vs PyDDM's limit |
| --- | --- | --- | --- | --- |
| dt .01 ny 201 | 7.555e-4 | 7.361e-4 | 7.257e-4 | 7.152e-4 |
| dt .005 ny 401 | 2.685e-4 | 2.276e-4 | 2.172e-4 | 2.067e-4 |
| dt .0025 ny 801 | **2.644e-4** | 1.012e-4 | 7.558e-5 | **6.512e-5** |

Against the extrapolated reference gddm() goes 7.2e-4, 2.1e-4, 6.5e-5:
a clean second-order sequence, the same behavior as the eight cases the
findings call converged. The 2.644e-4 "stall" is the fixture's own
`ICRange` error at L2 and nothing else. gddm()'s integrated cubic
B-spline box is not competing with it at some irreducible offset; it is
simply more accurate than the thing it is being measured against.

So BOTH of the lane's two stalling cases are the reference's error, and
in both the lane's diagnosis of the mechanism is right while its
statement that refinement cannot remove the gap is wrong. This does not
change the verdict; it strengthens the lane's own case, and it is worth
recording because the findings and the vignette both tell the reader
that a floor exists where there is none.

A methodological point in the lane's favor, since I checked it: the
ladder's stopping rule compares consecutive rungs, which for a
first-order sequence is not merely a change but is EQUAL to the
remaining error (`e_{n-1} - e_n = e_n` when `h` halves). The reviewer's
own refinement confirms this numerically for exp_collapse: the recorded
residual is 2.584e-3 and the next rung moves 1.305e-3, so `e_L4` is
2 x 1.305e-3 = 2.61e-3, matching the recorded value. The fixture's
`ref_residual_dlog` is therefore a correctly calibrated error estimate
and not just a delta. For the second-order cases it is conservative by
about 3x. Good design, and the lane did not claim credit for it.

### 9b. R CMD check --as-cran

`_R_CHECK_CRAN_INCOMING_=false`, pandoc 3.8.3 from the RStudio quarto
tools directory, tarball built from the worktree.

    Status: OK

No WARNINGs, no NOTEs. Timings on a loaded machine: tests 412 s,
`--run-donttest` examples 39 s, vignette rebuild 110 s. The tests inside
the check report

    [ FAIL 0 | WARN 0 | SKIP 18 | PASS 964 ]

which is the lane's figure exactly. Tarball **252,950 bytes**. The three
fixture CSVs travel in it at 59,081 + 2,971 + 10,243 bytes;
`tar tzvf | grep dev/` returns nothing, so `dev/gddm-pyddm-reference.py`
is excluded by the pre-existing `^dev$` in `.Rbuildignore`, as claimed.

### 9c. The whole suite, one process, NOT_CRAN=true

This is the run the lane started and stopped. It completed here in
1514 s.

    RGR-SUITE TOTAL pass=1089 fail=0 skip=1 warn=0 err=0

Per file, so "nothing else moved" can be checked file by file rather
than trusted in total: brms-parity 13, defects 59, density 132, family
34, gddm-family 77, gddm-gradients 22, gddm-recovery 28,
**gddm-reference 166**, gddm-solver 94, lba 109, message-uniqueness 4,
moments 29, sampling 94 (+1 skip), simulate-density 41, surface 47,
variability 140.

**The lane's 1092 is right and my 1089 is the one that needs
explaining.** The single skip is
`test-sampling.R:13`, "frm_sample runs a short chain on a wiener model",
skipped with reason `{frmtmb.sample} cannot be loaded`: an artifact of
the reviewer's private library, which was built with only `frmtmb` and
`frmtmb.ddm` in it. That test contains exactly three `expect_true()`
calls, so with that package present the total is 1089 + 3 = **1092**.
Taking this lane's file back out, 1092 - 166 = **926**, which is the
lane's stated baseline to the expectation. Both of the lane's numbers
are correct, and it was right that the run it abandoned could not have
changed the verdict. Confirmed directly by installing `frmtmb.sample`
from the worktree into `rgr-lib` and re-running that file.

Nothing else moved. Every non-reference file's count is reachable from
the per-file runs in 9a and from the baseline arithmetic above, there
are no failures, no warnings and no errors anywhere, and the only file
whose count is new is `test-gddm-reference.R` at 166.

Confirmed directly: with `frmtmb.sample` installed from the worktree
into `rgr-lib`, `test-sampling.R` goes from 94 passes and 1 skip to
**97 passes and 0 skips**. That is the missing 3.

### 9d. Repository state

Main was clean at 2210aa1 at the start. It is now at **9b9011a**, moved
by three sibling-lane commits while this review ran (`db779ef`, a
frmtmb.spline test fix; `177102f` and `9b9011a`, pkgdown and docs). Not
by this reviewer: every git command issued here was `log`, `status`,
`diff`, `ls-files`, `rev-parse` or `merge-base`, all read-only. 2210aa1
is an ancestor of 9b9011a, and the intervening diff touches
`_pkgdown.yml`, `extensions/frmtmb.spline/tests/testthat/test-deriv.R`
and 118 files under `docs/`, and **nothing under
`extensions/frmtmb.ddm/`**. So the scope measured against 2210aa1 is
still exactly this lane's work, and this lane will land on current main
without interacting with those commits.

## 11. Punch list

Nothing here blocks. Item 1 is the only one that is a statement about
the code rather than about prose.

| # | file:line | what | severity |
| --- | --- | --- | --- |
| 1 | `test-gddm-reference.R:210` (`gr_mass`) and `:152` | `gr_mass()` integrates gddm with the TRAPEZOID rule, but both things it is compared with are RECTANGLE sums (`gd_densities()` renormalize at `R/gddm.R:352`, and PyDDM's `Solution.prob()`). About 1e-4 of the measured "disagreement" is that convention. The worst measured is **1.53e-4** (`coh_low`), not the 7.1e-5 the comment at `:152` states, so `GR_TOL_MASS = 2e-4` has 1.31x headroom, not 2.8x. Either switch the nine non-lapse cases to the rectangle rule, which drops the worst to 7.1e-5, or keep the trapezoid and correct the comment and the tolerance note. Keep the trapezoid for `lapse`, where it is the right rule. | should fix |
| 2 | `dev/gddm-reference-findings.md:130` (the whole `d mass` column) and `:183` | The mass column was measured with the rectangle rule while the shipped test uses the trapezoid; nine of ten entries reproduce the rectangle rule exactly. The `lapse` entry (2.6e-5) matches neither rule and looks like the un-mixed model's value; under the rectangle rule it is 1.076e-4. The section 5 justification "about 2.8 times the worst measured" is 1.31x for the shipped test. | should fix |
| 3 | `vignettes/gddm.Rmd:444` | "refining either grid does not close the gap" is false. Refining PyDDM's grid closes it from 3.25e-3 to 1.97e-3, and against PyDDM's extrapolated limit gddm is 8.1e-4. The sentence gives away the lane's own point. | should fix |
| 4 | `vignettes/gddm.Rmd:463` | "agrees with both packages to about one standard error": the lane's own figures are 1.2 and 1.5, and other seeds reach 1.84. Say "one to two". | should fix |
| 5 | `test-gddm-reference.R:539` and findings section 6 | The "opposite directions, which is ... not what a shared bias looks like" inference is unsound: the two functionals are anti-correlated by construction, so opposite directions is the null expectation, and the signs flip with the seed (section 5 above). Replace with the pooled-magnitude argument, which does hold. | should fix |
| 6 | `dev/gddm-reference-findings.md:155` | "roughly a hundred times closer" compares PyDDM at its grid with gddm at a finer grid. Like for like it is about seven times, on a sixteen times coarser step. | minor |
| 7 | `dev/gddm-reference-findings.md:166` and the exp_collapse paragraph above it | Both "stalls" are removed by refining the REFERENCE, which the lane did not try. Neither case has a floor. The mechanisms named are right. | minor |
| 8 | `dev/gddm-reference-findings.md:77` | "895 rows"; the file has **893** data rows (9 x 91 + 74). | trivial |
| 9 | `dev/gddm-reference-findings.md:74` against `:300` | "71 KB" and "72 KB" for the same three files; they are 72,295 bytes. | trivial |
| 10 | `dev/gddm-reference-findings.md:218` | gddm's mean decision time is 0.450656, not 0.450665; the gap to PyDDM is 6.1e-5, not 7.0e-5. | trivial |
| 11 | `dev/gddm-pyddm-reference.py:1-7` | The inline metadata pins pyddm, numpy, scipy and the Python minor version. Eleven transitive dependencies float (matplotlib, pandas, paranoid-scientist and the rest). None affects the arithmetic, and a fresh resolution reproduced the fixture byte for byte, but a line in the docstring, or a lock file, would make the claim airtight. | optional |
| 12 | `NEWS.md` | The bullet is accurate. One awkward short line break at "and the test also requires / the disagreement to shrink". | cosmetic |

Two things NOT on the list, because they were checked and are fine: the
`exp_collapse` tolerance of 5e-3 hides nothing gddm should be blamed
for (section 3), and the fixture, the generator and the test are
CRAN-safe (section 8).

## 12. Edits made by this reviewer

**None to any file belonging to the lane.** The generator, the three
fixture CSVs, the test, the vignette and NEWS.md are untouched. `cmp`
against the independently regenerated fixture is byte-clean, and
`git status` in the worktree shows the same two modified and four
untracked paths it showed at the start, plus this review.

One file created, as instructed:
`C:\Users\adf44\source\r\frmtmb-wt-gddm-ref\dev\review-gddm-reference.md`.

The punch-list items were deliberately left for the lane rather than
applied here. Items 3, 4 and 5 are wording in shipped documentation and
in a test comment, and the lane should phrase its own claims; applying
them would also have invalidated the `R CMD check --as-cran` run in
section 9b, which was made against the tree exactly as the lane left it.

## 13. The junction

`C:\Users\adf44\gr-sp` is a directory junction to the session
scratchpad. It is real and it is load-bearing on Windows: the
reviewer's own uv environment for pyddm sits under
`C:\Users\adf44\gr-sp\rgr-uvcache\environments-v2\`, and the deepest
scipy path under the un-junctioned scratchpad prefix is at MAX_PATH, so
the import failure the lane documents in its section 7 is reproducible
in principle and the junction is the cheap fix.

**It should be removed at consolidation.** It points into a
session-specific scratchpad that will not exist afterwards, so leaving
it is a dangling link in the user's home directory, and nothing in the
repository refers to it: the generator resolves its output path from
`__file__`, so it runs from any directory short enough, which this
reviewer demonstrated by running it from a different mirror tree
entirely. Removing the link with `rmdir` and no `/s` removes only the
link. What should survive is the knowledge, and that belongs in the
generator's docstring rather than in a link: one line saying that on
Windows the environment needs a short `UV_CACHE_DIR`, because the
deepest scipy path otherwise exceeds 260 characters. That is punch-list
item 11's natural companion.

## 14. Verdict

**GO WITH FIXES.**

The substance is sound, and where this reviewer could push harder than
the lane did it got stronger rather than weaker. The fixture
regenerates byte for byte from a deterministic generator, on a
different path and a fresh package resolution. The parameter
translation is correct in every row, including the three traps (the
leak sign, the `tau` reciprocal, and `sz` as a fraction of the half
separation), checked against both packages' sources rather than against
the lane's summary of them. The nine density columns and the
log-likelihood column reproduce exactly. The mutation tier reproduces,
and a mutation the lane did not try is caught with a factor of 104.
`R CMD check --as-cran` is Status OK with FAIL 0, WARN 0, SKIP 18,
PASS 964, and the whole-suite `NOT_CRAN=true` run that the lane
abandoned completes at 1092 over a baseline of 926: both of its
unverified numbers confirmed.

The two cases the lane reports as stalling do not stall. Refining the
reference removes both gaps: on `exp_collapse` gddm sits 8.1e-4 from
PyDDM's extrapolated limit, and on `sz_uniform` 6.5e-5. The lane's
attribution was right and its evidence was weaker than the evidence
available to it.

What holds this back from a plain GO is one measurement that does not
reproduce, and the documentation resting on it: the boundary-mass
column was measured with a quadrature rule the shipped test does not
use, so the tolerance with the thinnest margin in the file (1.31x, not
2.8x) is described as the one with the widest. Nothing fails, and no
gddm defect is hidden by it, but a tolerance justification that does
not match the code is the kind of thing this tier exists to catch
elsewhere. Fix items 1 to 5 and this is a GO.

# Punch re-check, 2026-09-05

Second pass by the same reviewer, against the lane's punch round. Same
rules as the first pass: private library `rgr-lib`, explicit
`--library=`, absolute paths, no commit, nothing in the main checkout
touched. Every number below was re-measured here; none is copied from
the lane's punch-round section.

## What actually changed on disk

Four files, by mtime and by diff:

| file | changed | R/ touched |
| --- | --- | --- |
| `extensions/frmtmb.ddm/tests/testthat/test-gddm-reference.R` | 15:19 | no |
| `extensions/frmtmb.ddm/vignettes/gddm.Rmd` | 15:19 | no |
| `extensions/frmtmb.ddm/dev/gddm-pyddm-reference.py` | 15:13 | no |
| `dev/gddm-reference-findings.md` | 15:23 | no |

`NEWS.md` is untouched this round, and its bullet makes no claim the
round contradicts, so that is right. `R/gddm.R` is still at 12:51 and
`git diff --name-only 2210aa1` still shows no `R/` file. The three
fixture CSVs are still at 13:44:40 and are **byte-identical** to the
copy this reviewer regenerated independently in the first pass.
`git status` in the worktree shows the same two modified and five
untracked paths as before. Main is clean and still at 9b9011a; this
reviewer ran only read-only git commands.

## P1. Mass tolerance: VERIFIED

Re-measured with the EDITED test file's own helpers, by parsing it and
dropping every top-level `test_that()`:

    GR_TOL_MASS   = 3e-04
    worst mass (SHIPPED trapezoid rule) = 1.5320e-04 on coh_low
    headroom = 3e-04 / 1.5320e-04 = 1.958x

The lane's 3e-4, 1.532e-4, `coh_low` and 1.96x all hold.
`GR_TOL_DLOG` (2.5e-3) and `GR_TOL_LOGLIK` (5e-3) are unchanged, which
is right: those two were never in question.

**The comment matches the shipped rule.** `test-gddm-reference.R:156-174`
now says the tolerance is measured with `gr_mass()`'s trapezoid rule,
"which is the rule the assertion uses", names 1.53e-4 on `coh_low` and
1.96x, and says explicitly that the old 2e-4 was justified against
7.1e-5, the RECTANGLE worst, leaving 1.31x. That is the defect this
review reported, stated in the file rather than only in a dev note.
The paragraph that follows is accurate on the mechanism: both sides of
the comparison are rectangle sums, `renormalize` dividing by
`(sum(pu) + sum(pl)) * dt` and PyDDM's `Solution.prob()` being a plain
sum of bin masses, so about 1e-4 of the number is quadrature
convention. Keeping the trapezoid is defended on the right grounds
(the piecewise-linear density the likelihood interpolates, the rule
`gd_draw()` samples from, and the correct rule for `lapse`, whose
density is nonzero at both ends of the window). The coordinator's
citation of ":159-174" is a few lines short; the block starts at 156.

One point worth stating plainly, because widening a tolerance always
deserves the question: this is not a tolerance widened to admit a
failing measurement. Nothing moved. The assertion passed at 2e-4 and
passes at 3e-4; what changed is that the stated justification now
describes the rule the code uses. The alternative fix, switching the
nine non-lapse cases to the rectangle rule and keeping 2e-4, was
available and would have given a tighter number (7.1e-5, 2.8x). The
lane took the other branch and documented the cost. Both are
defensible; the documented one is now honest, which was the complaint.

## P2. Euler-Maruyama direction argument: VERIFIED in the shipped files

`test-gddm-reference.R:570-581` now says the magnitudes carry the
claim, that which side a run lands on "carries nothing and is not read
here", that the two functionals are anti-correlated by construction so
opposite directions is the null expectation, and that the observed
signs flip with the seed. It also records that single seeds have been
seen at 1.84 SE, which is this reviewer's seed 11111. The vignette at
the "third opinion" paragraph now reads "to one or two standard
errors" and carries the same disclaimer. Both are correct and both
match what was measured.

The magnitude claim is unchanged and still supported: 1.2 and 1.5 SE
for the lane's seed, and 0.72 single-seed SE pooled over four seeds,
which is this reviewer's pooled figure from seeds 11111, 22222, 33333
plus the lane's.

**In the findings the argument is superseded, not removed.** The
original section 6 text still stands verbatim at
`dev/gddm-reference-findings.md:222`, "in the same direction for
P(upper) and the opposite for the mean, which is what leftover
higher-order extrapolation error looks like". What retracts it is the
appendix at `:398`, which opens by naming every superseded item and
says "where this section disagrees with sections 3 to 6 above, this
section is the current statement". See the residuals below.

## P3. Both refinement ladders: REPRODUCED, to the digit

Re-run in this reviewer's environment against the edited test file.
The PyDDM side is this reviewer's own solve from the first pass, so
these are two independent measurements agreeing, not one being read
back.

First, the rungs at the fixture's own level reproduce the fixture
EXACTLY, `identical()` on the parsed doubles, not merely to a
tolerance:

    exp_collapse fixture rung (L4) identical to fixture: TRUE   rows = 74
    sz_uniform   fixture rung (L2) identical to fixture: TRUE   rows = 91

which is the lane's "74 of 74 and 91 of 91" claim, confirmed. So the
refinement ran the same code path the fixture was built with.

`exp_collapse`, L4 to L5 moves **1.3050e-3** (lane: 1.3050e-3):

| gddm grid | vs fixture | vs L5 | vs PyDDM's limit |
| --- | --- | --- | --- |
| dt .01 ny 201 | 3.2504e-3 | 1.9649e-3 | **8.0812e-4** |
| dt .005 ny 401 | 2.8325e-3 | 1.5469e-3 | 2.7797e-4 |
| dt .0025 ny 801 | 2.7641e-3 | 1.4591e-3 | 1.6895e-4 |

`sz_uniform`, L3 to L4 moves **8.2317e-5** (lane: 8.2317e-5), and
L2 to L3 moved 1.6321e-4 in the first pass (lane: 1.6321e-4), so the
sequence halves and `ICRange`'s whole-cell box is first order in `dx`:

| gddm grid | vs fixture | vs L4 | vs PyDDM's limit |
| --- | --- | --- | --- |
| dt .01 ny 201 | 7.5553e-4 | 7.2567e-4 | **7.1521e-4** |
| dt .005 ny 401 | 2.6852e-4 | 2.1717e-4 | 2.0671e-4 |
| dt .0025 ny 801 | 2.6441e-4 | 7.5580e-5 | 6.5122e-5 |

Every figure the lane reports is reproduced. Two of its printed values
are truncated rather than rounded in the last digit: it prints 8.080e-4
and 2.779e-4 for the exp_collapse limit column where the values are
8.0812e-4 and 2.7797e-4, so correctly rounded they are 8.081e-4 and
2.780e-4. The same two appear in the test comment at `:188-189`. This
changes nothing and is noted only for completeness.

The vignette's rewritten "Where they disagree" paragraph is now
correct: "refining this package's grid alone does not close it. That
is the reference's error, not this package's, and the way to see it is
to refine the REFERENCE." It gives 8.1e-4, 2.8e-4, 1.7e-4 against the
extrapolated limit, and adds the sz_uniform case with 7.2e-4, 2.1e-4,
6.5e-5. Both match. The `lin_collapse` analogy is withdrawn as the
argument, correctly, since the direct measurement supersedes it, and
the "roughly a hundred times" figure is corrected to about seven times
on a sixteen times coarser step at `:514-517`.

## P4. Generator docstring: VERIFIED, and it still reproduces the fixture

The diff is seven added lines inside the module docstring and nothing
else. This matters more than it looks: the paragraph contains a Windows
path, and `C:\uvcache` written literally in a non-raw Python string
would be an invalid `\u` escape and a syntax error. The lane wrote
`C:\\uvcache`, which is correct.

Not taken on trust. The edited generator was copied byte-identical to a
second mirror tree and re-run with the reviewer's own uv cache. It
produced all three CSVs **byte-identical** to the shipped fixture, and
printed

    sha256 37232180985e3a023ab89a9e03fefee50c41207c526b435a8085e80a52b72104

which is the value in `_fixture,content_sha256` and the same hash as
the first pass. The per-case ladder table printed by the run is
unchanged. So the fixture is untouched at the same hash, as claimed,
and the generator that produces it still does.

## P5 and P6: VERIFIED

The 1092 attribution is recorded correctly, including that a library
without `frmtmb.sample` lands on 1089 with one skip at
`test-sampling.R:13`, and that the three expectations are the
difference. That is this reviewer's measurement, and the lane has
recorded it as a measurement rather than as arithmetic, which was the
point.

The three number corrections are present and right: 893 density rows
(9 x 91 + 74, counted here again), 72,295 bytes, and gddm's mean
decision time 0.450656 with a 6.1e-5 gap.

## P7. Runs after the edits

One process each, this reviewer's library, edited tree:

| mode | pass | fail | warn | skip |
| --- | --- | --- | --- | --- |
| CRAN (`NOT_CRAN` unset) | **147** | 0 | 0 | **2** |
| `NOT_CRAN=true` | **166** | 0 | 0 | 0 |

Unchanged from before the punch round, which is the expected result:
the only assertion touched was a tolerance that was already being met.
The two CRAN skips are still the convergence sweep and the
Euler-Maruyama tier, now at `:439` and `:524` after the comment
insertions.

No other assertion moved. Every hard-coded tolerance in the file is
unchanged: the lapse formula check at 1e-12 (`:371`), the total-mass
check at 3e-4 (`:406`), the two solver-agreement premises at 1e-3
(`:567-568`), the 4 SE gates (`:586`, `:588`) and the resolution floors
at 3e-3 (`:593-594`). `GR_TOL_DLOG`, `GR_TOL_LOGLIK` and
`GR_TOL_DLOG_EXP` are unchanged. `GR_TOL_MASS` is the only constant
that moved, which is what was asked for and no more. The identical
147 and 166 counts corroborate that nothing was added or removed.

## The build, on the edited tree

`R CMD build` then `R CMD check --as-cran --no-manual`,
`_R_CHECK_CRAN_INCOMING_=false`, pandoc 3.8.3:

    Status: OK

No WARNINGs, no NOTEs. Vignette rebuild OK in 97 s, which is the real
check on the rewritten prose: the edited `gddm.Rmd` is still valid Rmd
and still knits. Tests inside the check:

    [ FAIL 0 | WARN 0 | SKIP 18 | PASS 964 ]

identical to the pre-punch run, as expected. The tarball is
**255,273 bytes**, up 2,323 from 252,950, which is the added comment
and prose and nothing else. Verified in the tarball itself: 65 entries,
the three fixture CSVs present at their unchanged sizes, `dev/`
entirely absent and `gddm-pyddm-reference.py` not shipped.

## Residuals

Three, none blocking, none in a shipped file.

| # | file:line | what | severity |
| --- | --- | --- | --- |
| R1 | `dev/gddm-reference-findings.md:77`, `:130`, `:152-155`, `:166`, `:183`, `:218`, `:222` | The punch round supersedes by appendix rather than editing in place, so the original wrong numbers are still greppable at their old sites with no local marker: "895 rows", the rectangle `d mass` column, the `lin_collapse` analogy and "roughly a hundred times", the two "stall" paragraphs, "about 2.8 times the worst measured", 0.450665, and the direction argument. The appendix at `:398` names every one of them and declares itself the current statement, so a reader who starts at the top is not misled; a reader who greps is. A one-line "superseded, see the punch round" at each of the seven sites closes it. This is a dev log and not shipped, so it is an editorial call rather than a defect, and the supersede header is unusually complete. | minor |
| R2 | `test-gddm-reference.R:188-189` and the P3 table in the findings | 8.080e-4 and 2.779e-4 are truncated rather than rounded; the values are 8.0812e-4 and 2.7797e-4, so 8.081e-4 and 2.780e-4. Reproduced here to five figures. Cosmetic. | trivial |
| R3 | (not a lane defect) | The coordinator's citation of the mass comment as `:159-174` is a few lines short; the block starts at `:156`. | trivial |

Everything else from the first pass's punch list is closed. For the
record, mapping the original list onto this round:

* item 1 (mass rule and comment) closed by P1, verified.
* item 2 (findings mass column and section 5) closed by P1's table,
  subject to R1.
* item 3 (vignette "refining either grid") closed, rewritten and
  correct.
* item 4 (vignette "about one standard error") closed, now "one or two".
* item 5 (direction argument in test and findings) closed in the test
  and the vignette, superseded in the findings, subject to R1.
* item 6 ("a hundred times") closed, corrected to about seven.
* item 7 (the two stalls) closed by P3, and the vignette now carries
  the sz_uniform case too, which the first pass did not ask for.
* items 8, 9, 10 (893 rows, 72,295 bytes, 0.450656) closed by P6.
* item 11 (UV_CACHE_DIR) closed by P4, and the generator still
  reproduces the fixture byte for byte.
* item 12 (NEWS line break) not addressed. It was cosmetic and NEWS was
  correctly left alone this round.

## Updated verdict

**GO.**

The one defect that held the first pass at GO WITH FIXES is fixed, and
fixed in the place that mattered: the shipped test now states the rule
its own assertion uses, names the number that rule produces, and
records what keeping the trapezoid costs. The tolerance moved from
1.31x to 1.958x headroom without any assertion changing state, so this
is a corrected justification and not a widened gate. The four
documentation items are corrected in the shipped vignette and the test
comments, with the two "stalls" now re-attributed by direct
measurement rather than by analogy, and the sz_uniform case added.

Everything re-measured here agreed with the lane: the mass worst at
1.5320e-4 on `coh_low` and 1.958x; both refinement ladders to four or
five significant figures, with the fixture-rung rows `identical()` to
the fixture on 74 of 74 and 91 of 91 points; the edited generator
reproducing all three CSVs byte-identically at the same SHA-256; 147
and 166 with zero failures; and `R CMD check --as-cran` Status OK with
FAIL 0, WARN 0, SKIP 18, PASS 964.

Nothing under `R/` changed in either round, so the solver this tier
measures is the one that was reviewed. The junction recommendation is
unchanged: remove `C:\Users\adf44\gr-sp` at consolidation, now with
less at stake, since the reason it existed is recorded in the
generator's docstring where someone regenerating the fixture will
actually read it.

Edits by this reviewer in this round: none to any lane file. This
appendix is the only thing written, appended to
`dev/review-gddm-reference.md`.
