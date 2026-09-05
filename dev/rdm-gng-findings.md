# RDM and go/no-go: working notes

Lane `wt-rdm-gng`, worktree `C:/Users/adf44/source/r/frmtmb-wt-rdm-gng`,
branch `wt-rdm-gng`, base 2210aa1 (frmtmb 0.51.0, frmtmb.ddm 0.2.0).
Written as the work happens; measurements are pasted where they were
taken, not summarized from memory.

## Reference environment

Private library
`.../scratchpad/rg-lib`: frmtmb 0.51.0 (from this worktree),
frmtmb.ddm 0.2.0, rtdists, EMC2 3.5.0 (CRAN binary, pulled WienR,
psych, corrplot, lpSolve, GPArotation, magic, matrixcalc).
R 4.6.1.

## What EMC2 actually computes

Read out of the installed namespace rather than the paper.

### RDM

`EMC2:::rRDM` draws `bs <- B + runif(n, 0, A)` and passes `bs` to
`rwaldt(n, k = bs, l = v)`, which is an inverse Gaussian with
`mu = k / l` and `lambda = k^2`. So:

* `bs` is the DISTANCE the accumulator has to travel, not the
  threshold. It is uniform on `(B, B + A)`.
* Therefore EMC2's `B` is the LBA's `k`: the gap between the top of
  the start-point range and the threshold, `b = A + B`.
* Unit diffusion after `pars[, c("A","B","v")] / s`, so `s` is the
  diffusion coefficient and dividing by it is the scale convention.

`EMC2:::dRDM` and `pRDM` are `dWald` / `pWald` (compiled), guarded by
`rt > t0 & !v < 0`: a negative drift is scored as a density of zero
rather than refused.

### Go/no-go DDM

`EMC2:::log_likelihood_ddmgng` is the whole model in six lines:

```r
isna <- is.na(dadm$rt)
ok <- attr(pars, "ok") & !isna
like[ok] <- model$dfun(dadm$rt[ok], dadm$R[ok], pars[ok, , drop = FALSE])
ok <- attr(pars, "ok") & isna
like[ok] <- pmax(0, pmin(1, (1 - model$pfun(dadm$TIMEOUT[ok], dadm$Rgo[ok],
                                            pars[ok, , drop = FALSE]))))
```

So a go row contributes the ordinary DEFECTIVE DDM density at the
observed boundary, and a no-go row contributes
`1 - F_go(TIMEOUT)`, one minus the defective distribution function of
the GO boundary at the deadline. That is the derivation the task
states, confirmed against the implementation: the no-go mass is the
lower-boundary mass accumulated by the deadline plus the mass still
diffusing at it.

`TIMEOUT` is a per-row column of the augmented data, so the deadline
is per-row data, which is what an addition term carries.

## RDM: the derivation, and what it measures

One accumulator is a unit-diffusion Wiener process with drift `v`
travelling a distance drawn uniformly on `(k, k + A)`, where `k` is the
gap between the top of the start-point range and the threshold. For a
fixed distance the first passage is inverse Gaussian; averaging over
the distance integrates in closed form because

    G(c) = v Phi((c - v t)/sqrt(t)) - phi((c - v t)/sqrt(t))/sqrt(t)

is an antiderivative in `c` of the inverse-Gaussian density (checked by
differentiating: dG/dc = c phi(x)/t^(3/2), which is the density). So
the density is `(G(k + A) - G(k))/A`, and it has EXACTLY the algebraic
shape of the LBA density with the drift standard deviation `s` replaced
by `1/sqrt(t)`. That is why `lba_phidiff()` and `lba_race_lpdf()` are
reused rather than re-derived: only the single-accumulator law changes.

The survival needs a second antiderivative, of the inverse-Gaussian
distribution function, and that one carries a `1/(2 v)`. The
singularity is removable: the bracket it divides is the change across
the distance range of a quantity that is identically one at zero
drift.

### Measured, package code, 648-point grid

Reference is an 80-point Gauss-Legendre average over the distance of
the POINTWISE cancellation-free density and survival. Pointwise the
survival is written `Phi(x) * -expm1(logE - logPhi(x))`, which never
subtracts, so the reference is trustworthy far below where the
subtractive forms die.

| reference value above | density max rel | survival max rel |
|---|---|---|
| 1e-3   | 1.55e-13 | 1.33e-12 |
| 1e-9   | 1.35e-12 | 1.26e-11 |
| 1e-15  | 3.01e-12 | 6.41e-11 |
| 1e-100 | 7.53e-12 | 6.40e-10 |
| 1e-250 | 1.27e-10 | 4.09e-08 |

No exact zeros in either, on any row.

### Against EMC2

On the 490 grid rows where `EMC2:::dWald()` is above 1e-12, EMC2 and
this family differ by up to 7.95e-4 relative. Adjudicated against the
non-subtracting reference on those same rows, this family is 3.0e-12
wrong and EMC2 is 7.95e-4 wrong: the gap is EMC2's subtractive
`pnorm(x2) - pnorm(x1)`, not a disagreement about the model.

`1 - EMC2:::pWald()` returns EXACTLY ZERO on 35 of 315 rows of a
comparable grid, the first at a survival near 1e-13. This family
returns positive values down to 1e-250 on the same rows. Same defect,
same shape, as `1 - rtdists::plba_norm()` in `lba.R`.

### The one place it is only good to parts per million

The `1/(2 v)` divides a bracket whose rounding error is absolute, so
the survival's relative error grows as the drift goes to zero. Measured
at `A = 0.5`, `k = 1`, `t` in `{0.5, 2}`:

| v | max rel |
|---|---|
| 1e-1 | 4.9e-16 |
| 1e-3 | 1.91e-13 |
| 1e-5 | 1.39e-11 |
| 1e-7 | 6.07e-10 |
| 1e-9 | 1.94e-07 |

A drift of 1e-9 is an accumulator that never finishes, where the
survival is 1 to eleven places anyway, so this is stated rather than
fixed. The log link on the drifts keeps a fit off zero.

## Go/no-go: the Wiener defective CDF this package did not have

`wiener()` declares no `lcdf` and its compat table says so. The no-go
branch needs one, so it is written here. Two series, as the density
has, and blended the same way.

Small time, the image sum of Blurton, Kesselmeier and Vollrath (2012),
derived by integrating the small-time density term by term: each image
term is an inverse-Gaussian density at a signed distance
`c_j = a (w + 2 j)`, so each integrates to an inverse-Gaussian
distribution function. The integration constant differs by the SIGN of
`c_j`, which is decided by `j` alone and not by any parameter, so the
branch is on the loop index and stays tape-safe.

Large time, one series for the no-go probability directly:

    nogo(t) = P_lower
              + 2 pi exp(v a (1 - w)) sum_k (-1)^(k+1) k sin(k pi w)
                       exp(-lambda_k t) / (v^2 a^2 + k^2 pi^2)

with `lambda_k = (v^2 a^2 + k^2 pi^2) / (2 a^2)` and `P_lower` the
gambler's-ruin probability of ever reaching the lower boundary. This
form matters: it computes the no-go probability DIRECTLY rather than as
`1 - F_upper`, so it does not cancel where a no-go trial is surprising.

`P_lower` is written `exp(-v a w) (1 - w) sinhc(v a (1 - w)) / sinhc(v a)`
rather than as the textbook ratio of exponentials, which is 0/0 at zero
drift; `sinhc(x) = sinh(x)/x` is smooth there, and `ddm_pos()` keeps its
argument off exactly zero.

### Verified at 260 bits

The two series are independent derivations, and at 260 bits they agree
with each other to between 8e-64 and 7e-77 depending on the row, the
spread being how much the small-time route's `1 - F_up` step cancels
there. (The 4.8e-78 first recorded here is a BEST row, not a worst
case; the reviewer's independent pair agreed to between 1.6e-64 and
6.0e-76, which is the same statement.) Against that truth, in double
precision:

| u = t / a^2 | small-time K=12 | large-time K=30 |
|---|---|---|
| (0, 0.01]    | 0        | 9.5e-06 |
| (0.01, 0.02] | 1.5e-16  | 1.4e-13 |
| (0.02, 0.05] | 9.1e-16  | 1.1e-14 |
| (0.05, 0.1]  | 1.0e-12  | 2.2e-15 |
| (0.1, 0.2]   | 1.8e-15  | 3.9e-16 |
| (0.2, 0.35]  | 1.2e-06  | 9.5e-16 |
| (0.6, 1]     | 6.2e-11  | 8.7e-16 |
| (1, 3]       | 6.8e-13  | 5.4e-16 |
| (3, 10]      | 1.3e-14  | 4.2e-16 |
| (10, Inf)    | 1.7e-04  | 2.7e-16 |

The small-time route's 1.2e-06 at u in (0.2, 0.35] is not truncation:
it is the `1 - F_upper` cancellation on a row whose no-go probability
is 1.9e-10.

Blend weight `lam = 0.5 (1 + tanh((log u - log u0) / us))`, measured:

| u0 | us | max rel over the grid |
|---|---|---|
| 0.03 | 0.12 | 1.03e-14 |
| 0.06 | 0.12 | 2.21e-15 |
| 0.10 | 0.12 | 7.62e-13 |
| 0.35 | 0.12 | 1.20e-06 |
| 0.06 | 0.30 | 8.87e-11 |

**SUPERSEDED at the punch round; see the "Item 1" section at the end of
this file.** Every number in that table is real but the GRID is too
narrow: its relative start point stops at 0.7, and the constant's worst
corner needs a start point strongly biased toward the go boundary. On a
1200-row grid reaching `w = 0.9`, `u0 = 0.06` is 7.75e-05, not
2.21e-15, and the shipped centre is now 0.02.

The reasoning that chose a low centre in the first place stands, and is
the reason the fix is one constant rather than a rewrite: the density's
own centre of 0.35 is the wrong one here, because the large-time route
is valid much further down for the CDF than it is for the density, and
taking it earlier is what avoids the `1 - F` cancellation. The error
was stopping at 0.06 on the evidence of a grid that could not see where
the small-time route collapses.

### WienR, hence EMC2, is the less accurate of the two

At `t = 2.5, a = 4, v = 5, w = 0.75` the no-go probability is
1.1883283058299342e-13 at 260 bits.

CORRECTED at the punch round; the first version of this section said
`WienR` was 4.4 percent wrong "at every `precision` setting from 1e-8
to 1e-16", and that is false at the top of the range. Measured:

| `WienR` precision | `1 - pWDM` | relative error |
|---|---|---|
| default | 1.1357581542e-13 | 4.4239 % |
| 1e-8 | 1.1357581542e-13 | 4.4239 % |
| 1e-10 | 1.1357581542e-13 | 4.4239 % |
| 1e-12 | 1.1357581542e-13 | 4.4239 % |
| 1e-14 | 1.1868284133e-13 | 0.1262 % |
| 1e-16 | 1.1879386363e-13 | 0.0328 % |

`WienR::WienerCDF()` behaves identically to `pWDM()`.

The substantive point survives and is stronger stated correctly:
`EMC2:::pDDM` and `EMC2:::dDDM` both default to `precision = 0.005`,
looser than every setting in that table, so EMC2 in practice sits at
the 4.4 percent end rather than at WienR's best. This family's
large-time route is 3.0e-15 at that point.

What this correction also exposed is FINDING 1 below: the SHIPPED
blend is 5.0e-09 at that same point, not the 1e-12 the notes claimed,
because the blend weight leaks onto the cancelling small-time route.
The large-time route was never the problem.

## Identity against EMC2, end to end

Run by `extensions/frmtmb.ddm/dev/rdm-gng-emc2-reference.R`, which is a
script and not a test on purpose: every EMC2 function that computes
either likelihood is internal, so a package test asserting on them
would be asserting on another package's private surface. The suite
uses exported references instead (statmod, WienR, RWiener) and this
script pins the EMC2 agreement.

One trap, recorded because it cost real time: `EMC2:::dWald` and
`EMC2:::pWald` are Rcpp and DO NOT RECYCLE. Passing `B` or `A` as a
scalar alongside a length-n `t` returns quiet garbage rather than an
error, and the first version of the race comparison read that garbage
as a 27.9 disagreement in log units.

### go/no-go, composed exactly as `EMC2:::log_likelihood_ddmgng` does

400 simulated trials, 290 go and 110 no-go, at
`mu = 1.0, bs = 1.4, ndt = 0.25, bias = 0.45, deadline = 1.5`:

| rows | max absolute log-likelihood difference |
|---|---|
| go (290)   | 6.11e-16 |
| no-go (110) | 1.33e-15 |
| total log likelihood | 1.71e-13 (mine -195.5478210393, EMC2 -195.5478210393) |

### RDM, composed exactly as `EMC2:::log_likelihood_race` does

`log d(winner) + sum log(1 - p(losers))`, 200 rows each:

| accumulators | rows | max absolute log-likelihood difference |
|---|---|---|
| 2 | 200 | 6.41e-13 |
| 3 | 200 | 1.03e-12 |
| 4 | 200 | 1.12e-12 |

### RDM single-accumulator density, adjudicated by statmod

statmod is the third party here: it knows nothing about either package,
and its inverse Gaussian averaged over the start point is what the
model says.

| `dWald` above | EMC2 vs mine | mine vs statmod | EMC2 vs statmod |
|---|---|---|---|
| 1e-3  | 4.46e-13 | 1.26e-13 | 4.22e-13 |
| 1e-8  | 6.65e-08 | 4.56e-13 | 6.65e-08 |
| 1e-12 | 2.14e-04 | 4.56e-13 | 2.14e-04 |

So where the two disagree, the disagreement is EMC2's.

### RDM survival

| `1 - pWald` above | max rel |
|---|---|
| 1e-3  | 4.88e-10 |
| 1e-6  | 5.27e-09 |
| 1e-9  | 1.34e-06 |
| 1e-12 | 1.11e-04 |

`1 - EMC2:::pWald()` returns EXACTLY ZERO on 35 of 315 grid rows; this
family returns none. The first is `t=15, v=2, k=0.3, A=0.1`, where this
family gives 2.165e-16.

The suite's own version of that finding is sharper and does not need
EMC2: written subtractively the survival STICKS at 8.88e-16 for
`t` in `{3, 5, 10, 20}` at `v=8, A=0.1, k=0.3`, reporting the same
number at every one, while the true value falls from 2.8e-44 to
9.4e-282. Wrong by more than 250 orders of magnitude, and the
tell is that it does not move.

## A defect the compat probe found, and the fix

`fitted()`, `predict(type = "response")` and
`residuals(type = "response")` were NOT refused by either family as
first written. frmtmb falls back to "the first primary dpar on the
response scale is the mean" when a family declares no `post$mean_fn`,
and for these two that fallback returns a DRIFT RATE:

* `rdm(3)` on data whose response times average 0.36 returned 3.51
  from `predict(type = "response")`, and `residuals(type = "response")`
  came back as a length-zero vector rather than refusing.
* `wiener_gng()` was worse, because its primary dpar is literally
  called `mu`: `fitted()` returned a constant 1.05 for data whose go
  response times average 0.6, with nothing to signal that the number
  was a drift.

Both families now declare a `post$mean_fn` that stops with a reason.
The refusal is the right answer on the merits and not merely a safer
one, especially for the go/no-go family: a go/no-go trial's outcome is
a PAIR, and the no-go rows have no response time to average at all.

Two consequences worth knowing. `residuals(type = "pearson")` now
reports the MEAN's refusal rather than the variance function's,
because pearson needs the mean first. And the compat rows for
`fitted`, `predict` and `residuals` were written from these runs
rather than from what the families were meant to do.

## Compat rows, measured rather than reasoned about

Every row was run by a probe that calls the feature and records what
came back. Several first guesses were wrong, and in both directions:

| pair | guessed | measured |
|---|---|---|
| rdm / quadrature | refused | WORKS, on a model with a random effect |
| rdm / REML | untested | WORKS |
| rdm / mixture | untested | REFUSED by frmtmb: components need a dpar called `mu` and this family's are `v1..vn` |
| rdm / residuals_osa | untested | refused, but inside RTMB with a type error, not by a sentence |
| rdm / predict | "vint() mandatory on newdata" | not true for `type = "link"`: no linear predictor reads it |
| rdm / dec() | "refused by the coercion" | was SILENTLY IGNORED; now refused by the family |
| wiener_gng / quadrature | refused | WORKS |
| wiener_gng / REML | untested | WORKS |
| wiener_gng / mixture | untested | assembles and runs, but the one case tried did not converge (false convergence, max gradient 7.7e13), so `conditional` |

The quadrature row is worth keeping in mind next to `wiener()`'s, which
records a refusal. frmtmb refuses `quadrature = TRUE` only when there is
no random-effect block to marginalize; wiener's row was measured on a
model without one. Both rows are right about what they measured.

### The dec() divergence from lba()

`lba()` does not refuse `dec()`. Measured: `rt | dec(two) + vint(choice)`
fits under `lba(3)` with the `dec()` term silently dropped. `rdm()`
refuses instead, which is a deliberate divergence between two sibling
race families. The reason is that a silent drop is the worse failure for
a model ported over from `wiener()`, where `dec()` is the spelling. If
the divergence is unwelcome the fix is to make `lba()` refuse too, not
to make `rdm()` ignore.

## Verification runs

Private library `rg-lib`: frmtmb 0.51.0 built from this worktree,
frmtmb.ddm 0.3.0, plus EMC2 3.5.0, WienR, statmod, rtdists, RWiener,
Rmpfr, numDeriv as references.

One trap in the harness, recorded because it produced six false
failures: `testthat::test_file()` runs in the caller's environment, so a
test that reaches an internal UNQUALIFIED (as `test-defects.R` does)
cannot see it. The runner has to pass
`env = new.env(parent = asNamespace("frmtmb.ddm"))`, which is what
`R CMD check` does. Nothing was wrong with the package.

| file | result |
|---|---|
| `test-rdm-gng.R` (new) | 149 pass, 0 fail, 0 skip |
| `test-simulate-density.R` | 69 pass (was 63 before the two new families) |
| `test-lba.R` | 109 pass |
| `test-surface.R` | 47 pass |
| `test-message-uniqueness.R` | 4 pass |
| `test-defects.R` | 59 pass |

All with `NOT_CRAN=true`, so the recovery and Monte Carlo tests really
ran rather than skipping.

### The core's test-compat.R with this extension loaded

267 pass, 1 fails, and the failure is PRE-EXISTING rather than mine:
`declared families exist in the family registry` checks every declared
family against `frmtmb:::family_registry`, which is core's list of
built-in families. A `custom_family()` never enters it. Measured, the
missing set is `wiener, gddm, lba, rdm, wiener_gng`: the three families
already shipped at 0.2.0 fail the same assertion, so the invariant does
not hold for ANY contributed family and my two rows did not break it.
Every other invariant in that file passes over the enlarged registry,
which is what proves the new rows resolve.

## Left out, and why

* **Across-trial variability on `wiener_gng()`.** `wiener()` offers
  `sv`, `sz` and `st0`. The go branch would inherit all three for
  free, because it is `wiener()`'s density. The no-go branch would
  not. The drift enters the DENSITY as an exponential-quadratic, which
  is what makes the `sv` integral exact and free; it enters the
  DISTRIBUTION FUNCTION through the eigenvalues of both series, where
  it is not, so `sv` would need a quadrature of its own. `sz` and
  `st0` could be reached with the existing Gauss-Legendre nodes, but
  shipping two of three would make `variability =` mean something
  different on this family than on `wiener()`, which is worse than not
  having it. EMC2's `DDMGNG` does carry all three, by integrating them
  numerically in compiled code.
* **A fitted mean for either family.** Refused rather than
  approximated; see the defect note above. For `rdm()` it is a missing
  closed form, for `wiener_gng()` it is a statement about the model.
* **`lcdf` for either family.** Neither declares one, so `cens()` and
  `trunc()` are refused by frmtmb. The go/no-go family DOES now have a
  Wiener distribution function, but it is the probability of no
  response by a deadline, not the response-scale CDF `cens()` needs.
* **A start-point-free `rdm()`.** `A` divides the density, so `A = 0`
  is a limit rather than a fitted value, exactly as for `lba()`. A
  separate pure shifted-Wald race would be numerically cleaner (no
  `1/A`, no `1/(2v)`) but it is a second family, not an argument.
* **`sd_v`-style scale argument on `rdm()`.** The diffusion
  coefficient is fixed at 1 and is not exposed, where `lba()` exposes
  `sd_v`. EMC2 exposes `s` and divides `A`, `B` and `v` by it, which is
  a reparameterization rather than a different model, so nothing is
  lost; a user who wants EMC2's `s` divides their own numbers.
* **`lba()` was left alone.** It ignores a `dec()` term rather than
  refusing it, which `rdm()` now refuses. Changing `lba()` would touch
  a file a sibling lane may also touch and would change a shipped
  family's behavior, so it is recorded here instead.

### The whole suite, one process

`NOT_CRAN=true`: **1100 pass, 0 fail, 0 error, 1 skip**, over 16 files.

| file | passing |
|---|---|
| test-brms-parity.R | 13 |
| test-defects.R | 59 |
| test-density.R | 132 |
| test-family.R | 34 |
| test-gddm-family.R | 77 |
| test-gddm-gradients.R | 22 |
| test-gddm-recovery.R | 28 |
| test-gddm-solver.R | 94 |
| test-lba.R | 109 |
| test-message-uniqueness.R | 4 |
| test-moments.R | 29 |
| **test-rdm-gng.R (new)** | **149** |
| test-sampling.R | 94 |
| test-simulate-density.R | 69 |
| test-surface.R | 47 |
| test-variability.R | 140 |

The 926 recorded at 0.2.0 is not directly comparable to 1100, because
1100 was taken with `NOT_CRAN=true` and the CRAN-skipped tests
elsewhere in the suite then run too. This lane added 155 of the
difference: 149 in the new file and 6 in `test-simulate-density.R`,
which went from 63 to 69 (both measured).

CRAN mode (`NOT_CRAN` unset): **959 pass, 0 fail, 0 error, 18 skip**,
over 15 files. `test-simulate-density.R` carries `skip_on_cran()` at
the top of the file, so it contributes nothing at all there, which is
why the file count drops by one; `test-rdm-gng.R` contributes 138 of
its 149.

Neither number is a clean successor to the 926 recorded at 0.2.0.
`test-simulate-density.R` did not exist at 0.2.0 (the development
section of NEWS.md introduces it), and 926 was not annotated with the
`NOT_CRAN` setting it was taken under. What this lane added is
measured directly instead: 149 new passing tests under `NOT_CRAN=true`
in a new file, plus 6 more in `test-simulate-density.R`, which went
from 63 to 69 with both numbers taken on this worktree.

### R CMD check --as-cran

`_R_CHECK_CRAN_INCOMING_=false`, core installed in `rg-lib`, pandoc
from the RStudio quarto tools directory on PATH, `NOT_CRAN=true` so
the check ran the heavy tests rather than skipping them.

```
* checking tests ...
  Running 'testthat.R' [31m] OK
* checking re-building of vignette outputs ... [135s] OK
Status: OK
```

**No errors, no warnings, no notes.** The vignette rebuild is the part
worth naming: it runs every fitted example in both new sections and
the live comparisons against statmod and WienR.

# Punch round, 2026-09-05

Against `dev/review-rdm-gng.md` (verdict GO WITH FIXES). Every number
below was re-measured on this worktree rather than copied from the
review; where the review and this lane disagreed, the review was right
every time.

## Item 4 (coordinator): the core trap, stated precisely

My first write-up said "core falls back to the first primary dpar",
which is right about `predict()` and loose about `fitted()`. The
mechanism is three sites, and the reviewer located all of them:

1. `frmtmb` `R/predict.R:639-642`, `mean_is_mu()` is the root:

   ```r
   mean_is_mu <- function(fam) {
     is.null(fam[["post"]]$mean_fn) ||
       identical(body(fam[["post"]]$mean_fn), quote(dpars[["mu"]]))
   }
   ```

   **A family that declares NO mean is read as one that declares "the
   mean is mu".** Missing is taken for agreement.
2. `R/predict.R:1168-1169`: because `mean_is_mu()` is TRUE,
   `type = "response"` never routes to `predict_mean_response()` and
   falls through to
   `dpar %||% if ("mu" %in% names(...)) "mu" else primary_dpars[1]`.
   That is where `rdm()`'s 3.51 came from.
3. `R/families.R:437-441`, `response_mean()`: the `else` branch is
   `dpars[["mu"]]`, which is the drift for `wiener_gng()` and NULL for
   `rdm()`, hence the length-zero residuals.

`fitted()` at `R/predict.R:1704` refuses only when BOTH
`!"mu" %in% names(dp)` and `is.null(post$mean_fn)` hold, so **a family
is protected by NOT having a dpar called `mu`, which is backwards**.
`rdm()` was saved from `fitted()` by an accident of naming;
`wiener_gng()`, whose drift is called `mu`, was not.

**Core left alone**, as instructed; the one-line fix is the
consolidation's.

### The refusals do not depend on which way that guard falls

This is the part the coordinator asked me to make sure of, and it is
measured rather than argued. `.../scratchpad/rg-meantrap.R` patches the
FIXED `mean_is_mu()` into the loaded `frmtmb` namespace and re-runs
both families:

| call | core as shipped | core fixed |
|---|---|---|
| `fitted(rdm fit)` | REFUSED, `rdm: ...` | REFUSED, `rdm: ...` |
| `predict(rdm, "response")` | REFUSED, `rdm: ...` | REFUSED, `rdm: ...` |
| `residuals(rdm, "response")` | REFUSED, `rdm: ...` | REFUSED, `rdm: ...` |
| `fitted(gng fit)` | REFUSED, `wiener_gng: ...` | REFUSED, `wiener_gng: ...` |
| `predict(gng, "response")` | REFUSED, `wiener_gng: ...` | REFUSED, `wiener_gng: ...` |
| `residuals(gng, "response")` | REFUSED, `wiener_gng: ...` | REFUSED, `wiener_gng: ...` |

Identical in both columns. A declared `mean_fn` whose body is not
`quote(dpars[["mu"]])` makes `mean_is_mu()` FALSE under either
definition, so the request routes into `response_mean()` and meets the
family's own `stop()` either way. Pinned in the suite: the refusal is
asserted to start with the family's own name, and `rdm()` is asserted
to have no dpar called `mu` while `wiener_gng()` does, so the two
cover both sides of the core guard.

## Item 2: the WienR sentence was overstated

Reproduced exactly. At `t = 2.5, a = 4, v = 5, w = 0.75` against a
260-bit reference of 1.1883283058299342e-13, `1 - WienR::pWDM()` is
4.4239 percent wrong at the default and at 1e-8, 1e-10 and 1e-12, then
0.1262 percent at 1e-14 and 0.0328 percent at 1e-16.
`WienR::WienerCDF()` behaves identically. So "at every precision
setting from 1e-8 to 1e-16" is false at the top of the range.

Corrected in `dev/rdm-gng-findings.md` (the section above),
`R/wiener-gng.R` roxygen, `vignettes/ddm.Rmd` and `NEWS.md`. The
substantive point is stronger stated correctly, and I verified the
fact that carries it: `EMC2:::pDDM` and `EMC2:::dDDM` both default to
`precision = 0.005`, looser than every setting in that table, so EMC2
in practice sits at the 4.4 percent end.

## Item 3: `lba()` silently dropped `dec()`

Reproduced. On 300 rows, `rt | dec(two) + vint(choice) ~ 1` and
`rt | vint(choice) ~ 1` under `lba(3)`:

* largest absolute difference in fixed effects: **exactly 0**
* log likelihoods identical
* **no warning at any point**
* the fitted object carries both `dec` and `vint1` in its aterm values

I checked the mirror case the reviewer flags too, and it is real:
`wiener()` fits `rt | dec(upper) + vint(extra) ~ 1` with a log
likelihood difference of exactly 0, accepting a `vint()` term it cannot
use. So this is not an `lba()` defect so much as one instance of a
general one: **frmtmb has no seam for "this family does not accept
that addition term."** `required_aterms` is a conjunction of what the
density NEEDS, with no complementary allow-list.

**Fixed for the two race families only.** `ddm_refuse_dec()` in
`R/ddm-shared.R` is one shared check, called by `lba_check_response()`
in `R/lba.R` and `rdm_check_response()` in `R/rdm.R`. One shared
check rather than two copies for two reasons: the message-uniqueness
test wants one template per condition, and two sibling families with
the same dpars should not be able to drift apart on this again. The
`lba` and `rdm` `dec()` compat rows were both rewritten, because the
`lba` one asserted a refusal by the two-level coercion that did not
happen.

`wiener()`'s `vint()` case is left alone and recorded: it is a
different family, a different term, and the durable fix is a core
allow-list rather than a third hand-written check.

## Item 6: the dev script claimed a comparison it never made

`dev/rdm-gng-emc2-reference.R` printed "EMC2 cannot score 0 of 200
rows at all ... e.g. n/a" on every accumulator count. The detection was
not broken; the grid is deliberately chosen to keep EMC2 inside its own
regime, so the count is legitimately zero and the line should not have
printed at all. It is now printed only when there is such a row, with
the zero case saying so plainly. The saturation itself is real and is
reported by the survival section, which uses a grid built to provoke
it.

## Item 5: EMC2 in Suggests, kept but flagged

The reviewer is right on the facts, and I verified them: `^dev$` is in
`.Rbuildignore`, the built tarball contains **0 files under `dev/`**,
and every EMC2 mention in shipped material (roxygen, vignette prose,
test comments) is prose. Nothing executable in the package needs EMC2,
so `dependencies = TRUE` pulls EMC2 plus GPArotation, magic,
matrixcalc, corrplot, psych and lpSolve for nothing.

**Kept anyway, and this is a judgement call worth flagging rather than
burying.** The lane brief asked for "the DESCRIPTION bump to 0.3.0
with EMC2 in Suggests" in as many words. It became unnecessary only
because the same brief also said never to use `:::` in the package
tests, which is what pushed the EMC2 comparison out to `dev/`. A
DESCRIPTION field cannot carry a comment, so the honest options were to
drop it or to say so here.

Dropping it is a one-line change and the consolidation should feel free
to make it. What it costs to keep: six extra packages under
`dependencies = TRUE`. What it buys: `dev/rdm-gng-emc2-reference.R`
runs without a separate install step for anyone checking the EMC2
claims that the roxygen and the vignette make.

## Item 1: the blend constant was wrong, and the claim with it

The review's FINDING 1 is correct and I reproduced it independently.
Grid: 1200 rows, `t` in {0.05, 0.1, 0.2, 0.5, 1, 2, 2.5, 5, 10, 15},
`v` in {-2, -0.5, 0.5, 1, 2, 5}, `a` in {0.8, 1.4, 2.5, 4}, `w` in
{0.25, 0.45, 0.5, 0.75, 0.9}. Reference is the 260-bit small-time image
sum, cross-checked on 30 rows against the 260-bit LARGE-time route,
which agree with each other to between 1.3e-18 and 2.6e-16 over that
subset and to between 8e-64 and 7e-77 on the pinned rows.

Neither route alone is usable over the whole grid: the small-time
route reaches **12.1 relative** at its worst and the large-time route
8.2e-03 at its worst, each outside its own regime.

### The sweep

| u0 | us | max rel | rows > 1e-12 | > 1e-9 | > 1e-6 |
|---|---|---|---|---|---|
| 0.06 | 0.12 | 7.75e-05 | 14 | 5 | 1 |
| 0.04 | 0.12 | 9.00e-08 | 7 | 1 | 0 |
| 0.03 | 0.12 | 7.45e-10 | 3 | 0 | 0 |
| **0.02** | **0.12** | **2.17e-12** | **1** | **0** | **0** |
| 0.015 | 0.12 | 2.17e-12 | 2 | 0 | 0 |
| 0.01 | 0.12 | 4.57e-11 | 4 | 0 | 0 |
| 0.005 | 0.12 | 3.24e-06 | 29 | 4 | 1 |
| 0.02 | 0.06 | 2.17e-12 | 1 | 0 | 0 |
| 0.06 | 0.20 | 4.47e-02 | 20 | 7 | 3 |

Independent of the reviewer and matching to the digit: 7.75e-05 at
0.06 and 2.17e-12 at 0.02, with 14 rows over 1e-12 against the
reviewer's 15 (one row's worth of grid difference).

`ddm_cdf_u0` is now **0.02**, at `R/wiener-cdf.R:67`. `us` stays 0.12;
0.06 gives the same worst case and changing two constants where one
will do is not an improvement.

Note the bottom of the table: 0.005 is WORSE than 0.02, at 3.24e-06.
The two routes cross near `u = 0.025` and below the crossing the small
route is the accurate one, so an over-low centre hands over too early.
That is why this is a measured constant with a two-sided optimum
rather than "as low as possible".

### The mechanism, at the reviewer's worst row

`t = 2.5, v = 5, a = 4, w = 0.90`, true no-go probability
8.2014416716939391e-16:

| centre | small route | large route | blend |
|---|---|---|---|
| 0.06 | 12.1 relative | 6.73e-15 | **7.75e-05** |
| 0.02 | 12.1 relative | 6.73e-15 | **8.78e-13** |

The series were never the problem. The small-time route is not gently
wrong there, it is 12 times wrong, because `1 - F_upper` has nothing
left to subtract from; at `u = 0.156` the 0.06 weight still gave it a
share of about 3e-05, and a small share of a hopeless number is the
whole error. **A smooth blend has no safe side unless both branches
degrade gently outside their regime**, and this one has a branch that
collapses, so the weight must saturate before the collapse rather than
after it. That is the general lesson and it is now written at the top
of `R/wiener-cdf.R`.

### Why the lane's own grid could not see it

The suite's no-go grid had `w` in {0.3, 0.5, 0.7} and kept only rows
with `ref > 1e-6`. The corner needs a large drift AND a start point
biased toward the go boundary AND a wide boundary, all together. Both
halves are fixed:

* the WienR grid now runs `w` in {0.25, 0.3, 0.5, 0.7, 0.9} and is
  STRATIFIED rather than filtered, with a tight tolerance where WienR
  is trustworthy (`ref > 1e-6`, 1e-9 relative) and a part-in-a-thousand
  tolerance below it, which is all WienR has left there;
* a new test pins six 260-bit constants directly, including the
  reviewer's exact point, at 1e-11 relative, with no Rmpfr at test
  time. It also asserts the MECHANISM: on that row the large-time route
  is within 1e-13 and the small-time route is more than 1 relative, so
  the blend is only right because the weight has saturated.
* a second new test pins `ddm_cdf_u0` and `ddm_cdf_us` themselves and
  checks the weight saturates at the pinned corner and does NOT
  saturate at a genuinely small-time row, so a later edit has to
  re-measure rather than nudge.

### What the claim should have been

The lane wrote that the blended probability "holds 1.0e-12 relative"
against a 260-bit reference. On the grid it was measured on that was
true; as a general claim it was not, and the honest version is now in
the roxygen, the vignette and NEWS: **2.2e-12 over 1200 points**
spanning start points to 0.9, one row worse than 1e-12 and none worse
than 1e-9.

## Punch-round verification runs

Six named files, one process each, `NOT_CRAN=true`:

| file | before | after |
|---|---|---|
| `test-rdm-gng.R` | 149 | **169** |
| `test-lba.R` | 109 | 109 |
| `test-surface.R` | 47 | 47 |
| `test-message-uniqueness.R` | 4 | 4 |
| `test-defects.R` | 59 | 59 |
| `test-simulate-density.R` | 69 | 69 |

0 failures, 0 errors, 0 skips in every one.

The 20 new tests in `test-rdm-gng.R` are the punch round's: the six
260-bit constants and the mechanism assertions around them, the two
blend-constant pins, the widened and stratified WienR grid, the
`lba()` / `rdm()` `dec()` agreement block, and the four assertions that
the mean refusals do not depend on which way core's `mean_is_mu()`
falls.

`test-lba.R` is unchanged at 109 despite `lba()` gaining a refusal,
which is the point: the new check fires only on a term no existing
test supplies.

Whole suite, one process, `NOT_CRAN=true`: **1120 pass, 0 fail,
0 error, 1 skip**, over 16 files, up from 1100 by exactly the 20 above.

The EMC2 identity is unchanged by the blend constant, re-run after it:
go rows 6.11e-16, no-go rows 1.33e-15, total log likelihoods differing
by 1.71e-13, and the racing-diffusion race 6.4e-13 / 1.0e-12 / 1.1e-12
at two, three and four accumulators. Agreement with `EMC2:::pDDM`
where the no-go probability is above 0.01 improved from 3.34e-13 to
1.32e-13 as a side effect of the new centre.

`R CMD check --as-cran` re-run on the punch-round source, same
environment as before (`_R_CHECK_CRAN_INCOMING_=false`, core in
`rg-lib`, pandoc from the RStudio quarto tools, `NOT_CRAN=true`):
tests OK in 18m, vignette rebuild OK in 79s, **Status: OK**, zero
errors, zero warnings, zero notes.

## What was NOT done, and why

* **Core is untouched.** `mean_is_mu()` at `R/predict.R:639-642`,
  `response_mean()`'s `else` branch at `R/families.R:437-441` and
  `fitted()`'s guard at `R/predict.R:1704` are the consolidation's,
  per the coordinator. This lane only made sure its own refusals do
  not depend on which way that guard falls, which is measured above.
* **`wiener()` still accepts a `vint()` term it cannot use**, silently
  and with a log-likelihood difference of exactly zero. Reproduced,
  recorded, left alone: it is a different family and a different term,
  and the durable fix is a core `permitted_aterms` allow-list rather
  than a third hand-written check. The `lba()` hunk was taken because
  `lba()` and `rdm()` are siblings with identical dpars whose
  disagreement was the actively harmful case.
* **`EMC2` stays in `Suggests`**, against the review's tidying
  suggestion, because the lane brief asked for it in as many words.
  The cost is measured and stated above; dropping it is a one-line
  change for the consolidation.

## Observation, not this lane's doing: main is mid-merge with a conflict

Noticed on the final read-only status check at the end of the punch
round, and reported rather than touched.

`C:/Users/adf44/source/r/frmtmb` is at **d34d309 "Merge branch
'wt-ce'"** with an unfinished merge in progress: `MERGE_HEAD` is
**e18045e "car(type = \"esicar\"): the exactly constrained intrinsic
CAR brms fits"**, and `git status` shows `UU NEWS.md`, a CONFLICT in
the core NEWS.md, alongside staged `R/covstruct.R`, `R/frame.R`,
`R/objective.R`, `R/compat.R`, `tests/testthat/test-car-spde.R` and a
tranche of new `dev/es-*` files.

**Not this lane.** Everything staged belongs to the esicar / CAR lane;
nothing under `extensions/frmtmb.ddm/` appears in the list at all.
Main was clean when this lane checked it earlier in the same session,
and every git command this lane has ever run against that repository
has been read-only (`status`, `log`, `rev-parse`, `reflog`). The
reflog shows main moved twice more since the review was written
(`316a28b` merging wt-gddm-ref, then `d34d309` merging wt-ce), so a
consolidation is working through the lanes and is part-way into the
esicar one.

Left exactly as found. Whoever owns that merge needs to resolve
`NEWS.md`; this lane still rebases onto main cleanly, and its own
files are untouched by it.
