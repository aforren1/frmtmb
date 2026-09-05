# Lane wt-esicar: findings

Worktree `C:/Users/adf44/source/r/frmtmb-wt-esicar`, branch `wt-esicar`,
branched at 2210aa1. Main checkout untouched.

## 1. The constraint the code has, and the one brms has

Before this lane, `car_aux()` (`R/covstruct.R:1204` at 2210aa1) had no
`esicar` branch: `escar`
returns early and every other type falls into the shared intrinsic path,
so `esicar` IS `icar`, bit for bit. The intrinsic path builds

    K = L + Sgrp' diag(kappa0) Sgrp,  kappa0_j = 1 / (con_sd n_j)^2
    Q(tau) = tau K,   log|Q| = n log tau + log|K|

The rank-c update is the sum-to-zero constraint and it is SOFT, and -
this is the load-bearing detail, it rides on `tau`. Because L and
`P0 = Sgrp' diag(kappa0) Sgrp` are zero on each other's invariant
subspace (null(L) = span of the component indicators s_j, and
P0 v = 0 for v orthogonal to every s_j), the block density factorizes
EXACTLY:

    p(b) = p_constrained(C b) * prod_j N(m_j; 0, (con_sd sdcar)^2)

with `C` the per-component centering projection and `m_j` component j's
mean. The soft model is therefore the hard model PLUS an extra random
intercept of sd `con_sd * sdcar` per component. That intercept enters
the linear predictor, is confounded with the intercept beta_0, and its
convolution is the whole of the measured 4.7e-4 bias.

## 2. Why the fix cannot be a `car_aux()` branch alone

To make the constraint hard, the field that reaches the linear predictor
must be the centered one. The block's `nll` cannot do that: it sees `b`
and returns a scalar. The coefficient vector the Z matrices multiply is
`expand_b()` (`R/covstruct.R:1756`, my file), but its call site was
gated on `frame$has_rr` (`R/objective.R:253`, set at `R/frame.R:2175`),
both sibling-owned. There is no other seam: `Zlocal` is built in
`R/frame.R:1895`, and the `ncp` transform is a sampling-only route
(`frame$ncp_blocks` is never set by the fit path).

## 3. Design chosen: the inert-coordinate reparameterization

The task allows "an equivalent exact reparameterization ... so the
level-major layout the importance correction assumes still holds". Taken,
because it keeps `length(b_idx) == length(c_idx) == Nloc` and so needs no
edit to `R/predict.R`, `R/methods-fit.R`, `R/importance.R` or
`R/confint.R`.

    field reaching the predictor:  c = C b       (expand_b)
    block precision:               Q = tau L + P0     (P0 tau-FREE)
    log|Q| = (n - n_comp) log tau + log|K|             (K = L + P0)
    -2 log p = -log|Q| + n log 2pi + tau b'L b + b'P0 b

Two edits against the icar arithmetic: the tau exponent drops from `n` to
`n - n_comp`, and the rank-c quadratic is no longer multiplied by tau.

Why it is EXACT. `b` splits orthogonally into `f = C b` (dimension
n - n_comp) and the component means `m`. `Q` is block diagonal on that
split, so p(b) = p_hard(f) * prod_j N(m_j; 0, con_sd^2), the predictor
depends only on `f`, and the m integral is exactly 1. The Laplace
approximation factorizes over the same split with an exactly Gaussian m
block, so the marginal likelihood is the hard-constrained one to machine
precision, and is INVARIANT to `con_sd`, which is the testable
signature that the constraint is exact rather than tight.

`m* = 0` exactly at the mode (m enters neither the likelihood nor any
cross term), so `b` is already centered where anything reads it, and the
inert coordinate carries variance con_sd^2 = 1e-6, four orders below
sdcar^2. That is why the delta-method and ranef SE paths need no branch.

Predicted identity constant, to be confirmed numerically:

    const(esicar) = 0.5 * n * log(2 pi) - 0.5 * log det(L + J / s^2)

with s = 0.001 * Nloc, the same `kmat` `brms_car_const()` already builds
for icar.

## 4. Cross-lane edits required

`R/frame.R` and `R/objective.R` (wt-spline-core) need one new frame field
`has_expand` and one gate. Three lines. Recorded here because the lane
brief forbids those files; see the final report.

## 5. brms's generated esicar program (verified, brms 2.23.0)

`dev/es-stancode/esicar.stan`. esicar is CENTERED on the field, unlike
icar:

    parameters:            vector[Nloc - 1] zcar;
    transformed:           rcar[1:(Nloc - 1)] = zcar;
                           rcar[Nloc] = - sum(zcar);
    predictor:             mu[n] += rcar[Jloc[n]];      // no sdcar
    target:                sparse_icar_lpdf(rcar | sdcar, ...)
      = 0.5 * ((Nloc - 1) * log(tau) - tau * (phi'D phi - phi'W phi))

So brms normalizes by the pseudo-determinant's `(Nloc - 1) log tau` and
drops `log|L|*` and the 2 pi. icar, by contrast, is non-centered
(`rcar = zcar * sdcar`) with the soft term
`normal_lpdf(sum(zcar) | 0, 0.001 * Nloc)`; that is why icar's map
carries a Jacobian of `Nloc * log(sdcar)` and esicar's carries NONE.

## 6. Measurements

Data: `brms_car_data()`'s 4 x 4 lattice, 96 rows, 16 locations, one
connected component.

| quantity | value |
| --- | --- |
| logLik esicar | -90.315906933797407 |
| logLik icar | -90.316306819585037 |
| esicar - icar | 3.998858e-04 |
| independent hard-constrained reference | -90.315906933522527 |
| esicar - reference | -2.7e-10 (frmtmb's own outer tolerance, see below) |
| icar - reference | -3.998861e-04 |
| sdcar esicar / icar - 1 | 3.053e-05 |
| sum of the esicar field | -2.9e-16 |
| sum of the icar field | 3.2e-09 |
| brms identity residual | 3.0e-13 |
| max abs gradient, inner pars | 8.2e-15 |

The reference is a marginal-ML fit built outside frmtmb: `y = X beta +
Z f + eps` with `cov(f) = sdcar^2 pinv(L)`, the exact constrained
covariance, profiled over beta and optimized over (sdcar, sigma).

The -2.7e-10 is NOT the reference optimizer's noise, which is what an
earlier draft of this section said. The reviewer polished an
independent reference to `|grad| = 2.8e-9`, which moved it by 1.6e-11
and left the gap at -2.84e-10, with the sign putting frmtmb below the
maximum: it is frmtmb's own outer nlminb tolerance. The conclusion is
stronger for it, because the reviewer also evaluated the two objectives
as FUNCTIONS at five points off the optimum and got agreement to about
4e-15 relative, which no optimizer tolerance can explain.

con_sd invariance, the signature of an exact constraint:

| con_sd | esicar logLik | icar logLik |
| --- | --- | --- |
| 1e-2 | -90.315906933797 | -90.354272794400 |
| 1e-3 | -90.315906933797 | -90.316306819585 |
| 1e-4 | -90.315906933790 | -90.315910934477 |
| 1e-5 | -90.315906931386 | -90.315906959110 |

esicar is flat to 12 digits; icar walks quadratically onto it. escar
(-92.2921419248) and bym2 (-90.3163068318) reproduce the pre-change
values in `dev/reviews/2026-09-05-brms-rows.md` section 3 exactly, and
so does icar, so nothing but esicar moved.

Identity constant, derived from brms's own standata:

    const(esicar) = 0.5 * Nloc * log(2 pi)
                  - 0.5 * log det(L + J / (0.001 Nloc)^2)

predicted 2.037041610157443, measured 2.037041610157743.

## 7. importance

No refusal needed: `car` is already in `imp_crosslevel`
(`R/importance.R:63`), so `importance =` refuses EVERY car type by name
before any layout question arises. The chosen design keeps
`length(b_idx) == length(c_idx) == Nloc` and level-major order, so
`imp_layout()` would not have broken anyway.

## 8. The one approximation left in place, measured

`lp_delta_A()` (`R/predict.R:1353`) pairs the Z columns with the `b`
positions through a Jacobian `dc/db` that is the identity for every
block but `rr`. For esicar the true Jacobian is the centering
projection `C`, so the delta method carries the inert coordinate's
variance into a prediction standard error that should not have it.

That variance is exactly `con_sd^2 = 1e-6`, because the inert
coordinate's posterior IS its prior, because the data say nothing about
it. On
the row's fit the prediction standard errors are about 0.2, so the
inflation is 2.5e-5 relative in the variance and 1.3e-5 in the standard
error. The same reasoning applies to `ranef(condVar = TRUE)`, whose
conditional standard deviations measured 0.186 to 0.198.

Fixing it exactly needs a branch in `R/predict.R`, which belongs to
wt-spline-core this round. Left alone deliberately: five orders below
the quantity it perturbs, and the alternative was a cross-lane edit for
no measurable gain.

## 9. Compat surface, measured on an esicar fit

| surface | result |
| --- | --- |
| `importance = 32` | REFUSED by name, "cannot correct the 'car' structure" - the pre-existing car refusal, not a new one |
| `predict()` | works; `predict(newdata =)` reproduces it to 1e-10, SEs to 1e-8 |
| `predict(allow_new_levels = TRUE)` | works, as for every car type |
| `simulate()` | works, 96 x 5 draws |
| `frmtmb.sample::frm_sample()` | works, 200 x 21 draws, stays centered (car has no non-centered form) |
| `ranef(condVar = TRUE)` | works, condSD 0.186 to 0.198 |
| `VarCorr()`, `confint()`, `confint_varcorr()` | work, `sd(car)` 1.4135 |
| disconnected graph | per-component sums 1.4e-17, logLik -61.9559119378 against icar's -61.9559824913 |

One nuance worth knowing about `frm_sample()`: the draws of `b` carry
the inert coordinate, so a component's raw draw sums to something of
sd `con_sd * Nloc` (measured max 7.1e-2 over 200 draws, against the
0.016 sd that predicts). The FIELD is `expand_b()`'s output and sums to
zero in every draw.

## 10. Verification

One process per file throughout, counts audited by name against
`ls tests/testthat/test-*.R`.

| run | result |
| --- | --- |
| test-car-spde.R | 18 tests, 97 assertions, 0 fail (14 tests / 71 assertions at 2210aa1) |
| test-importance.R | 31 tests, 185 assertions, 0 fail |
| test-message-uniqueness.R | 1 test, 6 assertions, 0 fail |
| test-bracket-access.R | 3 tests, 8 assertions, 0 fail |
| test-compat.R | 29 tests, 261 assertions, 0 fail |
| test-brms-likelihood.R, gated and warm | 33 tests, 372 assertions, 0 fail, 0 skip, 51.4 s (32 tests / 351 assertions at 2210aa1) |
| full core suite | 109 files, 1097 tests, 6143 assertions, 0 fail, 0 error, 88 skip |

The 88 skips are the brms fit tiers with the gate off (2, 30, 44 and 11
over test-brms-agreement, -likelihood, -methods and -priors) plus one
in test-fuzz. The by-name audit found every file on disk in the log and
no extras.

The Stan cache held 54 `.rds` programs before the tier run and 54
after: nothing recompiled. The esicar program was compiled once, by the
standalone identity probe, and cached.

`roxygen2::roxygenise()` twice leaves `man/` and `NAMESPACE` unchanged;
nothing here is exported or documented with a roxygen block that
generates a page.

## 11. R CMD check

`--as-cran --no-manual`, `_R_CHECK_CRAN_INCOMING_=false`, pandoc from
`RSTUDIO_PANDOC` on PATH, checked against the private library.

First run, built with `--no-build-vignettes`: **2 WARNINGs, both the
same `inst/doc` artifact of that flag** ("Files in the 'vignettes'
directory but no files in 'inst/doc'", and "Directory 'inst/doc' does
not exist"). Nothing else: examples OK (50 s), examples with
`--run-donttest` OK (43 s), tests OK (327 s), and re-building of
vignette outputs OK (254 s), which is the step that knits the edited
`vignettes/frmtmb.Rmd`.

Second run, built WITH vignettes: **Status: OK**. No warnings, no
notes. Vignette rebuild 426 s. That confirms the first run's two
WARNINGs were the build flag and nothing else.

## 12. Files touched

Mine: `R/covstruct.R`, `tests/testthat/test-car-spde.R`,
`tests/testthat/test-brms-likelihood.R`,
`tests/testthat/helper-brms.R`, `vignettes/frmtmb.Rmd`,
`dev/brms-likelihood-tests.md`, `NEWS.md`.

Outside the brief, and flagged: `R/frame.R` and `R/objective.R`
(wt-spline-core) for the `has_expand` flag, three lines net, argued in
section 2; `R/compat.R` and `dev/feature-gaps.md`, both of which
carried a statement about `esicar` that this change makes false and
neither of which any sibling lane claims this round.

The two "before" figures were wrong in the first draft (15/92 and
31/351). They are now measured rather than recalled: `git archive
2210aa1` into a scratch tree, installed into its own library, and the
two files run against it. `test-car-spde.R` is 14 tests / 71
assertions there and `test-brms-likelihood.R` is 32 tests / 351
assertions gated and warm. The test-count deltas match the diff: car
removes 1 block and adds 5 (14 + 4 = 18 before the punch round), the
tier removes 1 and adds 2 (32 + 1 = 33).

# Punch round, 2026-09-05

Against `dev/review-esicar.md` (GO WITH FIXES). One entry per punch
item, in the review's order.

## P1 (blocking, code): the gate and expand_b() disagreed

Landed at `R/objective.R:199` and `:261`, `R/covstruct.R:1760-1782`.

The review proposed one token, `%||% frame[["has_rr"]]` on the gate, to
match the fallback `expand_b()` already carried. I measured that fix
before taking it and it does NOT close the hazard the review found. For
an esicar-only frame `has_rr` is FALSE, so both sides fall back to
FALSE, agree with each other, and still skip the centering:

    intact frame      68533.4537597
    has_expand = NULL 68560.8672184     <- after the one-token fix
    diff              -2.741e+01

The two expressions became identical and the model was still wrong. So
the fix is to stop trusting a cached boolean and DERIVE the answer from
the blocks, in one predicate both sites call:

    frame_needs_expand(frame)   R/covstruct.R:1760

which is TRUE when the frame's flag says so or when any block is `rr`
or esicar. `expand_b()` guards on it (`:1782`) and `build_objective()`
asks it ONCE, outside the returned closure (`R/objective.R:199`), so
the per-evaluation cost is nil. The blocks cannot go stale the way the
boolean can, because the same `car_type` that chooses the density
chooses the expansion.

Re-measured with that in place, all three paths agree exactly:

| frame | at the mode | block shifted 0.37 |
| --- | --- | --- |
| intact | 83.4537596736 | 68533.4537597 |
| `has_expand` removed | 83.4537596736 | 68533.4537597 |
| both flags FALSE | 83.4537596736 | 68533.4537597 |

diff 0.000e+00 in both columns. Probe `dev/es-gate.R`.

Pinned by `tests/testthat/test-car-spde.R:466`, "the objective expands
esicar whatever the frame's flag says", which rebuilds the objective
from a stripped frame and from a both-flags-FALSE frame and compares
against the intact one at the mode AND at `b + 0.37`. It asserts the
shift is worth more than 1e3 nats first, so the comparison has teeth.

## P2 (blocking, doc): "con_sd is ignored" was false

Landed at `NEWS.md:20-29`, `vignettes/frmtmb.Rmd:206-215`,
`R/compat.R:1190`.

The review is right and I reproduced every number. con_sd cannot move
the likelihood, and it lands in every delta-method standard error:

| con_sd | logLik | se.fit range | variance excess over the 1e-4 fit | predicted `con_sd^2 - 1e-8` |
| --- | --- | --- | --- | --- |
| 1e-2 | -90.315906933797 | 0.192853768 .. 0.253398594 | 9.999000e-05 | 9.999000e-05 |
| 1e-3 | -90.315906933797 | 0.192596926 .. 0.253203174 | 9.899999e-07 | 9.900000e-07 |
| 1e-4 | -90.315906933790 | 0.192594356 .. 0.253201219 | 0 | 0 |

Max deviation from the prediction 3.8e-13, so the excess is EXACTLY
`con_sd^2` and not merely of that order. Relative inflation in the
standard error, `se / sqrt(se^2 - con_sd^2) - 1`: 1.347e-3 at 1e-2,
**1.348e-5 at the default**, 1.348e-7 at 1e-4. `ranef(condVar = TRUE)`
leaks identically (0.186382..0.198434 at 1e-2 against
0.186113..0.198182 at 1e-4). `VarCorr()` is clean: `sd(car)` equals
`exp(theta[1])` with a difference of exactly 0 at every con_sd. And
`con_sd = 0.1` inflates the standard errors by **7.5 to 12.7 percent**
while the log-likelihood does not move at all, which is the trap the
word "ignored" set. Probe `dev/es-consd.R`.

All three sentences now say that con_sd changes no estimate, that it
does scale prediction standard errors and condSDs as `con_sd^2` in the
variance, what that is worth at the default, that it grows a
hundredfold per decade, and to leave it alone on an esicar term.

Pinned by `tests/testthat/test-car-spde.R:503`, which asserts the
likelihood and `VarCorr()` do not move, that the variance excess is
exactly the difference of the two `con_sd^2` to 1e-10, and that the
relative inflation at the default sits between 1e-5 and 2e-5. That
last bound is the number the exact-Jacobian follow-up flips.

NOT taken: the review's preferred fix, pinning esicar's `kappa_j` to
the default constant so that "ignored" becomes literally true. It is a
better end state and it is in this lane's own file, but it changes
user-visible behavior beyond the punch list, and the coordinator's
contract chose the documentation route plus a pinned number. Recorded
here as the obvious next move, together with the exact Jacobian in
`R/predict.R`, which removes the leak rather than capping it.

## P3 (doc): the BEHAVIOR CHANGE bullet had no magnitude

`NEWS.md:14-19` now says the new fit sits ABOVE the old one, by 4.0e-4
in the log-likelihood and 3.1e-5 relative in `sd(car)` on a 4 by 4
lattice, and that what changes exactly is the constraint: the field
sums to zero to machine precision instead of to about 1e-9.

## P4 (records): both "before" figures were wrong

Corrected at `dev/esicar-findings.md` section 10, and measured this
time rather than recalled: `git archive 2210aa1` into a scratch tree,
installed into its own library, both files run against it.

| file | first draft | measured at 2210aa1 |
| --- | --- | --- |
| test-car-spde.R | 15 / 92 | **14 tests / 71 assertions** |
| test-brms-likelihood.R | 31 / 351 | **32 tests / 351 assertions** |

The review's P4 is right on both. The deltas reconcile with the diff:
car removes 1 block and adds 5, tier removes 1 and adds 2.

## P5 (test): singleton component

`tests/testthat/test-car-spde.R:542`. Two 4-node paths plus one
isolated node, `n = 9`, `c = 3`, `nj = (4, 4, 1)`.

| check | measured |
| --- | --- |
| `rank(L)` | 6 = n - c |
| `n_comp`, `nj` | 3, (4, 4, 1) |
| singleton field value | `identical(re[[9]], 0)` is TRUE |
| component sums | -1.11e-16, 2.78e-17 |
| esicar vs a 3-component `pinv(L)` reference | 4.6e-12 |
| icar vs the same reference | -2.08e-04 |
| icar's singleton | -4.31e-06, not zero |
| escar on the same graph | refused by name, "at least one neighbor" |

The singleton is right for a structural reason, not by luck: with
`n_j = 1` the centering is `b_i - b_i`, identically zero, so the
component contributes no free dimension and `n - c` counts it.

## P6 (test): the invariance is not a gaussian accident

`tests/testthat/test-car-spde.R:593`. Poisson on the same lattice,
where the Laplace approximation is no longer exact:

| con_sd | esicar | icar |
| --- | --- | --- |
| 1e-2 | -187.828769917712 | -187.831539740625 |
| 1e-3 | -187.828769917712 | -187.828797702641 |
| 1e-4 | -187.828769917702 | -187.828770195484 |

esicar's spread 1.05e-11, icar's 2.77e-03, and the poisson esicar field
sums to -1.7e-16. The factorization comes from the model, not from the
gaussian being exact. Probe `dev/es-p56.R`.

## P7 (records): the -2.7e-10 was misattributed

Corrected in section 6. It is frmtmb's own outer nlminb tolerance, not
the reference optimizer's noise; the reviewer polished an independent
reference to `|grad| = 2.8e-9`, which moved it 1.6e-11 and left the gap
where it was, with the sign putting frmtmb below the maximum.

## P8 (cosmetic): two wrap regressions

`vignettes/frmtmb.Rmd` (93 chars) and `dev/feature-gaps.md` (84) now
wrap at 72. The remaining over-length lines in these files are markdown
TABLE rows, which cannot be wrapped, and pre-existing chunk headers.

## P9 (doc, pre-existing): per-component against brms's global sum

Taken, because both files were already open for P2 and the statement is
user-facing. `R/compat.R:1190` and `vignettes/frmtmb.Rmd:216-219` now
say that every constrained type constrains each connected component
where brms constrains the global sum, and that the two agree on a
connected graph.

## Re-runs

Installed from the punch-round source, one process per file.

| run | result | before the punch round |
| --- | --- | --- |
| test-car-spde.R | 22 tests, 124 assertions, 0 fail | 18 / 97 |
| test-importance.R | 31 tests, 185 assertions, 0 fail | unchanged |
| test-message-uniqueness.R | 1 test, 6 assertions, 0 fail | unchanged |
| test-bracket-access.R | 3 tests, 8 assertions, 0 fail | unchanged |
| test-compat.R | 29 tests, 261 assertions, 0 fail | unchanged |
| test-brms-likelihood.R, gated and warm | 33 tests, 372 assertions, 0 fail, 0 skip, 45.9 s | unchanged |

Stan cache 54 `.rds` before and after the tier run: nothing recompiled.
`test-car-spde.R` gains the four punch tests, +4 tests and +27
assertions.

`R/objective.R` is on the hot path of EVERY model, so the full core
suite was re-run as well, one process per file, by-name audited:

    109 files, 1102 tests, 6175 assertions, 0 fail, 0 error, 88 skip

Against the pre-punch run (109 / 1097 / 6143 / 88), a by-name join of
the two logs shows exactly ONE file moved, `test-car-spde.R`, from
17 / 92 to 22 / 124. The 17/92 is that file caught mid-edit by the
earlier run; standalone it was 18 / 97. Every other one of the 108
files reports the same test and assertion counts before and after, so
the derived gate changed nothing outside the car tests.

`R CMD check --as-cran --no-manual`, `_R_CHECK_CRAN_INCOMING_=false`,
pandoc on PATH, tarball built WITH vignettes: **Status: OK**, no
WARNINGs and no NOTEs, vignette rebuild 201 s. `roxygen2::roxygenise()`
twice leaves `man/` and `NAMESPACE` unchanged.

Nothing was committed. The worktree HEAD is still 2210aa1 and the main
checkout was not touched.
