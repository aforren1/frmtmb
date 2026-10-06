# Lane gpby: `gp(x, by = f)` and exact `gp()` prediction at new positions

Branch `wt-gpby`, worktree `frmtmb-wt-gpby`, base 9e902909 (frmtmb
0.67.0). Library `wt-gpby-lib`; the "before" arm is `rellib-r5`. Nothing
is committed.

Every number below comes from a script in `dev/gpby-*.R` and the log
beside it (`dev/gpby-*.txt`, gitignored, local). Counts were pasted from
the logs, not typed.

## Punch round 1 (review `dev/reviews/2026-10-05-gpby.md`)

One entry per item. Scripts are `dev/gpby-p1-*.R`, logs beside them.
The reviewer's oracle stands as the record for the band itself: with
the hyperparameters held at the truth, 200 seeds, the simultaneous band
covers 192/200 on the lane against 188/200 on base, and its `"Sigma"`
equals the closed-form posterior covariance to 6.6e-7 relative.

### B1. `frm_curve_deriv(order = 2)` differenced the extra covariance

**Fixed by building the derivative from its sources.** New export
`frm_extra_cov_deriv()` (`R/extra-deriv.R`) returns the covariance of
the `order`-th derivative of the part that is not coefficient
uncertainty. An exact `gp()`'s kriging residual takes the squared
exponential kernel's closed form (`d^o k = (-1)^o He_o(u) / l^o k`,
`d^o d^o' k = (-1)^o He_2o(u) / l^2o k`, probabilists' Hermite
polynomials), minus `g K^-1 g'`. An unseen level differences the rows
`M` of its draw's design with the stencil and forms `(L M) S (L M)'`,
so no covariance is ever differenced. frmtmb.spline's
`frm_curve_deriv()` and `frm_curve_feature()` take it from there
(`sp_extra_deriv()`, `R/curve-cov.R`); the stencil matrix `L` is gone
from `R/curve-deriv.R`. Order 1 goes through the same function, and
the feature code's value and peak standard errors read the extra
variance at the root rows and `frm_extra_cov_deriv(order = 1)`.

The review's constructions on the final build:

| | before | after | reference |
|---|---|---|---|
| exact gp past data, order-2 `.se` at 3.0, 4.2, 4.8, 5.5 | 0, 0, 0, 0.875 | 0.284941, 0.691987, 0.830147, 0.958588 | eps 6e-3: 0.284936, 0.691976, 0.830138, 0.958580 |
| same, simultaneous crit on [3.5, 5.5] | 19.34 | 2.4994 | |
| same, order-1 `.se` | 0.115711, 0.454307, 0.615145, 0.650162 | 0.115697, 0.454314, 0.615238, 0.650127 | |
| `y ~ t + (1 + t \| g)` unseen level, order-2 `.se` | 0 to 1.495 | 0 to 0 | 0 |
| same, order-1 `.se` / closed form | [0.999849, 1.000159] | [1.000000, 1.000000] | 1 |

Logs `dev/gpby-p1-deriv2b-lane.txt` (the reviewer's
`dev/gpby-rev-deriv2b.R`) and `dev/gpby-p1-deriv2-lane.txt`
(`dev/gpby-rev-deriv2.R`; its GP half stencils `extra_cov` itself in
the reviewer's own code and stops, since `extra_cov` is no longer the
matrix it indexed). The order-1 `.se` moved by up to 1.5e-4 relative:
the order-1 stencil on the extra covariance was itself off by up to
0.4 percent of that part at the default step (the review's step table,
`dev/gpby-rev-deriv2.txt`), and the closed form replaced it.

Tests (`extensions/frmtmb.spline/tests/testthat/test-deriv.R`, two new
`test_that` blocks): the review's gp construction against a hand
closed form of `v1`, `v2` from the fitted theta (`.se^2 >= v2`, stable
to 1e-4 under a 4 times larger step, a 21-point simultaneous band with
`.se > 0` and crit below twice the pointwise one, and
`frm_extra_cov_deriv()`'s diagonal equal to `v1`, `v2` to 1e-8 of the
prior derivative variance, `s2 / l^2` and `3 s2 / l^4`: inside the
data `v1` is 1.7e-6 and both formulas lose it to the cancellation
against the prior variance, measured 8.4e-12 absolute,
`dev/gpby-p1-dbg1.R`); and the unseen-level construction (order 1
against `vcov(t, t) + VarCorr` to 1e-6, order 2 below 1e-6 of the
order-1 `.se`, the extra part below 1e-12 of its square). On the
pre-fix build: pass=39 fail=5 error=2
(`dev/gpby-p1-deriv-prefix.txt`); now pass=48 fail=0.

### m1. The fit-route chain that sticks is not a flat-prior chain

Both fixtures of `test-gp-by-draws.R` sample from the formula route
now, `frm_sample(bf(...), family = gaussian(), data = d, chains = 1,
iter = 300, refresh = 0, seed = 4)`, and the fixture's fit is `ds$fit`.

**Whether the fit route's sticking is a defect: yes, filed, and the
lane's first explanation was wrong.** Since frmtmb.sample 0.43.0 the
fit route applies brms's defaults too; on the review's construction
(`y ~ gp(x)`, 60 points, data seed 5) `prior_summary()` is the same
four rows on both routes, `student_t(3, 0.5, 2.5)` Intercept,
`student_t(3, 0, 2.5)` sigma and sdgp, `inv_gamma(1.804861, 0.108106)`
lscale (`dev/gpby-p1-m1.txt`). The fit route at `chains = 1,
iter = 600, seed = 4` has step size NaN, acceptance 0, 300 of 300
transitions divergent and every draw at one point, on rellib-r5 and on
the lane alike; the formula route on the same data and seed samples,
acceptance 0.97 (lane) and 0.92 (base), `dev/gpby-p1-m1b.R`,
`dev/gpby-p1-m1b-{lane,base}.txt`. It is pre-existing and not gp-by;
it is in `dev/test-backlog.md`, Open - medium, with the repro.

### m2. The kriging draw and the nugget

The draw now uses exactly the law the fit gives the latent values: the
smooth conditional covariance plus `gp_nugget * sd^2` as white noise
per DISTINCT position (`gp_krig_factor()`, `R/predict.R`). Two rows at
one position draw one value. The nugget value is unchanged at `1e-6`.

#### Nugget decision (decided by the user, 2026-10-06: keep 1e-6)

`dev/gpby-p1-nugget.R <nugget>` fits five designs with
`assignInNamespace("gp_nugget", ...)` before `frm()`;
`dev/gpby-p1-nugget-brms.R <nugget>` reruns `dev/gpby-brms-epred2.R`
(1500 draws of brms's own `y ~ gp(x, by = f)`) with that nugget. brms
2.23.0 jitters `K` by `1e-12` in the fit and adds `1e-8` to the
predictive variance; `1e-8` is measured as the middle value.

Fits (logs `dev/gpby-p1-nugget-{1e-6,1e-8,1e-12}.txt`):

| design | 1e-6 | 1e-8 | 1e-12 |
|---|---|---|---|
| `gp(x)`, 60 rows | conv 0, pdHess, logLik -63.511402, kappa(K) 1.4e7 | conv 0, pdHess, -63.511402, 1.4e9 | false convergence, max grad 8.1e-2, -63.511127, 1.4e13 |
| `gp(x, by = f)` | conv 0, pdHess, -30.762849 | conv 0, pdHess, -30.762862 | false convergence, no pdHess, -30.762624 |
| `gp(x1, x2)`, 2-D | conv 0, pdHess, -55.840858 | conv 0, pdHess, -55.840862 | conv 0, pdHess, -55.840862 |
| `gp(x)`, 200 rows 0.01 apart | conv 0, pdHess, -59.806187, kappa 9.9e7 | conv 0, pdHess, -59.806188, 9.9e9 | false convergence, no pdHess, -59.808650, 1.0e14 |
| `gp(x)`, 10 position pairs 1e-9 apart | conv 0, pdHess, -20.543721 | conv 0, pdHess, -20.543721 | false convergence, no pdHess, -20.542947 |

`chol(K)` succeeded at every fitted optimum, all three values. On
`gp(x)`, 60 rows, `se.fit` at 2.55, 6.5, 7.5 is 0.172941, 0.414199,
0.511395 at `1e-6`; 0.172941, 0.414134, 0.511420 at `1e-8`; and
0.173453, 0.363052, 1.889262 at `1e-12`, the last from the
unconverged fit. Two rows `1e-9` apart take draws whose difference has
sd 5.7e-4, 5.7e-5 and 5.8e-7 (the nugget's `sqrt(2 nugget) sd`): at
`1e-6` that is 0.27 percent of the field's conditional sd past the data
and larger than the conditional sd inside it (4.6e-4).

Against brms, at the same 1500 draws, per row (x, level): the
sampled field's sd across draws, frmtmb's over brms's, kriging draw on:

| row | 1e-6 | 1e-8 | 1e-12 |
|---|---|---|---|
| 2.55 a | 1.0002 | 1.0000 | 1.0000 |
| 6.50 a | 0.9946 | 0.9952 | 0.9950 |
| 7.50 a | 0.9784 | 0.9857 | 1.0161 |
| 2.55 c | 1.0011 | 1.0011 | 1.0001 |
| 6.50 c | 0.9817 | 0.9832 | 1.0185 |
| 7.50 c | 0.9829 | 0.9901 | 1.0179 |

and the per-draw conditional sd of the field given brms's own draw,
medians in units of the row's posterior sd (brms, then frmtmb):

| row | brms | 1e-6 | 1e-8 | 1e-12 |
|---|---|---|---|---|
| 6.50 a | 6.28e-3 | 7.37e-2 | 2.87e-2 | 6.38e-3 |
| 7.50 a | 1.07e-1 | 3.63e-1 | 2.35e-1 | 1.09e-1 |
| 6.50 c | 5.54e-4 | 2.21e-2 | 5.73e-3 | 4.68e-4 |
| 7.50 c | 1.30e-2 | 1.38e-1 | 6.23e-2 | 1.22e-2 |

Reading: per draw, only `1e-12` reproduces brms's conditional law
(the nugget conditions on noisy latent values and widens it, 3 to 40
times at `1e-6`); across draws all three are within 2.2 percent of
brms's spread, which is what an interval is made of. `1e-12` fails to
converge on four of five maximum-likelihood fits and gives an
unconverged `se.fit` of 1.89 against 0.51; `1e-8` converges on all
five, agrees with `1e-6` on logLik to 1.3e-5 and on `se.fit` to 5e-5
relative, and halves the per-draw gap to brms. The choice is the
user's.

**Decided 2026-10-06: the user keeps `1e-6`**, as the review
recommended (its r1, 180 fits: `1e-6` converges on 45 of 45 per
design set, `1e-8` warns on 6 of 45 and `1e-12` converges on none).
Recorded at the 0.68.0 release (`dev/round-20261005.md`) and in
`dev/round-handoff.md`. If per-draw parity with brms is wanted later,
the place for it is a sampled fit with brms's `1e-12` jitter, not a
maximum-likelihood fit.

### m3. Rows at one unseen position

Fixed in code: the draw is made once per distinct position and copied
(`gp_krig_factor()`'s `idx`). Rows 2 and 3 of the review's construction
(`dev/gpby-p1-krigdraw.R`, 300 draws): max |diff| 0.000e+00,
`identical()` TRUE, sd of the row 0.454. It was 1.498e-08.

### m4. The cost of the kriging draw

The per-draw eigen decomposition was replaced first by a dense
Cholesky, then by a pivoted Cholesky that builds only the kernel
columns it pivots on and stops once every remaining conditional
variance is below `gp_krig_tol = 1e-12` of the prior variance; what
remains joins the white part, so each row's variance is still
`extra_var`. The covariance, the draw and `extra_var` now share one
arithmetic, `P = Ks R^-1` with `R' R = K`, stored by `pred_design()`
in `krig$P` and its row sums in `krig$rs`; with two formulas the
pivots divided one's rounding by the other's, 5.0e-12 against
5.2e-13 now. Accuracy of the factor against `gp_krig_cov()`
(`dev/gpby-p1-krigfactor.R`):

    past     rows 300 | rank  10 | max |C - S| / sd^2 5.194e-13 | diag max rel 1.010e-11 | 0.80 ms per factor
    inside   rows 300 | rank  17 | max |C - S| / sd^2 3.785e-14 | diag max rel 1.053e-09 | 0.80 ms per factor
    neardup  rows   6 | rank   3 | max |C - S| / sd^2 2.040e-16 | diag max rel 6.675e-11 | 0.00 ms per factor

`posterior_epred()` at 300 unseen rows over 300 draws, CPU seconds,
three runs each, with a control of 300 rows at observed positions
(`dev/gpby-p1-krigdraw.R lane|lane dense|base`):

| arm | 300 unseen rows, CPU s | median | control, median |
|---|---|---|---|
| lane, pivoted factor | 0.75 0.78 0.88 | 0.78 | 0.18 |
| lane, dense Cholesky of `gp_krig_cov()` | 4.89 5.83 5.76 | 5.76 | 0.24 |
| base (conditional mean, no draw) | 0.59 0.73 0.63 | 0.63 | 0.38 |

The review measured 7.85 s against 0.47 s for the eigen version. The
machine is shared with other lanes, so CPU time is reported beside the
wall time in the logs and only same-session numbers are compared; an
earlier session put the dense draw at 3.44 s against 0.71 s pivoted,
before `gp_krig_cov()` was rebuilt for memory (m5). Timed one at a
time, 300 times, the steps of a dense draw cost 2.31 s, of which
`exp()` of the `n x n` kernel is 1.03 s, the Cholesky 0.38 s and the
kriging product 0.34 s (`dev/gpby-p1-krigparts.txt`); the pivoted
factor does none of the three at `n x n`.

Counts: 900 calls of `gp_krig_draw()` over three runs of 300 draws,
one factor each, 0 on the control; rank 7 to 15, median 10, over the
300 draws. Per draw the factor allocates `n x r` and vectors of `n`;
the `n x n` kernel and its `n^3` factorization are gone. Test: "the
kriging draw's factor is the conditional covariance"
(`tests/testthat/test-gp-multidim.R`, anisotropic and isotropic 2-D,
past, inside, a repeated and a near-duplicate position): `L L' + white`
equals `gp_krig_cov()` to `4 * gp_krig_tol * sd^2`, and the repeated
position draws one value. It errors on rellib-r5, which has no
`gp_krig_factor()`.

### m5. Memory of `extra_cov`

`extra_cov` is built from blocks (`extra_cov_assemble()`): nothing is
allocated when no source contributes (an empty `sparseMatrix`), a
sparse matrix when the blocks cover at most a quarter of it, dense
otherwise. `gp_krig_cov()` fills its result in column blocks of about
4 MB, upper triangle mirrored, so the result is the only `n x n` it
holds (it was seven: `outer()` alone makes three).
`frm_extra_cov_deriv()` returns the same kinds. Peaks on a 2000-row
grid, gc's max Vcells since a reset less what was in use then; gc
samples at collections, so a peak carries uncollected garbage and
moves between runs (`dev/gpby-p1-mem.R`, `dev/gpby-p1-mem2.R`,
`dev/gpby-p1-mem3.R`, one `n x n` is 30.5 MB):

| call | base | lane, review | lane, final |
|---|---|---|---|
| `frm_lp_basis(gp, extra_cov = TRUE)` | (no argument) | 386 MB peak | 73 to 106 MB own, result 31 MB |
| `frm_lp_basis(s(x), extra_cov = TRUE)` | (no argument) | 2 dense `n x n` | 3.4 MB own, empty sparse |
| `frm_curve(gp, nsim = 1000)` | 192 to 199 MB | 398 MB peak | 223 to 255 MB own |
| `frm_curve(s(x), nsim = 1000)` | 191 MB | 332 MB peak | 199 MB own |
| `frm_curve_deriv(s(x), nsim = 1000)` | 467 MB | | 467 MB own |

The gp curve's remaining excess over base is the kriging covariance
the band now carries and its sum with `A V A'`.

### m6. Band numbers on the final build

`dev/gpby-crit.R lane` (`dev/gpby-p1-crit-lane.txt`):

    GRID [1.0, 5.0] r in [0.00026, 0.00045] | frm_curve crit 2.66513 | A V A' self-standardized 2.66127 | ratio 1.00145
    GRID [5.5, 6.5] r in [0.00011, 0.00207] | frm_curve crit 2.28036 | A V A' self-standardized 2.27772 | ratio 1.00116
    GRID [7.0, 12.0] r in [0.01255, 0.69521] | frm_curve crit 2.43076 | A V A' self-standardized 2.33046 | ratio 1.04304

The first two grids move between builds (recorded 2.6654, 2.2812; the
review's build 2.6698, 2.2880) with the `A V A'` column beside them,
which no extra-covariance change reaches (2.66473 to 2.66127): on a
near-singular covariance the eigen basis of `crit_of()` rotates under
rounding-level changes of the design (the kriging weights now come
through `chol(K)`), and with it which 20000 draws are taken. The grid
past the data, where the band's change matters, is 2.43076 against the
review's 2.4308.

`dev/gpby-cov-band.R lane` over seeds 1 to 400, then
`dev/gpby-cov-summarise.R` (`dev/gpby-p1-cov-summary.txt`):

    == arm lane (lib: C:/Users/adf44/source/r/wt-gpby-lib/frmtmb C:/Users/adf44/source/r/wt-gpby-lib/frmtmb.spline )
      seeds 1 to 400, 400 distinct, 400 fitted, 0 failed
      simultaneous whole-curve coverage 324 / 400 = 0.8100, binomial mcse 0.0196, Wilson 95% (0.7687, 0.8454)
      pointwise per-point coverage 0.8866
      critical value: median 2.8038, range (2.5797, 3.0783)
    == paired on 400 seeds: crit lane / base median 1.0383, range (0.9962, 1.0994); covered by lane only 13, by base only 0; pointwise identical on 400

frmtmb.spline's NEWS carries 2.431.

### m7. BREAKING marks in NEWS

Core `NEWS.md`, Breaking changes: `get_prior()` lists no class `"sd"`
rows for a `gp()` term; `summary()$gp` lengthscale rows take brms's
names (`lscale(gpxzx)` for `lscale(gpxz[1])`; the `sdgp(gpx)` row was
already brms's, `dev/gpby-p1-sumgp-{lane,base}.txt`). New features:
`frm_extra_cov_deriv()`, and `extra_cov`'s sparse, empty and dense
forms (the `extra_white` element is gone). frmtmb.spline `NEWS.md`,
Breaking changes: `"Sigma"` is the whole covariance, `A V A'` plus
`extra_cov`; `frm_curve_deriv()` leaves out the nugget and takes the
extra part's derivative from `frm_extra_cov_deriv()`. frmtmb.sample
`NEWS.md`: the draw's nugget and its cost.

### m8. The class-`"sd"` refusal's separator

On `y ~ gp(x)` with `set_prior(class = "sd")` the message now ends
"This model has no random-effect standard deviations. A gp() term's
standard deviation is class "sdgp", as in brms", and with a group term
beside it "... (no resp, dpar or nlpar). A gp() term's standard
deviation is ..." (`dev/gpby-p1-m8.txt`).

### m9. Absolute correlation bounds

`test-gp-by-draws.R`: a variance ratio within four Monte Carlo
standard errors, `4 * sqrt(2 / (n - 1))`; a correlation of two
different fields within `4 / sqrt(n)`; one field's correlation within
`4 * (1 - rho^2) / sqrt(n)` of `rho` read off `extra_cov`.
`test-gp-by.R` (the scaled range of a sub-GP's inputs): within eight
ulps of the 1 it is scaled to.

### m10. The `test-fd-chain.R` reference

The reference is now the per-coefficient route at step 1e-7
(`fdc_plain()`), truncation about 1e-9 and rounding about 2e-9, with
the bound `fdc_tol <- 1e-7` derived from those in the comment. The
review's fixture, `gp(x) + (1 | g)` at a seen level, is in the second
test: the coefficient part of the chain route's `Est.Error`,
`sqrt(Est.Error^2 - (dP/deta)^2 extra_var)`, against `fdc_plain()`
within 1e-6. pass=11 fail=0.

### m11. Check C on `gp()` in sigma and in a multivariate model

`tests/testthat/helper-brms.R`: brms's dpar-suffixed GP data names
(`Kgp_sigma_1`, `Dgp_`, `Igp_`, `dmax_`), a block's prefix matched by
its dpar or by its response (`brms_prefix_is()`), and response-
suffixed latent names (`zgp_y_1_1`) in `brms_inner_pat`. New gated
test in `test-gp-by.R`: `ysig ~ x, sigma ~ gp(x, by = f, k = 6)`,
`mvbf(y ~ gp(x, by = f, k = 8), y2 ~ gp(x, k = 6))` and both with
`by = f`, `set_rescor(FALSE)`, each `brms_lp_check(joint = TRUE)` on
the review's data. Gated `test-gp-by.R`: pass=82 fail=0.

### Found, not mine: `disc ~ 0 + gp(x)`

Not local to gp: `disc ~ 0 + (1 | g)` stops the same way, while
`disc ~ 0 + s(x)` (it keeps a fixed column), `sigma ~ 0 + gp(x)` and
`sigma ~ 0 + (1 | g)` work. `est$betad` is NULL when `disc` has no
fixed column, and `lp_eta_design()` multiplies by it
(`dev/gpby-p1-disc.R`, `dev/gpby-p1-disc2.R`, `.txt`). Filed in
`dev/test-backlog.md`, Open - medium.

### For consolidation

The fixes lane's `s()` basis (`diagonal.penalty = TRUE`) moves the
truncation error of references on smooth fixtures. On the merged tree,
rerun: `test-fd-chain.R` (the `fdc_tol` bound on the `s()` cells),
`test-predict-re-uncertainty.R` (the three `sqrt(eps)` bounds),
`test-emmeans.R`, `test-brms-suite-emmeans.R`, `test-prior-compat.R`,
`test-brms-priors.R`, `test-get-prior-route.R` and frmtmb.eam's
`test-family.R`. The new tests of this round use `gp()` and group
terms only and do not reach `s()`.

### Tests and checks on the final build

Affected files one per process, gated (`dev/gpby-p1-tests.sh`,
`dev/gpby-p1-tests-b/`): test-brms-likelihood 543, test-emmeans 56,
test-fd-chain 11, test-gp-by 82, test-gp-multidim 52,
test-predict-re-uncertainty 78, test-prior-compat 198, frmtmb.sample
test-gp-by-draws 20, frmtmb.spline test-curve 68, test-deriv 48,
test-difference 61, test-new-levels 57; all fail=0 error=0 warn=0.

Plain tier of every file, `dev/gpby-suite.sh plain lane`
(`dev/gpby-suite-plain-lane/SUMMARY.txt`):

    tier plain arm lane, 341 files
    frmtmb files=204 pass=14173 fail=0 error=0 skip=171 warn=0
    frmtmb.coupling files=11 pass=542 fail=0 error=0 skip=5 warn=0
    frmtmb.eam files=29 pass=1743 fail=0 error=0 skip=3 warn=0
    frmtmb.latent files=10 pass=359 fail=0 error=0 skip=2 warn=0
    frmtmb.learn files=15 pass=429 fail=0 error=0 skip=13 warn=0
    frmtmb.ode files=11 pass=547 fail=0 error=0 skip=1 warn=0
    frmtmb.sample files=46 pass=2319 fail=0 error=0 skip=4 warn=0
    frmtmb.spline files=15 pass=566 fail=0 error=0 skip=1 warn=0

Gated tier of frmtmb.sample, frmtmb.spline and core's
`test-brms-suite-*.R` (`dev/gpby-p1-gated.sh`,
`dev/gpby-p1-gated/SUMMARY.txt`):

      BAD RESULT test-brms-suite-emmeans.R pass=10 fail=1 error=0 skip=0 warn=0
      BAD RESULT test-brms-suite-methods.R pass=161 fail=1 error=0 skip=0 warn=0
    frmtmb files=10 pass=431 fail=2 error=0 skip=0 warn=0
    frmtmb.sample files=46 pass=2535 fail=0 error=0 skip=1 warn=0
    frmtmb.spline files=15 pass=566 fail=0 error=0 skip=1 warn=0

The two failures are the two stale "now HOLDS" verdicts of section 5,
unchanged.

`R CMD check --as-cran` on the final tree, `dev/gpby-check.sh`, logs
`dev/gpby-p1-check-{core,frmtmb.spline,frmtmb.sample}.log`: core
Status: 1 NOTE, frmtmb.spline Status: 1 NOTE, both "Skipping checking
math rendering: package 'V8' unavailable"; frmtmb.sample Status: OK.

## 1. What brms 2.23.0 does, read off its own code and output

`dev/gpby-brms-explore.R` (log `dev/gpby-brms-explore.txt`) runs
`stancode()`, `standata()` and `default_prior()` on nine designs, and
`dev/gpby-brms-newlevel.R` (log `dev/gpby-brms-newlevel.txt`) runs
`posterior_epred()` on brms's stored `brmsfit_example6`.

* `gp(..., by = NA, k = NA, cov = "exp_quad", iso = TRUE, gr = TRUE,
  cmc = TRUE, scale = TRUE, c = 5/4)`. **`iso = TRUE` is the default**:
  `gp(x, z)` declares `vector[1] lscale_1`. frmtmb's default was
  `iso = FALSE` with a comment calling that brms's; it was wrong.
* A factor `by` is `Kgp` sub-GPs, one per column of
  `model.matrix(~ 0 + byval)` (`cmc = TRUE`) or `~ 1 + byval`
  (`cmc = FALSE`). Sub-GP `j` sees only rows `Igp_j = which(Cgp != 0)`,
  is multiplied by `Cgp`, and is scaled (`dmax`), centered and bounded
  (`L`) over those rows alone; `gr = TRUE` groups them to distinct
  positions (`Jgp`). Each has its own `sdgp_1[j]` and `lscale_1[j]`.
  A numeric `by` is one GP multiplied by `Cgp` on every row.
* Names: `sdgp_<sfx1>`, `lscale_<sfx2>`, `zgp_<sfx1>[i]`, with
  `sfx1 = "gp" + rename(vars) + rename(by) + level` (`gpxfa`,
  `gpxfIntercept` under `cmc = FALSE`) and, for a non-isotropic term,
  `sfx2 = outer(sfx1, covars, paste0)`, levels varying fastest
  (`lscale_gpxzfax`, `lscale_gpxzfbx`, ...). `summary()$gp` rows are
  the same names with `sdgp(...)` and `lscale(...)`.
* Defaults: `sdgp` takes `def_scale_prior`, `student_t(3, 0, 2.5)` on
  this data; `lscale` takes `inv_gamma` per sub-GP from
  `def_lscale_prior()`, tuned to that sub-GP's scaled distances.
* **A by-level the fit never saw is refused**, by `validate_newdata()`:
  "New factor levels are not allowed. Levels allowed: '0', '1'. Levels
  found: '0', '1', '2'", with or without `allow_new_levels = TRUE`. A
  missing `by` column is refused too.
* **Prediction at a new position draws the field**: on fit6, the same
  five draws of `posterior_epred()` at three new `Age` values under two
  seeds differ by up to 1.567822 (`dev/gpby-brms-newlevel.txt`), so
  brms draws from the conditional (`brms:::.predictor_gp_new()`, an
  `rmulti_normal()` of `K** - v'v + 1e-8`), not its mean.

## 2. `gp(x, by = f)` in core

**What changed.** `parse_gp_call()` (`R/parse.R`) takes brms's argument
list, exact names only (`aa$c` would partial-match `cmc`). The frame
(`R/frame.R`, the gp section of the component builder) builds one gp or
hsgp block per sub-GP from `gp_sub_terms()` (`R/gp-by.R`): its rows,
its multiplier, its own `dmax` (1 under `scale = FALSE`), center and
boundary over its own rows (distinct ones under `gr = TRUE`), brms's
`sfx1`/`sfx2` and the written term it belongs to (`gp_brms`). The
objective, the Laplace step and the sampler need no change: the blocks
are independent. Prediction (`pred_design()`, `R/predict.R`) reads the
multiplier off newdata with `gp_by_mult()`, which refuses an unseen
level with brms's message. `gp_brms_terms()` and `gp_brms_values()`
give brms's names and values to `variables()`, `hypothesis()` (through
`hyp_env_vals()`) and `summary()$gp` (`summary_gp_frame()`).

**Likelihood, against brms's Stan program.** Check C of
`dev/brms-likelihood-tests.md`: brms's `log_prob` at frmtmb's estimates
(flat priors) against frmtmb's joint density, and brms's gradient in its
latent `zgp`. `dev/gpby-lpcheck.R`, log `dev/gpby-lpcheck.txt`:

    RES y ~ gp(x, by = f, k = 8) | const -2.842e-14 | max_grad 2.776e-15 | ours -76.853269
    RES y ~ gp(x, by = w, k = 8) | const 0.000e+00 | max_grad 9.992e-16 | ours -83.245140
    RES y ~ gp(x, by = f, k = 8, cmc = FALSE) | const -2.842e-14 | max_grad 9.215e-15 | ours -76.998211
    RES y ~ gp(x, by = f, k = 8, gr = FALSE) | const -2.842e-14 | max_grad 1.887e-15 | ours -76.863579
    RES y ~ gp(x, z, by = f, k = 5, iso = FALSE) | const -2.842e-14 | max_grad 2.554e-14 | ours -105.394090
    RES y ~ gp(x, z, by = f, k = 5) | const 0.000e+00 | max_grad 1.416e-14 | ours -105.809134

The same check is a test (`test-gp-by.R`, gated) for the first two. The
translator in `tests/testthat/helper-brms.R` now maps `Kgp > 1`
(`brms_gp_term()`, per-sub-GP `dmax_<i>_<j>`, `zgp_<i>_<j>`), and
declares one lengthscale for an isotropic multi-dimensional term, which
it did not (the old rule repeated it per dimension and only ever ran in
one dimension).

The EXACT form cannot pass check C, for the nugget divergence already
recorded (frmtmb `1e-6` on the correlation, brms `1e-12` absolute). Its
validation is structural against brms's `standata()` (each sub-GP's
`Igp`, `Xgp` scaled to unit range, and for the Hilbert-space form the
basis rows and `slambda`, all equal to 1e-12) and against the closed
form: the gaussian marginal likelihood of independent sub-GPs with
frmtmb's nugget equals `-logLik()` to 1e-8 relative, and an `optim()`
of it from a perturbed start finds no better value (factor `by`,
numeric `by`, `cmc = FALSE`; `test-gp-by.R`).

**Names.** `variables(y ~ gp(x, by = f))` is `b_Intercept, sdgp_gpxfa,
sdgp_gpxfb, sdgp_gpxfc, lscale_gpxfa, lscale_gpxfb, lscale_gpxfc,
sigma`, and `summary()$gp` rows are `sdgp(gpxfa)` ... `lscale(gpxfc)`;
the lengthscale divided by each level's own largest distance, to 1e-12.
A plain `gp(x)` fit now lists `sdgp_gpx` and `lscale_gpx`; on 0.67.0
`variables()` listed neither.

**Priors.** Classes `"sdgp"` and `"lscale"` are native (`R/priors.R`,
`resolve_gp_class()` in `resolve_priorlist()`), with `coef` brms's sub-GP
name, on brms's natural scales with the log-Jacobian. The lengthscale
density sits at `theta - log(dmax)` through a constant entry `shift`
(`entry_offset()`, and `prior_draw_to_internal()` for `frm_simulate()`).
The MAP objective of a fit with a `normal()` on one sub-GP's lscale and
an `exponential()` on another's sdgp is the likelihood less the
hand-computed log densities to 1e-8 relative (`test-gp-by.R`). Class
`"sd"` no longer reaches a gp block, as in brms. `get_prior()` lists
brms's rows. frmtmb.sample's defaults are brms's: `dev/gpby-priors.R`
(log `dev/gpby-priors.txt`) compares every `sdgp`/`lscale` row of brms's
`default_prior()` with `get_prior(route = "sample")` on nine designs:

    ROWS 63 compared, 0 differ

(brms writes an empty string on a per-coefficient row with no prior and
frmtmb `(flat)`; the comparison reads both as flat.) `def_lscale_prior()`
calls `nleqslv`, which frmtmb.sample does not import;
`brms_lscale_prior()` solves the same two equations by Newton's method
with step halving, and every inv_gamma above agrees to brms's six
printed decimals.

**`conditional_effects()`.** A `y ~ gp(x)` fit had "No plottable
predictors found for dpar 'mu'" on 0.67.0 (`dev/gpby-ce.R`). `gp()`
covariates and the `by` variable are now effects and a multi-variable
term adds its pairs, as brms's `get_all_effects_type(x, "gp")` does:
`x`, `f`, `x:f` for `gp(x, by = f)`; the `x:f` display is each level's
own curve (`test-gp-by.R`).

**Sampling.** `frm_sample()` on `gp(x, by = f)` runs, and its draws carry
`sdgp_gpxfa`, `lscale_gpxfa`, ... on brms's scales (`gp_brms_natural()`,
exported for frmtmb.sample beside `brms_par_labels()`), equal to
`exp(theta)` and `exp(theta) / dmax` to 1e-12 (`test-gp-by-draws.R`).
`zgp` is not mapped: frmtmb samples the field itself, `b[i]`, which is
not brms's standardized `zgp`.

## 3. The joint covariance of the predicted field (item 2)

**The seam.** `frm_lp_basis(extra_cov = TRUE)` returns `extra_cov`, the
`n x n` covariance over the rows of everything in `extra_var`, and
`extra_white`, the nugget's share. `pred_design()` keeps each exact
sub-GP's kriging pieces (`krig`: the unseen rows, their positions, the
weights, the cross kernel, in correlation units); `gp_krig_cov()` forms
`sd^2 (k(X, X) - Xw Ks')` times the multipliers; `lp_extra_cov()` adds
every new level's `M S M'` over the rows its `extra_var_blocks()` key
reaches (the `new_key` grouping), so one unseen level is one draw and two
are two. Its diagonal is set to `lp_extra_var_vec()`, so a band and
`se.fit` read the same number to the bit. Two rows at one position get
the same residual exactly (`S[same] <- s`), which is what makes a
difference at one position cancel to the bit rather than to a rounding.
Rows at an observed position now carry the indicator exactly instead of
kriging weights that equal it to roundoff.

**(a) The simultaneous band.** `frmtmb.spline` draws from
`A V A' + extra_cov` and standardizes by its diagonal.
`dev/gpby-crit.R` rebuilds the R9 table of
`dev/reviews/2026-09-08-diffcurve.md` (whose data generator was not
kept), `y ~ fac + gp(x)`, 60 points on [0, 6], noise 0.2,
`nsim = 20000`, seed 1. Logs `dev/gpby-crit-base.txt` and
`dev/gpby-crit-lane.txt`:

    base GRID [1.0, 5.0]   r in [0.00026, 0.00045] | frm_curve crit 2.66436 | A V A' self-standardized 2.66473 | ratio 0.99986
    base GRID [5.5, 6.5]   r in [0.00011, 0.00207] | frm_curve crit 2.28636 | A V A' self-standardized 2.28660 | ratio 0.99990
    base GRID [7.0, 12.0]  r in [0.01255, 0.69521] | frm_curve crit 2.02211 | A V A' self-standardized 2.33046 | ratio 0.86769
    lane GRID [1.0, 5.0]   r in [0.00026, 0.00045] | frm_curve crit 2.66538 | A V A' self-standardized 2.66473 | ratio 1.00024
    lane GRID [5.5, 6.5]   r in [0.00011, 0.00207] | frm_curve crit 2.28116 | A V A' self-standardized 2.28660 | ratio 0.99762
    lane GRID [7.0, 12.0]  r in [0.01255, 0.69521] | frm_curve crit 2.42491 | A V A' self-standardized 2.33046 | ratio 1.04053

Extrapolating, the shipped critical value was 0.868 of the reviewer's
lower bound (it was 0.851 on the review's data, 1.17485 inverted), and
the new one is 1.0405 of it, above it as the review predicted, since a
residual that decorrelates faster than the mean pushes the maximum up.

**Coverage.** `dev/gpby-cov-band.R`: per seed, the true field is one
draw of the model's own GP (sd 1, lengthscale 1) at 60 observed
positions on [0, 6] and 25 grid points on [4, 9] jointly;
`y = 0.5 + f + N(0, 0.2^2)`; the band covers when the truth is inside
at all 25 points. 400 seeds per arm, 4 processes of 100, `nsim = 4000`
with the seed as the simulation seed. `dev/gpby-cov-summarise.R`, log
`dev/gpby-cov-summary.txt`, pasted verbatim:

    == arm base (lib: C:/Users/adf44/source/r/rellib-r5/frmtmb C:/Users/adf44/source/r/rellib-r5/frmtmb.spline )
      seeds 1 to 400, 400 distinct, 400 fitted, 0 failed
      simultaneous whole-curve coverage 311 / 400 = 0.7775, binomial mcse 0.0208, Wilson 95% (0.7342, 0.8155)
      pointwise per-point coverage 0.8866
      critical value: median 2.6976, range (2.5598, 2.8123)
    == arm lane (lib: C:/Users/adf44/source/r/wt-gpby-lib/frmtmb C:/Users/adf44/source/r/wt-gpby-lib/frmtmb.spline )
      seeds 1 to 400, 400 distinct, 400 fitted, 0 failed
      simultaneous whole-curve coverage 324 / 400 = 0.8100, binomial mcse 0.0196, Wilson 95% (0.7687, 0.8454)
      pointwise per-point coverage 0.8866
      critical value: median 2.7996, range (2.5787, 3.0753)
    == paired on 400 seeds: crit lane / base median 1.0379, range (0.9958, 1.0983); covered by lane only 13, by base only 0; pointwise identical on 400

Paired, the change only adds coverage: 13 seeds the new band covers and
the old did not, none the other way, and the pointwise band is
identical on all 400 (it already carried the variance). The rest of the
shortfall from 0.95 is not the band's: the pointwise band covers 0.887
per point, which is plug-in maximum likelihood on two GP
hyperparameters from 60 points, and a band built on those estimates
inherits it. This design was chosen to be the band's own failure mode
(half the grid past the data); I did not run a design that separates
the band from the plug-in error, such as known hyperparameters, and
the coverage numbers should be read as before against after, not as
calibration. The lane arm was run twice; the first run overlapped a
core reinstall into the lane library, so it was discarded and rerun
on a stable library, with identical counts.

**(b) A difference at different exact positions** is computed from one
call on the stacked grid: `E11 + E22 - E12 - E21`. `sp_same_latent()`
is deleted. Against the closed form built from the fitted kernel by
hand, `test-difference.R` asserts agreement to 1e-6 at shifts 0.05 and
1, and on the mirrored grid (kriging variances equal to 1.11e-16,
different draws) where the residual's share is real. A shift of 1e-10
differs from the one-position difference by exactly the two nuggets,
`2e-6 sd^2` in variance (to 1e-3 relative): the nugget is white noise
and does not cancel between two positions. The by-factor smooth beside
a `gp()`, which the old predicate refused, differences to the
coefficient part alone to 1e-10. The same seam removes the refusal of
every route under `allow_new_levels = TRUE` (`sp_new_level_stop()` is
deleted): two different unseen levels add their variances to 1e-10, one
unseen level on both sides cancels, a random intercept leaves the
derivative's standard error unchanged (1e-6) and a random slope adds its
variance to it (1e-4 against the slope variance read off `extra_var` at
three values of `t`).

The derivative and a peak's standard error use `extra_cov - extra_white`:
the nugget is white noise and a finite difference of it divides
`2e-6 sd^2` by `4 e^2`, which with `e = 1e-6` of the range is enormous.

**emmeans.** `emm_basis_grid()` adds each part's `extra_cov` to `V`, and
`emm_check_part()` no longer refuses an exact `gp()` at an unseen
position. A mean at `mean(x)`, which no observation holds, has
`frm_linpred()`'s standard error to 1e-8, and the contrast between two
levels at that position has `sqrt(vcov()["fB", "fB"])` to 1e-8
(`test-emmeans.R`; refused on 0.67.0).

**(c) The finite-difference `Est.Error`.** `fit_fd_se()` takes
`b_chain`: for a smooth, `gp()` or `hsgp()` block the derivative of the
output in the block's coefficients is `(d out / d eta) Z`, and
`d out / d eta` on every row is ONE pair of evaluations with that
predictor's `eta` moved by the step (`fit$eta_shift`, read in
`lp_eta_design()`). The same derivative carries `extra`, the kriging and
new-level variance, which this route left out. The batched, chain and
per-coefficient columns now mix in one quadratic form (before, one
unbatched coefficient sent every coefficient to the per-coefficient
route). `dev/gpby-fdcost.R` counts every call of the fitted value,
estimate included, on the cells of `dev/resmooth-batchcost.R` with
their seeds (logs `dev/gpby-fdcost-base.txt`, `dev/gpby-fdcost-lane.txt`):

    base A gp(x) n_b=160, newdata 3 NEW positions, NA evals  330 | 1.8700 s | Est.Error[1, 2, 1] 0.08182400254
    base B gp(x) n_b=160, in sample, NA               evals   12 | 0.0256 s | Est.Error[1, 2, 1] 0.09957996612
    base C gp(x) n_b=160, newdata 3 observed, NA      evals   12 | 0.0391 s | Est.Error[1, 2, 1] 0.09957996612
    base D s(x,k=8)+(1|g) 40 lv, in sample, NULL      evals   26 | 0.1525 s | Est.Error[1, 2, 1] 0.09224894051
    base E s(x,k=8)+(1|g) 40 lv, in sample, NA        evals   24 | 0.1587 s | Est.Error[1, 2, 1] 0.03776644314
    base F hsgp(x) k=12, in sample, NA                evals   34 | 0.0484 s | Est.Error[1, 2, 1] 0.07602082
    base CONTROL x+(1|g) 40 lv, in sample, NULL       evals   12 | 0.0525 s | Est.Error[1, 2, 1] 0.05414304746
    lane A gp(x) n_b=160, newdata 3 NEW positions, NA evals   12 | 0.1156 s | Est.Error[1, 2, 1] 0.08182412751
    lane B gp(x) n_b=160, in sample, NA               evals   12 | 0.0222 s | Est.Error[1, 2, 1] 0.09957996612
    lane C gp(x) n_b=160, newdata 3 observed, NA      evals   12 | 0.0203 s | Est.Error[1, 2, 1] 0.09957996612
    lane D s(x,k=8)+(1|g) 40 lv, in sample, NULL      evals   16 | 0.0750 s | Est.Error[1, 2, 1] 0.09224894051
    lane E s(x,k=8)+(1|g) 40 lv, in sample, NA        evals   14 | 0.0563 s | Est.Error[1, 2, 1] 0.03776644314
    lane F hsgp(x) k=12, in sample, NA                evals   12 | 0.0184 s | Est.Error[1, 2, 1] 0.07602082
    lane CONTROL x+(1|g) 40 lv, in sample, NULL       evals   12 | 0.0459 s | Est.Error[1, 2, 1] 0.05414304746

330 is resmooth's 329 plus the estimate's own call. Cells B to F and the
control print the same `Est.Error` to 10 digits on both arms, so the
chain rule is the per-coefficient derivative there (the test asserts
1e-6 relative). Cell A moves in the seventh digit, which is the kriging
variance at three positions inside the data; past the data it is not
small: on `test-fd-chain.R`'s fixture at `x = 2.6` and `3.2`, 0.67.0's
`Est.Error` is 8 percent short of `sqrt(plain^2 + (dp/deta)^2 extra_var)`
with `dp/deta` written out for the cumulative logit, and the lane's
agrees with it to 1e-5. The wall clock in cell A, 1.87 s to 0.116 s, is
still 2.3x 0.64.0's 0.0497 s (the count is 12 against 0.64.0's 11 for
the same reason as the 330: one call more than that record counted).
The analytic route `ord_prob_se()`, which `conditional_effects()` takes
for an ordinal family, had the same missing term and now adds it.

**Against brms's posterior predictive.** `dev/gpby-brms-epred2.R` fits
`y ~ gp(x, by = f)` in brms (2 chains of 1500, 750 warmup, seed 2026,
the test data) and maps every one of its 1500 draws into frmtmb's
parameters (`theta = log sdgp, log(lscale * dmax)`, `b = L zgp` with
brms's own kernel). frmtmb.sample's per-draw predictor now draws the
field at an unseen position from its conditional law
(`fit$krige_draw`, set by `draws_fit_at()`). Log
`dev/gpby-brms-epred2.txt`, the spread across the same 1500 draws:

    row 1 x 2.55 f a | 0.39523 | 0.39526 | 1.0001 | 0.39523 | 1.0000
    row 2 x 6.50 f a | 0.57846 | 0.57648 | 0.9966 | 0.56727 | 0.9806
    row 3 x 7.50 f a | 0.93147 | 0.90733 | 0.9741 | 0.80167 | 0.8607
    row 4 x 2.55 f c | 0.38052 | 0.38069 | 1.0004 | 0.38066 | 1.0004
    row 5 x 6.50 f c | 0.48870 | 0.48671 | 0.9959 | 0.46664 | 0.9549
    row 6 x 7.50 f c | 0.63475 | 0.63774 | 1.0047 | 0.57599 | 0.9074

(columns: brms `posterior_epred()` less its intercept | frmtmb with the
draw | ratio | frmtmb with the conditional mean only, the 0.67.0 route |
ratio). With the draw the spread is 0.974 to 1.005 of brms's; without it
0.861 to 1.000. Draw by draw the two conditional LAWS are not the same,
and the cause is the nugget divergence again: in units of brms's
posterior sd at the row, the conditional sd at `x = 7.5`, `f = a` is a
median 0.107 in brms and 0.363 in frmtmb, and the conditional means
differ by a median 0.157. A `1e-12` nugget lets a long-lengthscale
field be extrapolated from components that `1e-6` treats as noise. The
two move spread between the mean and the draw and agree on the total.
The test pins frmtmb's own law instead: 1500 repetitions of one draw
have the conditional variance `frm_lp_basis()` reports to 0.15 (four
relative standard errors of a variance from 1500 draws) and the
correlation of two positions of one field to 0.1
(`test-gp-by-draws.R`). On 0.67.0 the variance across repetitions was
exactly zero.

A side fix found here: kriging solved `K` with the sd in it, so a
sampled draw whose sd underflows (seen on a fit-route chain that stuck;
Punch round 1, m1) made the
solve exactly singular. The weights are computed in correlation units
now.

## 4. Tests seen to fail on rellib-r5

Logs `dev/gpby-base-tests/*.txt` (lane logs `dev/gpby-core-tests/`,
`dev/gpby-spline-tests/`). Behavioral failures on base, by test:

* `test-gp-by.R` "variables() lists a gp() term's sdgp_ and lscale_":
  `variables(f0)` was `b_Intercept, sigma`. "conditional_effects()
  draws a gp() term's covariates": "No plottable predictors found". The
  other blocks of the file fail on base at `gp(x, by = f)` itself, the
  weak form.
* `test-gp-multidim.R`: the default `gp(x1, x2)` fit had logLik
  different from the `iso = TRUE` fit and three `confint_varcorr()`
  terms instead of two.
* `test-fd-chain.R`: 250 and 34 evaluations against 12 and 12; and the
  kriging term, `Est.Error` 0.07994 relative short of the closed form,
  not above the plain route.
* `test-emmeans.R`: "emmeans cannot use this reference grid: an exact
  gp() term is predicted at a position the fit did not see".
* frmtmb.spline `test-curve.R` (new block): `sqrt(diag(Sigma))` 0.448
  relative away from `.se`, and the critical value 2.02 below the bound
  2.28. `test-difference.R`: both new blocks refused by the
  `sp_same_latent()` message. `test-new-levels.R`: "a difference across
  two unseen levels" refused; the "every route" block errors on
  `extra_cov` (weak form; its behavioral failure is the record in
  `dev/reviews/2026-09-29-splinecurve.md`, Finding A).
* frmtmb.sample `test-gp-by-draws.R` "a plain gp() draws the field past
  its positions too": the variance across 1500 repetitions of one draw
  was 1.00 relative short of the conditional variance, that is zero.

## 5. Ledger rows

`dev/brmsport-verdicts.tsv` is not edited. The fixture fit6 of both
`helper-brms-suite.R` files is now brms's own formula,
`volume ~ Trt + gp(Age, by = Trt, gr = TRUE)`. The gated run (section
6) shows exactly two stale verdicts, both "now HOLDS":

* `brmsfit-methods:991` (`lscale_volume_gpAgeTrt0`,
  `lscale_volume_gpAgeTrt1` in `variables(fit6)`), recorded "cannot
  transfer ... needs gp(Age, by = Trt)".
* `emmeans:42` (`emmeans(fit6, "Age", by = "Trt", epred = TRUE)`),
  recorded "cannot transfer ... the kriging covariance between two grid
  points is not available".

No other row that reads fit6 changed its outcome.

## 6. Suites

`dev/gpby-suite.sh` runs every test file of core and the seven
extensions, one R process per file, through `dev/gpby-runtest.R` (which
sums `failed` and `error` and prints `BAD` lines). Sibling extensions
I did not change load from `rellib-r5` with this lane's core first,
which the `lib:` lines show. Summaries generated by the script, pasted:

Plain tier (`NOT_CRAN=true`), `dev/gpby-suite-plain-lane/SUMMARY.txt`:

    tier plain arm lane, 341 files
    frmtmb files=204 pass=14161 fail=3 error=0 skip=170 warn=0
      BAD RESULT test-predict-re-uncertainty.R pass=75 fail=3 error=0 skip=0 warn=0
    frmtmb.coupling files=11 pass=542 fail=0 error=0 skip=5 warn=0
    frmtmb.eam files=29 pass=1743 fail=0 error=0 skip=3 warn=0
    frmtmb.latent files=10 pass=359 fail=0 error=0 skip=2 warn=0
    frmtmb.learn files=15 pass=429 fail=0 error=0 skip=13 warn=0
    frmtmb.ode files=11 pass=547 fail=0 error=0 skip=1 warn=0
    frmtmb.sample files=46 pass=2319 fail=0 error=0 skip=4 warn=0
    frmtmb.spline files=15 pass=556 fail=0 error=0 skip=1 warn=0

The three failures were three assertions that two finite-difference
routes agree to 8 or 100 ulps. The chain rule is a different central
difference of the same derivative, and the measured gaps are 4.6e-11,
6.8e-11 and 2.7e-10 relative, the square of the 1e-5 step; the
defects those tests pin were 0.676, 0.52 and several-fold. The bound
is now `sqrt(.Machine$double.eps)` relative, with the reason in each
comment; rerun alone, `test-predict-re-uncertainty.R` pass=78 fail=0.

Gated tier (`FRMTMB_BRMS_FIT_TESTS=true` too),
`dev/gpby-suite-gated-lane/SUMMARY.txt`:

    tier gated arm lane, 341 files
    frmtmb files=204 pass=16451 fail=2 error=0 skip=14 warn=0
      BAD RESULT test-brms-suite-emmeans.R pass=10 fail=1 error=0 skip=0 warn=0
      BAD RESULT test-brms-suite-methods.R pass=161 fail=1 error=0 skip=0 warn=0
    frmtmb.coupling files=11 pass=542 fail=0 error=0 skip=5 warn=0
    frmtmb.eam files=29 pass=1743 fail=0 error=0 skip=3 warn=0
    frmtmb.latent files=10 pass=359 fail=0 error=0 skip=2 warn=0
    frmtmb.learn files=15 pass=500 fail=0 error=0 skip=2 warn=0
    frmtmb.ode files=11 pass=547 fail=0 error=0 skip=1 warn=0
    frmtmb.sample files=46 pass=2535 fail=0 error=0 skip=1 warn=0
    frmtmb.spline files=15 pass=556 fail=0 error=0 skip=1 warn=0

The two failures are the stale verdicts of section 5, which are
expected. The skips are `test-drmtmb-agreement.R` (13, its own gate,
`FRMTMB_DRMTMB_FIT_TESTS`), `test-fuzz.R` (1) and the extensions'
`test-scale.R` files (the scale tier's gate). Run with its gate,
`test-drmtmb-agreement.R` gives pass=131 fail=0 skip=0
(`dev/gpby-core-tests/test-drmtmb-agreement.txt`). The plain tier ran
before the three tolerances were changed and the gated tier after.

## 7. Not done, and found but not fixed

* **`cov = "matern32"`, `"matern52"`, `"exponential"`.** brms's other
  kernels are refused by name. The exact form needs the kernel in
  `gp_corr()` and `gp_cross_cov()`; the Hilbert-space form needs brms's
  spectral densities (`spd_gp_matern32()` and the rest, in
  `R/formula-gp.R`). Not started.
* **The exact `gp()` nugget** (`1e-6` on the correlation against brms's
  `1e-12`) is unchanged: it is the maintainer's decision recorded in
  `dev/brms-likelihood-tests.md`, and section 3 measures what it does to
  a prediction at a new position.
* **`zgp` is not a draw name.** frmtmb samples the field, `b[i]`, and
  brms its standardized `zgp`; mapping one to the other is the change of
  variables `L^-1 b`, which a draw column cannot be renamed into.
* **The draws' column order** puts each sub-GP's `sdgp` beside its
  `lscale` (template order); brms puts every `sdgp` first. Names match,
  order does not, because the sampler reads columns by position.
* **Core `predict()`** (the maximum-likelihood predictive draws) holds a
  `gp()` curve at its mode and draws neither its coefficients nor the
  kriging residual (`predict_b_drawer()` says so on purpose). brms's
  `predict()` draws both; core's interval past the positions is
  therefore narrower than brms's. Not changed here.
* **A fit-route chain on `gp(x)` can stick.**
  `frm_sample(fit, chains = 1, iter = 300, seed = 2)` on
  `y ~ gp(x)` (40 points, seed 9) returned every draw at `sdgp = 0`,
  `lscale = 3.1e56`, `sigma = 2.9e79`; seed 4 sampled. The reason
  given here first, that the fit route keeps flat priors, is wrong:
  Punch round 1, m1, measures the same priors on both routes and a
  step size of NaN on the fit route, and files it.
* A nonlinear body still refuses a contributing exact `gp()` at an
  unseen position (`lp_basis_nl()`); `extra_cov` is returned there as
  the diagonal the route allows, zero.

## 8. R CMD check

`dev/gpby-check.sh`, `R CMD check --as-cran` once per changed package,
built and checked in `dev/gpby-check/<pkg>/`, `R_LIBS` the lane library,
then `rellib-r5`, then the user library,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH:

* frmtmb: `Status: 1 NOTE` (the HTML manual's "Skipping checking math
  rendering: package 'V8' unavailable", the expected one).
* frmtmb.spline: `Status: 1 NOTE`, the same V8 note.
* frmtmb.sample: `Status: OK`.

## 9. Versions

Core needs a MINOR bump: new arguments to `gp()` and `frm_lp_basis()`,
new prior classes, a new sampling-API export, and two breaking defaults.
frmtmb.sample and frmtmb.spline each need a minor bump and a floor on
the development version of core (`gp_brms_natural()`, prior classes
`"sdgp"`/`"lscale"`, `krige_draw`; `frm_lp_basis(extra_cov = TRUE)`).
Their `DESCRIPTION` floors are not edited, per the brief.

## At the 0.68.0 release (2026-10-06)

The consolidation closed the review's remaining findings
(`dev/round-20261005.md`, "Small changes"):

- **r1, the nugget**: kept at `1e-6` by the user (above).
- **r2, the high-rank factor**: `gp_krig_factor()` now takes the dense
  route (the whole conditional covariance over the distinct positions
  and LAPACK's pivoted Cholesky, `gp_krig_factor_dense()`) once its
  pivoted loop passes a quarter of the positions. Measured on the
  release build, `dev/rel068-krig-rank.R` (`dev/rel068-log/krig-rank.txt`),
  interleaved arms, each the minimum of 5 rounds of blocks over 1.2 s,
  with a control that times the lane's function against a copy of
  itself:

  | length scale | rows | rank | lane | release | release / lane | control |
  |---|---|---|---|---|---|---|
  | 1.48 (fitted) | 1000 | 23 | 0.0054 s | 0.0059 s | 1.10 | 0.88 |
  | 0.05 | 100 | 100 | 0.0026 s | 0.0012 s | 0.45 | 0.98 |
  | 0.05 | 400 | 400 | 0.0709 s | 0.0195 s | 0.28 | 1.02 |
  | 0.05 | 1000 | 499 | 0.2925 s | 0.1612 s | 0.55 | 1.03 |

  The low-rank row runs the same code on both arms; its 1.10 is inside
  the control's spread at 5 ms (0.88). The law is unchanged: `L L' +
  diag(white)`, scaled, equals `gp_krig_cov()` to 8.8e-13 of `sd^2` on
  the dense route and 8.5e-13 on the pivoted one.
- **r3, stale text**: frmtmb.spline's `frm_curve` x `gp` row said 2.425;
  the release build gives 2.43076 (`dev/gpby-crit.R` with the release
  library, `dev/rel068-log/crit.txt`), so it says 2.431, as the NEWS
  does. The two flat-prior phrases above are corrected.
- **r4, the fixture chain**: the formula-route chain that
  `test-gp-by-draws.R` samples is poorly mixed (on the m1b construction,
  `iter = 600`, one chain: 215 transitions over the maximum tree depth
  and R-hat 1.56, `dev/gpby-p1-m1b-lane.txt`). Its assertions read the
  kriging law at one fixed draw, not the posterior, so they are not at
  risk; a later test that reads posterior summaries off that fixture
  must not trust it. Filed in `dev/test-backlog.md`.

