# Lane optima: fits that stop short of the maximum likelihood

Round of 2026-10-06, base 4f5ea39f (frmtmb 0.68.1), worktree
`frmtmb-wt-optima`, private library `wt-optima-lib`. "Base" is
`rellib-r6`. Every number below is generated from a log in
`dev/optima-log/` (gitignored) by the script named beside it;
`dev/optima-record.R` prints the main blocks again.

## Final check, c1: a weight at 1

A simplex at a vertex has a weight at 1, the complement of weights at
0, and the same degenerate delta method: on `ls ~ mo(income) * age`
(`optima_mo_data(1234)`) `moincome:age1[2]` is 1 - 5e-10 with a Wald
error of 1.79e-4 (`dev/optima-p3-vertex-probe.R`). `summary_mo_frame()`
now marks a weight within `mo_face_tol` of 0 or 1, so the face
attribute, the `NA` error and interval, the footnote ("A weight at 0
or 1 (...)") and `check_laplace()`'s naming cover every weight of a
simplex at a vertex. New tests: core `test-optima.R` "a simplex at a
vertex has no delta-method interval at 1", frmtmb.sample "check_laplace()
names a weight at 1 instead of judging it". On the punch-2 build:
test-optima.R pass 97 fail 3, test-mo-simplex-draws.R pass 24 fail 3
(`dev/release/optima-p3-before-files/`). On this build: pass 100 and
pass 27, fail 0, under both BLAS (`optima-p3`, `optima-p3-ob`).

## Punch round 2

### B2. `check_laplace()` on a `mo()` fit

The first round broke it: `draws_outer_cols()` listed D weights where
`fit$opt$par` has D - 1 coordinates, and it stopped with "length(ml)
== length(keep) is not TRUE" on every `mo()` fit (the review's
`dev/optima-rev2-laplace.R`, 4 of 4). Equal counts would not have been
enough: 106 of 200 simplexes of `ls ~ mo(income) * age` (seeds 1 to
100) have fitted coordinates on another sheet of the chart than the
one draws are carried back to.

- The simplex is compared on the weight scale: its rows are brms's
  `simo_` weights, `ml` and `wald_se` the weights and delta-method
  errors `summary()` reports (`summary_mo_frame()`, now exported to
  frmtmb.sample), `post_mean` and `post_sd` the posterior of the
  `simo_` draws. The coordinates leave both sides.
- What is added, stated in `?check_laplace`: the simplex is sampled in
  softmax coordinates under `dirichlet(1)`, the flat density on the
  simplex, under which the weights' posterior is the likelihood the
  fit maximized. Without a density there, the coordinates themselves
  are the prior. A MAP fit keeps its own prior: under `.diagnostic`
  the retape resolves the fit's own prior again, where the first round
  kept `fit$obj` and its improper chart.
- A weight at 0 has no Wald error (`NA`) and is named in one message
  instead of being flagged.

The review's `dev/optima-rev2-laplace.R` (it reads `wt-optima-lib`
first; stdout to `dev/optima-log/p2-rev2-laplace2.txt`): on ML and
MAP fit1 the simplex rows have |z_shift| at most 0.06 and `sd_ratio`
1.02 to 1.08; on the weak case the two interior weights are flagged
(Wald SE 0.69 against posterior sd 0.23, `sd_ratio` 0.33 to 0.34) and
`simo_mox1[3]` is named at the boundary. Its first run on this build
found a bug in that message: on a fit with no weight at a face it said
"simo_ sits at 0", because `paste0("simo_", character(0))` is
`"simo_"`. Fixed; the test asserts no boundary message on fit1.

Other functions that pair the fit's coordinates with draw columns: a
grep of `draws_outer_cols(`, `opt$par`, `last.par`, `cov.fixed` and
`sdr_of(` in frmtmb.sample finds `check_laplace()` alone. The other
`draws_outer_cols()` callers (`summary()`, `print()`, `mcmc_plot()`,
`pairs()`) read draws by name. Every reader that hands a draw to the
model goes through `draws_internal_matrix()`, the inverse the review
verified to 1e-14, under `laplace = TRUE` too. `log_posterior()` reads
the stanfit's own `lp__`. `frm_sample(start =)` is for the formula
route only. `bridge_sampler()` refuses all draws, as before.

The review's sampler script (`dev/optima-p2-rev2-sample.R`, its copy
with output to `dev/optima-log/p2-rev2-sample.txt`), against brms
2.23.0: names identical on both cases, 0 divergences on both. Weak
case: weight means 0.355/0.351/0.294 (brms 0.357/0.351/0.292), sds
0.236/0.239/0.220 (brms 0.236/0.233/0.221), |z_mean| at most 0.44.
Strong case: means 0.651/0.215/0.134 (brms 0.652/0.214/0.134),
|z_mean| at most 0.70, R-hat at most 1.0004.

### `as_tmbstan()`

The backlog entry is corrected: every `mo()` fit, ML or MAP, samples
the fit's own tape there, and a weakly identified simplex collapses to
the barycenter (sd 6e-18 to 1e-17, coordinates to 7.1e17, reproduced
in `p2-rev2-laplace2.txt`). By contract the route adds nothing. It now
says so in one message on a fit with `mo()` terms, and `?as_tmbstan`
has the sentence.

### n1, n2, n3

- n1. `summary()` prints `NA` for the error and interval of a weight
  at a face, and a line under the table says why (weak case:
  "A weight at 0 (mox1[3]) is on the boundary of the simplex, where
  the delta method has no standard error or interval to give;
  frmtmb.sample's frm_sample() gives a posterior one.").
  `summary_mo_frame()` carries the face rows in attribute `face`.
- n2. Each term's last weight goes right after its other weights, so
  with several `mo()` terms each term's weights are adjacent, in
  brms's order (`draws_to_natural()`, the block flag `adjacent`; the
  inverse finds columns by name).
- n3. Core DESCRIPTION is `0.68.1.9000`; frmtmb.sample imports
  `frmtmb (>= 0.68.1.9000)`. The punch-1 note is corrected below.

### Tests

frmtmb.sample `test-mo-simplex-draws.R` has four new tests
(`check_laplace()` on fit1 and on the weak case, column order, the
`as_tmbstan()` message). On the punch-1 build: pass 10, fail 5,
error 2 (`dev/release/optima-p2-before-files/`). On this build:
pass 24, fail 0 (`optima-p2c`). Core `test-optima.R` asserts the face
`NA`s and the footnote: pass 95 (`optima-p2b`).

Suites, one file per process, gates on, from
`dev/optima-suite-sum.R` (logs `dev/optima-log/suite-p2-ref.txt`,
`suite-p2-ob.txt`):

| tier | package | files | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|---|
| reference BLAS | frmtmb | 216 | 17356 | 1 | 0 | 0 | 0 |
| reference BLAS | frmtmb.sample | 49 | 2640 | 0 | 0 | 1 | 0 |
| OpenBLAS | frmtmb | 216 | 17355 | 2 | 0 | 0 | 1 |
| OpenBLAS | frmtmb.sample | 49 | 2640 | 0 | 0 | 1 | 1 |

Every non-pass is in punch round 1's runs too: the ledger row
brmsfit-methods:326 (both BLAS), and under OpenBLAS the Malingering_2
constant, the singular-convergence warning in `test-ordinal-mixture.R`
and an R-hat warning in `test-brms-shapes-draws.R`.

R CMD check --as-cran (`dev/optima-check/`): frmtmb 2 NOTEs (V8
unavailable, as before; "Version contains large components
(0.68.1.9000)", from n3). frmtmb.sample: OK.

### Files touched in punch round 2

`DESCRIPTION`, `NEWS.md`, `NAMESPACE`, `R/methods-fit.R`,
`R/sampling-api.R`, `man/frmtmb-sampling-api.Rd`,
`tests/testthat/test-optima.R`; `extensions/frmtmb.sample/`:
`DESCRIPTION`, `NEWS.md`, `R/sample.R`, `R/draws-brms.R`,
`man/check_laplace.Rd`, `man/as_tmbstan.Rd`,
`tests/testthat/test-mo-simplex-draws.R`; `dev/test-backlog.md`,
`dev/optima-findings.md`, `dev/optima-p2-rev2-sample.R` (new).

## Punch round 1 (review `dev/reviews/2026-10-07-optima.md`)

### B1. The sampler: softmax coordinates under `dirichlet(1)`

The review was right: frmtmb.sample put no density on `zeta`, so the
chart was the prior, and the change of chart changed every `mo()`
posterior; the earlier "The sampler" section below read as "unchanged"
and is corrected there. Fixed in frmtmb.sample, not documented:

- The maximum-likelihood chart stays for `frm()`. The sampler retapes
  the model with every simplex read in softmax coordinates
  (`mo_chart_frame()`, `mo_simplex(z, "softmax")`), one sheet with no
  fold, and adds brms's default `dirichlet(1)`: the uniform density on
  the simplex times the softmax's Jacobian `prod(w)`, so
  `-sum(log w) - lgamma(D)` (`mo_simplex_nlp()`). It is the flat
  density on the simplex, so `prior = "flat"` keeps it.
- At the boundary: starts move from the fit's chart to softmax
  coordinates (`mo_sample_pars()`, a weight at a face starts at 1e-8),
  the draws move back to the fit's chart (`mo_draws_to_chart()`), and
  `draws_to_natural()` reports them as brms's `simo_<label>1[k]`, D
  weights in place of the D - 1 `zeta<j>_k` columns, with the inverse
  for every reader that hands a draw to the model (`draws_mo_cols()`).
  This is the rename lane surface found (frmtmb.sample named the draws
  `zeta1_*` where brms says `simo_*`); it is done here, not in lane
  surface.
- `simo` prior rows stay refused: `set_prior()` has no Dirichlet
  density, so that is filed (dev/test-backlog.md, "Filed by lane
  optima after punch round 1").

Against brms 2.23.0, the review's two cases with its data code, 4
chains of 2000 (1000 warmup) on each side, seed 1
(`dev/optima-p1-sample.R lane` and `lane brms`, compared by
`dev/optima-p1-sample-cmp.R`, log `dev/optima-log/p1-sample-cmp.txt`):

```text
== weak: divergences frmtmb 0 brms 0
       weight mean_frmtmb mean_brms  z_mean sd_frmtmb sd_brms    z_sd  rhat_frmtmb rhat_brms
 simo_mox1[1]      0.3586    0.3533  0.9091    0.2382  0.2336  1.1506  1.002 1.0008
 simo_mox1[2]      0.3440    0.3512 -1.4260    0.2349  0.2355 -0.1421  1.002 1.0000
 simo_mox1[3]      0.2974    0.2955  0.3470    0.2251  0.2202  1.2189  1.000 0.9997
== strong: divergences frmtmb 0 brms 0
 simo_moincome1[1]      0.6519    0.6533 -1.6956   0.03485 0.03589 -1.3405  0.9999 1.001
 simo_moincome1[2]      0.2144    0.2136  0.6693   0.04178 0.04271 -0.8993  1.0007 1.001
 simo_moincome1[3]      0.1337    0.1331  0.4826   0.03885 0.04013 -1.2173  1.0033 1.002
largest |z| over every weight's mean and sd: 1.7
```

`z` is the difference over the joint Monte Carlo standard error. Every
weight's mean and sd agree with brms within 1.7 of it, with 0
divergences on both sides, and R-hat on `w` at most 1.0033 (review:
base every draw at a vertex with 44 divergences; the first lane build
every draw at the barycenter with `zeta` to 1.4e18; the chart under
`dirichlet(1)` 50 divergences and R-hat 8.47). The draws carry brms's
names: `b_Intercept bsp_moincome sigma simo_moincome1[1]
simo_moincome1[2] simo_moincome1[3] lp__`.

A model without `mo()` takes the same code path it took: three models
(`y ~ x + (1 | g)`, poisson, cumulative with `cs()`), default priors
and `prior = "flat"`, 2 chains of 600, each arm pinned to core 0
(`dev/optima-p1-nomo.R`, `-cmp.R`): draws identical on 6 of 6, max
|diff| 0.

### Minors

- m1. NEWS and the comments of `R/links.R` say "in value" now and give
  the fit movement the review measured: probit estimates by up to 9e-10
  relative (log-likelihood 3.4e-13), cloglog and softit by up to
  3.6e-12 (the adjoint's `f'(e) + 1 - 1`).
- m2. Text the change made false, rewritten: the SE warning's clause
  ("a standard deviation at 0 does this, and so do a mo() simplex split
  by a category no row is in and parameters that enter only through a
  combination"; `R/se-check.R`, `se_lost_clauses()`), the `check_se`
  documentation (`R/fit.R`, `man/frmtmb_control.Rd`),
  `vignettes/diagnostics.Rmd`, the comments of `R/se-check.R` (header,
  `se_flat_tol`'s, `jc_nonest()`'s), `test-brms-likelihood.R`,
  `test-brms-agreement.R` and `test-brms-names.R`. For consolidation:
  lane setier's `R/se-check.R` keeps the old clause (its line 1220);
  whichever lands second takes this wording.
- m3. The RTMB derivative curve is the review's (below, item 2), and
  `dev/upstream-bugs.md` RTMB-3 carries it.
- m4. `escape_stationary()` now settles the objective's state as
  `mo_search()` does (`obj_settle()`, new, used by both). Test
  "a skew normal escape reports the modes of its optimum"
  (`test-optima.R`, seed 20 of the review's
  `dev/optima-rev-escape.R`): seen to fail on rellib-r6 (modes 1.0531
  against 1.0764 at the optimum).
- m5. `frmtmb_control(mo_search = FALSE)` turns the search off, and
  NEWS and its documentation say the cost grows with the model (3.3 to
  3.8 times on the review's 5000-row model). A cheaper trigger was not
  added: the measurement under "Cost" below found none that keeps the
  199 of 200 (a |z| screen would skip little and needs a Hessian; the
  gradient screen at `b = 0` never signals).
- m6. `summary()` reports "Monotonic Simplex Parameters" as brms does,
  rows `moincome1[1]`, ..., with delta-method standard errors through
  the outer covariance and Wald intervals held in `[0, 1]`
  (`summary_mo_frame()`); `confint()` and `vcov()` keep the chart
  coordinates, documented in `?frm`. The draws are renamed as above.
- m7. The claim-7 construction is a test ("a lost direction's small
  loading stays out of a kept row", `test-optima.R`); the review saw it
  fail on the mutant without the zeroing line (null row 0.00189).
- m8, m9. Seed 194 and cox()'s softmax baseline are in
  dev/test-backlog.md; the cox change is not local (its sampler would
  need the same treatment).
- m10. A search that replaces `opt` keeps the escape's
  `stationary_escape` record.

### Reruns

- The `mo()` study on the final build (`dev/optima-mo-study.R`,
  `mo-lane4`): 199 of 200 at the exact maximum to 1e-6, max gap
  0.04888, codes all 0, no SE lost, no warning, evaluations 85701;
  paired with the first round's lane run, 0 logLiks differ. The plain
  `ls ~ mo(income)` (`mo-main-lane4`): 200 of 200, max gap 1.6e-9,
  evaluations 15167, as before.
- Seen to fail on rellib-r6 (`dev/release/optima-p1b-base-files/`,
  `optima-p1c-base-files/`): frmtmb.sample's new
  `test-mo-simplex-draws.R` (pass 0, fail 3, error 2: no `simo_`
  columns; on the weak case 26 divergences, R-hat 2.23 and the
  draws at the vertex), and in `test-optima.R` the escape-state test
  (modes 1.0531 against 1.0764) and the summary test (no `mo` block).
  The loading test passes on base, as it should: it restores coverage
  of code that exists there, and the review saw it fail on the
  mutant.
- Core and frmtmb.sample, one file per process, every gate on, lane
  library first (logs `dev/release/optima-p1-ref-files/`,
  `optima-p1-ob-files/`; `dev/optima-suite-sum.R`):

```text
reference BLAS: 265 files, 0 without a RESULT line
        frmtmb   216 17352 pass  1 fail  0 err  0 skip  0 warn
 frmtmb.sample    49  2626 pass  0 fail  0 err  1 skip  0 warn
OpenBLAS 0.3.26: 265 files
        frmtmb   216 17351 pass  2 fail  0 err  0 skip  1 warn
 frmtmb.sample    49  2626 pass  0 fail  0 err  1 skip  1 warn
```

  The reference failure is the expected ledger row
  (`brmsfit-methods:326` "now HOLDS"); the skip is frmtmb.sample's
  scale tier. Under OpenBLAS, besides that row: Malingering_2's stated
  constant, the "singular convergence (7)" of `test-ordinal-mixture.R`
  and the R-hat warning of `test-brms-shapes-draws.R`, all three the
  same on rellib-r6 under the same emulator (first round,
  `optima-obbase`). `test-perf.R` passed on both this time.

- R CMD check --as-cran, built from the final source
  (`sh dev/optima-check.sh frmtmb` and `frmtmb.sample`): frmtmb
  **Status: 1 NOTE** (the HTML manual's skipped math rendering, V8);
  frmtmb.sample **Status: OK**.

frmtmb.sample now needs the frmtmb that exports `mo_simplex()`,
`mo_coords()`, `mo_chart_frame()` and `mo_frame_terms()`; its
DESCRIPTION floor is set in punch round 2 (frmtmb (>= 0.68.1.9000)).

## Summary

| item | outcome |
|---|---|
| 1. `mo()` below the maximum | fixed: 116 of 200 seeds at the exact maximum to 1e-6 on base, 199 of 200 on the lane |
| 2. probit log-odds underflow | fixed, and the same defect fixed in cloglog and softit |
| 3. `cs()` ordinal mixtures die on a NaN gradient | fixed: 11 errors in 40 fits on base, 0 on the lane |
| 4a. vigport defect 7, `fit_loss2` | measured, not fixed |
| 4b. fixes F11 seed 8 | measured, not fixed |
| 4c. multi-start for ordinal mixtures | measured, not built |

Behavior changes (NEWS, "Breaking changes"): the `mo()` simplex
coordinates (and, since punch round 1, the sampler's `dirichlet(1)`
on a simplex and its `simo_` names), and the probit's log-odds (probit
fits move by up to 9e-10 relative in their estimates). The version
bump is a minor one for both packages. CORRECTED in punch round 2:
frmtmb.sample's floor does change, to the frmtmb that exports the
`mo_*` helpers and `summary_mo_frame()`; in the worktree core is
0.68.1.9000 and frmtmb.sample imports `frmtmb (>= 0.68.1.9000)`.

## Item 1: `mo()` fits below the maximum

### The exact maximum

`dev/optima-mo-study.R` fits brms_monotonic's own data code,
`ls ~ mo(income) * age`, seeds 1 to 200. Given the sign of each scale
coefficient the increments `b * D * w` share that sign, so the model is
least squares under sign constraints: four convex quadratic programs
(quadprog), whose best is the global maximum with no start and no
tolerance. (The test file uses the equivalent enumeration over signs
and supports, which is exact by the KKT conditions.) The script's
`|exact max a - b|` line compares the two arms' maxima: 0.

### What was wrong: two classes of miss

On base (`dev/optima-mo-sum.R dev/optima-log/mo-base2`):

<!-- generated: dev/optima-record.R, block 1 -->
```text
== dev/optima-log/mo-base2 : 200 seeds, distinct 200 range 1 200
  gap = exact max - logLik: <= 1e-6 116 ; > 1e-6 84 ; > 1e-3 61 ; > 1e-2 53 ; > 0.5 17 ; max 1.948 ; min 3.78e-11
  optimizer code: 0=183 1=17
  misses (> 1e-2) in the exact maximum's sign pattern: 21 ; in another pattern: 32
  min simplex weight < 1e-6 at the fit: 162 ; any SE not finite: 18 ; any warning: 43
  objective evaluations: total 36191 median 171 max 317
```

- 21 misses are in the right sign pattern: the softmax plateau. A step
  of 0 is at an infinite softmax coordinate, the gradient toward a
  vertex is of the order of the weight left, and nlminb stops on it.
- 32 misses are in another sign pattern. For a fixed sign of `b` the
  increments range over a convex cone, and a Gaussian (or any
  likelihood concave in the linear predictor) has one maximum per
  cone; the interaction's coefficient is near 0 on these data (no true
  interaction, `age` uncentered), so the two cones of its sign hold two
  local maxima and the fit lands in either. nanse's "38 at another
  local maximum after the escape" are these.

### Alternatives measured

Each arm is the same 200 seeds, counted against the exact maximum.
Prototypes: `dev/optima-mo-flip.R` (post-fit refits on the fit's own
tape; summaries `dev/optima-mo-flip-sum.R`). They ran on an
intermediate build whose `mo_simplex()` chose the softmax or the
sphere chart by an environment variable, `FRMTMB_MO_PARAM`, since
removed, so those arms do not rerun on the final build; the final
build's own arms are `mo-lanens` (chart alone, by
`OPTIMA_NOSEARCH=1`) and `mo-lane3`. The softmax control arm of the
intermediate build reproduced base to the bit (`mo-psoft` against
`mo-base`: 0 of 200 logLiks differ).

<!-- generated: dev/optima-record.R, prototype arms -->
```text
== dev/optima-log/flip-soft : 200 seeds, distinct 200
  gap_fit    <= 1e-6 116 | > 1e-6  84 | > 1e-2  53 | > 0.5  17 | max 1.948
  gap_flip   <= 1e-6 169 | > 1e-6  31 | > 1e-2  16 | > 0.5   5 | max 0.9501
  gap_vertex <= 1e-6 126 | > 1e-6  74 | > 1e-2  39 | > 0.5  15 | max 1.948
  gap_both   <= 1e-6 169 | > 1e-6  31 | > 1e-2  15 | > 0.5   5 | max 0.9501
== dev/optima-log/flip-sphere : 200 seeds, distinct 200
  gap_fit    <= 1e-6 174 | > 1e-6  26 | > 1e-2  23 | > 0.5   6 | max 1.157
  gap_flip   <= 1e-6 181 | > 1e-6  19 | > 1e-2  17 | > 0.5   5 | max 1.157
== dev/optima-log/cone-sphere : 200 seeds, distinct 200
  gap_flip   <= 1e-6 190 | > 1e-6  10 | > 1e-2   8 | > 0.5   2 | max 0.775
== dev/optima-log/conescreen-sphere : 200 seeds, distinct 200
  gap_flip   <= 1e-6 196 | > 1e-6   4 | > 1e-2   3 | > 0.5   0 | max 0.3691
```

| arm | at the maximum to 1e-6 |
|---|---|
| base (softmax) | 116 |
| softmax + nanse's vertex escape (`gap_vertex`) | 126 |
| softmax + sign flip from the barycenter | 169 |
| softmax + both | 169 |
| sphere chart alone | 174 |
| sphere + sign flip from the barycenter | 181 |
| sphere + bounded search of the other sign's cone | 190 |
| + start at the vertex the gradient at `b = 0` favors | 196 |
| + back from the chart's pole (the lane) | 199 |

The vertex escape is the weakest: it only treats the first class and,
in softmax coordinates, lands on the next plateau. A reparameterization
treats the first class at its cause; the second class needs a search,
since no local method crosses from one cone to the other.

### The change

`R/objective.R`, `mo_simplex()`: the simplex is the squared
coordinates of a point of the unit sphere, and the `D - 1` free
coordinates are those of stereographic projection from the pole
`-(1, ..., 1) / sqrt(D)`. `zeta = 0` is the uniform simplex, as the
softmax's zero was, so every start is where it was. Every face
`w_j = 0` is at the finite point `u_j = 0`, where the likelihood is
quadratic in `u_j` with curvature twice the derivative along the
face's normal: a face that is the maximum is an ordinary minimum of
the objective with a finite Hessian, and a face the likelihood should
leave is a saddle. `mo_coords()` is the inverse, on the sheet with
every `u_j >= 0`, which lies inside the unit ball and is the
minimum-norm preimage. `R/predict.R`'s `mo_col_values()` reads the
same map.

`R/fit.R`, `mo_search()`, run after `escape_stationary()` in
`fit_assembled()`:

1. A term whose coordinates left the unit ball is moved to
   `mo_coords()` of its own simplex and refitted. The chart's one
   point at infinity is the pole, whose image is the barycenter; with
   `b` near 0 the simplex barely moves the likelihood, nlminb takes a
   long step, and creeps toward the pole on a gradient that decays
   with the distance (seed 38 of the chart-only arm ended at
   `zeta2 = (-4601, 3066)`, 0.23 below the maximum,
   `dev/optima-mo-seed.R`).
2. A term whose simplex has a weight below `mo_face_tol` (1e-6) is
   refitted with its coefficient's sign held at the other sign by a
   bound, the simplex started at 0.9 of the vertex along which the
   gradient in `b` at `b = 0` favors that sign most, then released
   and refitted without the bound if the held fit is better. An
   interior simplex with `b != 0` is a stationary point of the
   unconstrained increments, which a concave likelihood has once, so
   it is not searched.

The objective's state is restored afterwards (`obj_state_save()`), or
set to the kept optimum: `parList()` reads the random effects' inner
modes from the last point evaluated, and the first version of the
search left a `(1 | g)` fit's modes at another point (found by
`test-brms-likelihood.R`'s check C row 3, joint density off by 0.31;
pinned by `test-optima.R`, "reports the modes of its optimum").

### Result

<!-- generated: dev/optima-record.R, blocks 1 and 2 -->
```text
== dev/optima-log/mo-lane3 : 200 seeds, distinct 200 range 1 200
  gap = exact max - logLik: <= 1e-6 199 ; > 1e-6 1 ; > 1e-3 1 ; > 1e-2 1 ; > 0.5 0 ; max 0.04888 ; min 1.137e-13
  optimizer code: 0=200
  misses (> 1e-2) in the exact maximum's sign pattern: 0 ; in another pattern: 1
  min simplex weight < 1e-6 at the fit: 192 ; any SE not finite: 0 ; any warning: 0
  objective evaluations: total 85701 median 337 max 2357
== paired on 200 seeds (base -> lane)
  logLik b - a: improved > 1e-6 84 ; worse < -1e-6 1 ; max gain 1.948 ; max loss -0.04888
  seeds at the maximum in a (gap <= 1e-6): 116 ; of those, |logLik b - a| max 0.04888
== dev/optima-log/mo-lanens (chart alone) : 200 seeds
  gap = exact max - logLik: <= 1e-6 174 ; > 1e-6 26 ; > 1e-3 25 ; > 1e-2 23 ; > 0.5 6 ; max 1.157
  misses (> 1e-2) in the exact maximum's sign pattern: 1 ; in another pattern: 22
  objective evaluations: total 41812 median 200.5 max 313
== paired (chart alone -> lane)
  logLik b - a: improved > 1e-6 25 ; worse < -1e-6 0 ; max gain 1.157 ; max loss 0
  seeds at the maximum in a (gap <= 1e-6): 174 ; of those, |logLik b - a| max 3.632e-07
```

- 199 of 200 reach the exact maximum to 1e-6 (base 116), every
  optimizer code is 0 (base: 17 at code 1), no standard error is NaN
  (base 18 seeds), and no fit warns (base 43).
- One seed is worse than base: seed 194, 0.0489 below, at the other
  sign's maximum, which base had found. The search's vertex start did
  not reach that cone's optimum there.
- The search never lost: from the chart-only fit it improved 25 seeds
  and moved the 174 already at the maximum by at most 3.6e-7.

The plain monotonic model, `ls ~ mo(income)`, seeds 1 to 200
(`dev/optima-mo-main.R`, the same exact construction over two signs):

```text
  base  gap <= 1e-6: 200 of 200; max gap 5.02e-08; codes 0=200; searched 0; evaluations total 11467, median 54.5
  lane3 gap <= 1e-6: 200 of 200; max gap 1.62e-09; codes 0=200; searched 101; evaluations total 15167, median 74.0
```

The plain fits were at the maximum and still are: paired by seed the
lane's logLik minus base's runs from -9.9e-10 to 4.9e-8, none past 1e-6
(`dev/optima-mo-main-cmp.R`). 101 of them have a weight below 1e-6
and pay for a search that gains nothing.

### Cost

Counted objective and gradient evaluations (`fit$opt$evals`, load
independent): 36191 to 85701 on the interaction (2.37 times), 11467 to
15167 on the plain model (1.32 times). The chart alone costs 41812
(1.16 times base); the search is the rest. A cheaper trigger was
looked for and not found: on the chart-only fits
(`dev/optima-mo-zscreen.R`, `dev/optima-zscreen-sum.R`), the 23 terms
the search improved had a marginal |z| of the coefficient of at most
1.27, but the terms it did not improve had a median of 1.2, so a
|z| screen would skip little, and it needs a Hessian first. The
gradient screen at `b = 0` with the other parameters held never
signals a better cone: as a "skip unless favorable" rule it skipped
every search (`coneskip-sphere`: 0 extra evaluations, 174 of 200, the
chart alone), because the other parameters still fit the current
sign. As a choice of start vertex it is what took 190 to 196. Wall clock:
`dev/optima-mo-time.R` (see "Timing" below).

### The sampler

CORRECTED in punch round 1. The first version of this section said
the sampler was unchanged, and it was not: frmtmb.sample put no
density on a simplex, so the chart itself was the prior, and the
change of chart changed what `frm_sample()` drew. The review measured
it (`dev/optima-rev-sample1-implied.R` to `-sample5-summary.R`): a
flat softmax coordinate is `1 / prod(w)` on the simplex and a
weakly identified simplex sat at a vertex (44 divergences); a flat
chart coordinate is not integrable at the barycenter and the same
simplex sat there with sd 0 and `zeta` up to 1.4e18; on `fit1` the
chart's R-hat on `zeta` was 1.10 and 1.37 while the weights agreed
with brms. The one-chain comparison that stood here (fit1 only) could
not see the weak case. The sampler now samples the softmax under
`dirichlet(1)`; see "Punch round 1" at the top.

### brms

brms samples the simplex under `dirichlet(1)` and reports a posterior,
not an ML point, so its posterior mode is the ML point only up to the
priors on the intercept and sigma; it was not run for this item. The
ported suite's monotonic fixture `fit1` of `brmsfit-methods` (a model
with `mo(Exp)`) now converges: `brmsfit-methods:326` "now HOLDS",
recorded as "cannot transfer" because "FIXTURE 1 DOES NOT CONVERGE
(nlminb code 1, NaN standard errors)". The ledger row flips; it is not
edited here (regenerated at consolidation).

## Item 2: the robust log-odds in the far tails

### Probit

`R/links.R`, `probit$logit_eta`, now
`pnorm(eta, log.p = TRUE) - pnorm(-eta, log.p = TRUE)`.
`dev/optima-probit.R` against Rmpfr at 256 bits through RTMB's tape
(log `dev/optima-log/probit-lane.txt`):

```text
== the old form's finite range : 2001 points on [ -38 , 38 ]
  old       value: non-finite    0, max rel err 4.32e-12 | derivative: non-finite   24, max rel err 5.56e-14
  installed value: non-finite    0, max rel err 5.72e-16 | derivative: non-finite    0, max rel err 9.28e-14
== wide : 4001 points on [ -200 , 200 ]
  old       value: non-finite 3232 | derivative: non-finite 3250
  installed value: non-finite    0, max rel err 6.94e-16 | derivative: non-finite    0, max rel err 3.58e-12
== far points : 14 points on [ -10000 , 10000 ]
  old       value: non-finite   14 | derivative: non-finite   14
  installed value: non-finite    0, max rel err 4.61e-08 | derivative: non-finite    0, max rel err 1e-06
old vs new on [-38, 38]: 483 of 2001 values differ, max relative 4.32e-12
```

The 483 differing values are the old form's error, not the new one's:
the old form was 4.3e-12 off where it was finite, the new one 5.7e-16,
and the old form's derivative was already NaN on 24 of 2001 points
inside `[-38, 38]`.

RTMB's `pnorm(log.p = TRUE)` derivative itself degrades smoothly from
about `|x| = 1e3`; the review measured the curve against the
Mills-ratio series (`dev/optima-rev-links.R`): relative error 4.8e-11
at 1e3, 1.3e-9 at 1e4, 2.7e-8 at 3e4, 3.2e-7 at 1e5, 2.0e-5 at 1e6,
1.5e-3 at 9.4e6, 0.34 at 1e8, and exactly 0 (relative error -1) from
1e9 on, with the derivative non-finite on 137 of 2001 points within
0.1 percent of 4.46e9 (this lane's coarser scan counted 270 of 2001
over a wider window, `dev/optima-pnorm-scan.R`). The value is exact
to 1e154. A degenerate component of
`mixture(cumulative("probit"), acat())` on the review's data reached
4.46e9 and died on that NaN. So past `|eta| = 1000` the log-odds is
read at +-1000 and continued as the tail's leading term
`+-(eta^2 - 1000^2) / 2` (the "far points" row: 4.6e-8 in value, 1e-6
in derivative, at 1e4). The clamps are exact inside `[-1000, 1000]`.
Upstream: RTMB's `pnorm(log.p = TRUE)` derivative in the far tail, for
the user to file.

### Cloglog and softit (found on the way)

With the probit fixed, the review's
`mixture(cumulative("probit"), sratio("cloglog"), acat())`
(`test-nonfinite-gradient.R`) ran on into sratio's cloglog far tail
and died with "NA/NaN gradient evaluation" (`dev/optima-mix3.R`):
its log-odds `log(-expm1(-t)) + t`, `t = exp(eta)`, has a derivative of
`1 / t` times `t` on the tape, and `1 / t` overflows once `t` is
subnormal. The softit's `log(log1p(exp(eta)))` has the same defect.
Below `eta = -40` both log-odds are `eta` to double precision (the
rest is under `exp(eta) / 2 < 2.2e-18`), and that is what they are now,
read at -40 with slope one. `dev/optima-cloglog.R` against Rmpfr
(log `dev/optima-log/cloglog.txt`):

```text
cloglog  eta in [-745, -708]:  old derivative non-finite 952 of 1001; new 0, max rel err 0
         eta in [-1000, -745]: old value non-finite 1000, derivative 1001; new 0, max rel err 0
         eta in [-40, 6.5]:    old vs new values differ on 0 of 1001
softit   the same counts; eta in [-40, 30]: old vs new values differ on 0 of 1001
```

## Item 3: `cs()` on a cumulative component of an ordinal mixture

`dev/optima-csmix.R` repeats the release review's construction
(`dev/relrev-csfail.R`, 20 seeds, n = 300); base reproduces its
counts exactly.

<!-- generated: dev/optima-record.R, cs() block -->
```text
  arm base
    cum                    ok=20
    sratio                 ok=20
    mix_cum_cum            nan_grad=4 nonconv=16
    mix_cum_sratio         nan_grad=7 nonconv=11 ok=1 warn=1
    mix_sratio_sratio      nonconv=6 ok=3 warn=11
    mix_cum_sratio_nocs    nonconv=5 ok=14 warn=1
  arm lane6
    cum                    ok=20
    sratio                 ok=20
    mix_cum_cum            nonconv=20
    mix_cum_sratio         nonconv=18 ok=1 warn=1
    mix_sratio_sratio      nonconv=6 ok=3 warn=11
    mix_cum_sratio_nocs    nonconv=5 ok=14 warn=1
  fits with a logLik on both arms: 109 ; identical: 109
  controls (no cs() on a cumulative component of a mixture): 80 fits; logLik identical 80
```

### Causes

1. **nlminb's reported point** (`dev/optima-csmix-debug.R`, seed 3):
   43 of 89 trial points of the first run had a NaN objective (a row's
   `cs()` thresholds crossed in its own category). PORT stopped with
   "false convergence (8)" and returned `par` at the last, rejected,
   trial: the objective there is NaN while `objective` is 391.683 from
   another point. frm()'s restart began at that `par`, its first
   gradient was NaN, and the fit died. This is general, not
   mixture-specific: any fit whose run ended on a rejected trial read
   its estimates from it. `dev/optima-parcheck.R` counts the fits of
   the four mixture designs whose estimates sit at a non-finite
   objective: 3 on base (1 `mix_cum_sratio`, 2 `mix_cum_cum`, besides
   their 11 errors), 0 on the lane. (6 fits on each arm differ from
   their reported objective by 5.7e-14 to 3.4e-13, which is rounding.)
2. **The touch** (seed 12, `dev/optima-csmix-where.R`): after the first
   fix, nlminb's line search, backing off the crossed side, converged
   onto a point where one row's two thresholds are the same double.
   The component's log-density there is `-Inf`, the mixture's stays
   finite through the other component, and its gradient is NaN.

### The change

- `R/fit.R`: `nlminb_trial_fn()` remembers the last and the best
  finite point; `nlminb_best_par()` replaces `par` by the best point
  when `par` is the rejected last trial, and `run_optimizer()` then
  evaluates the objective there so its state (and the inner modes)
  follow.
- `R/families.R`: in a mixture (`gap_floor`),
  `ord_cumulative_logpmf_cs()` adds `ord_touch_nan(a - b)`, which is 0
  with a zero derivative for a positive gap and NaN at a zero gap, so
  the touch is as undefined as the crossing past it. It is spelled
  `log(g + abs(g)) - log(2 * abs(g))` because the tape folds `0 * x`
  and `x - x` (`dev/optima-guard-check.R`), and with `abs()` so plain
  numeric evaluation warns about nothing (a first spelling with
  `log(2 * g)` let "NaNs produced" escape on 2 fits).

**A floor was tried and rejected.** Holding a crossed row's own category
at `ord_log_interior()`'s floor made the objective finite everywhere,
and raised `mix_cum_sratio`'s log-likelihood by up to 115 units (seed
1: -394.68 to -279.44): a crossed component's OTHER categories read
probabilities summing past one, and the floor let the optimizer use
that (`csmix-lane3`, with "a step function of its predictors" warnings
on every fit). Crossing has to stay undefined. The floor also exposed
that `ord_log_interior()`'s minimum form loses its cap when `b - a` is
positive (rtmb-pitfalls item 11); the cs path now never passes it a
crossed pair, and its documentation says so.

### What remains: the wall

The non-convergence is honest. `dev/optima-csmix-wall.R` on the lane:

```text
code 1 fits: 18 ; of them with min gap < 1e-3: 18
code 0 fits: 2 ; min gaps: 0.509 0.893
```

All 18 fits that stop with "false convergence (8)" have a row whose two
thresholds in its own category are within 2.3e-13 of each other, with
max |gradient| 0.15 to 9.3: the maximum of these data lies on the
boundary of the region where every row's category is open, which a
barrier-style objective reaches but cannot certify. brms has the same
geometry (brms-22, crossing `cs()` thresholds).

## Item 4

### 4a. vigport defect 7, `fit_loss2`

`dev/optima-loss2.R` (the ClarkTriangle data, cached in
`dev/optima-log/`; 40 starts jittered about the default's estimates):

```text
default: logLik -358.545112 code 1 max|grad| 0.00311 (false convergence (8) | warned) evals 267
logLik over starts: max -355.4707197 ; code 0 max -358.5132284
```

13 of 40 starts end at code 0, 10 of them at -358.5132 and 3 between
-358.5232 and -358.5178; the default is 0.0319 below the best and says
so ("false convergence"). One start reached -355.47
at a point with max |gradient| 9.99e7, not an optimum. `restarts = 3`
and `restarts = 10` stay at -358.545112; `optimizer = "optim"` ends at
-358.7135 with code 0 (`dev/optima-loss2-restarts.R`). The lane's
changes leave the default identical. Not fixed: the default does warn,
the better optimum is 0.03 log-likelihood units away, and no cheap
change reached it.

### 4b. fixes F11 seed 8

`dev/optima-f11.R` on base and lane alike: default -589.023726 with
sigma's smoothing log-SD at -0.80; from that SD started at -10 it ends
at -588.881888 (49 evaluations); 4 of 20 jittered starts find it, 0.1418
above the default. A one-evaluation screen (each smoothing SD moved to
-10 with the rest held, `dev/optima-f11-screen.R`, seeds 1 to 40) does
not see it: on seed 8 it fires on the other SD and gains 8e-8, and on
38 of 40 seeds it fires with gains of 1e-9 to 1.3e-5. A refit from
the boundary per smoothing SD would cost one refit per SD on every
smooth fit for a defect fixes measured on about 1 percent of fits. Not
fixed.

### 4c. multi-start for ordinal mixtures

`dev/optima-ordmix-starts.R` (20 seeds each; 10 starts jittered by
N(0, 0.5) about the default start):

```text
cum_cum    seeds 20; best of 10 above the default by > 1e-3: 9; of those degenerate: 5; non-degenerate: 4; default degenerate: 6
cum_sratio seeds 20; best of 10 above the default by > 1e-3: 11; of those degenerate: 7; non-degenerate: 4; default degenerate: 6
```

Best-of-10 beats the default on 20 of 40 fits, but on 12 of those 20
the better point is degenerate (the component warning fires), which
confirms ordmix's 16 of 22. 8 of 40 fits have a better non-degenerate
optimum. An option that ranks by likelihood picks the degenerate point
more often than not, and one that ranks by anything else is not
maximum likelihood. Not built; `start =` and `frm_allfit()` reach the
same points on request.

## Tests

New: `tests/testthat/test-optima.R`; `tests/testthat/helper-mo-flat.R`.
Seen to fail on rellib-r6 (`dev/release/optima-base-t2-files/`):

```text
test-optima.R: a mo() fit reaches the exact maximum: relative gap 0.006002 (seed 12), 0.004658 (54), 0.000857 (38)
  a cs() ordinal mixture at its crossing wall warns, not errors: "NA/NaN gradient evaluation"
  a cumulative probit stays defined past |eta| = 38.2: objective and gradient not finite
  the cloglog and softit log-odds keep a derivative far below: -Inf values, NaN derivatives
  the chart and search tests: mo_simplex / mo_search absent (the weak form)
test-se-check.R: "a mo() simplex on its boundary keeps a finite Hessian": solve(H) fails
test-nonfinite-gradient.R: "the three-component mixture ... ends finite": gradient not finite
```

The nlminb test is the cs() mixture of seed 3 on its own tape; its
first spelling, a toy objective with a NaN wall, passed on base and was
replaced.

Changed, because they read the softmax plateau as a fixture:
`test-se-check.R` (seeds 7, 71 and 12 of brms_monotonic's data; 71 and
12 are now `mo_flat_data()`, a monotonic predictor whose middle
category is never observed, so its one coordinate does not enter the
likelihood and its Hessian row is exactly zero on every seed),
`test-diagnostics-ux.R` (seed 71, the same replacement),
`test-numerical-robustness.R` (pinned the probit underflow at 40),
`test-nonfinite-gradient.R` (its fixture was the probit underflow);
`helper-brms.R`, `test-brms-agreement.R`, `test-mo-terms.R` read the
simplex through `mo_simplex()` instead of a softmax.

Lost coverage: seed 12's concave direction with a 0.022 loading on a
kept parameter cannot be built from a mo() fit any more; the test that
pinned it ("a kept parameter's small loading") now uses two exactly
flat coordinates and checks the bands that do not read them. A concave
direction with a small loading needs another construction (for lane
setier).

### All eight suites, both BLAS builds

One file per process, the lane's core first in the library path (every
log's `lib:` line names `wt-optima-lib`), `NOT_CRAN`,
`FRMTMB_BRMS_FIT_TESTS`, and for the OpenBLAS pass also
`FRMTMB_DRMTMB_FIT_TESTS` and `FRMTMB_FUZZ` (the reference pass ran
those two files again with the gates on: drmTMB agreement 131 pass,
fuzz 2 pass, 0 fail). Summaries: `dev/optima-suite-sum.R optima-ref`
and `optima-ob` (logs `dev/release/optima-ref-files/`,
`dev/release/optima-ob-files/`).

```text
reference BLAS (optima-ref): 355 files, 0 without a RESULT line
             pkg files  pass fail err skip warn
          frmtmb   216 17206    2   0   14    0
 frmtmb.coupling    11   542    0   0    5    0
      frmtmb.eam    29  1743    0   0    3    0
   frmtmb.latent    10   360    0   0    2    0
    frmtmb.learn    15   501    0   0    2    0
      frmtmb.ode    11   549    0   0    1    0
   frmtmb.sample    48  2616    0   0    1    0
   frmtmb.spline    15   590    0   0    1    0
OpenBLAS 0.3.26 (optima-ob, dev/optima-openblas.sh): 355 files
          frmtmb   216 17337    4   0    0    1
 frmtmb.coupling    11   542    0   0    5    0
      frmtmb.eam    29  1743    0   0    3    2
   frmtmb.latent    10   360    0   0    2    0
    frmtmb.learn    15   501    0   0    2    0
      frmtmb.ode    11   549    0   0    1    0
   frmtmb.sample    48  2616    0   0    1    1
   frmtmb.spline    15   590    0   0    1    0
```

Every extension skip is its `test-scale.R` (the scale tier, gated by
`FRMTMB_SCALE_TESTS`, not run). The failures and escaped warnings:

- `test-brms-suite-methods.R`, both BLAS: `brmsfit-methods:326` "now
  HOLDS", the expected stale ledger row (item 1, "brms").
- `test-perf.R`, both BLAS: the wall-clock bound the backlog already
  files ("Filed at the 0.68.0 release", low). Alone: lane 1 of 2 runs
  failed (1.72 s against 1.00), OpenBLAS 1 of 1 (1.69 s), and base
  rellib-r6 1 of 2 (2.23 s against 2.00) under the same load.
- `test-optima.R`, OpenBLAS: the new test asserted a positive search
  gain on seed 54, which is 3.6e-11 with the reference BLAS and 0 with
  OpenBLAS (`dev/optima-seed54.R`). The assertion is now `>= 0`; the
  file passes 79 of 79 on both BLAS (`optima-ref-t3`, `optima-ob-t3`).
- OpenBLAS only, all pre-existing: run on rellib-r6 with the same
  emulator (`optima-obbase`), `test-bcm-latent-mixtures.R` fails the
  same assertion with the same number (Malingering_2's stated constant,
  measured -0.0001295824987), `test-ordinal-mixture.R` lets the same
  "singular convergence (7)" escape (the backlog's item at :751),
  frmtmb.eam's `test-sampling.R` the same two ESS warnings and
  frmtmb.sample's `test-brms-shapes-draws.R` the same R-hat warning.

### Timing

`dev/optima-mo-time.R`, 50 fits each, three interleaved rounds per
arm, on a machine shared with other lanes (log
`dev/optima-log/mo-time.txt`):

```text
TIME arm base interaction 1.81 s main 1.84 s control 1.81 1.57 s
TIME arm lane interaction 2.58 s main 1.47 s control 2.30 2.05 s
TIME arm base interaction 1.45 s main 1.45 s control 1.94 1.89 s
TIME arm lane interaction 1.33 s main 1.61 s control 1.65 1.75 s
TIME arm base interaction 1.94 s main 1.82 s control 1.54 2.41 s
TIME arm lane interaction 1.53 s main 0.86 s control 2.34 2.46 s
```

The control moves by 1.6 times between processes, more than any arm
difference, so the wall clock cannot separate the arms at this size:
a fit of this model is about 30 ms, most of it per-fit setup rather
than objective evaluations. The evaluation counts above are the
measure.

### R CMD check

`sh dev/optima-check.sh frmtmb` (`R CMD check --as-cran`, pandoc and
TinyTeX on PATH, the remote incoming check off, `R_LIBS` the lane's
library, rellib-r6, then the user library): **Status: 2 NOTEs**, both
the expected ones: examples over 5 s on a loaded machine
(`frm_hazard_reads` 14.84 s elapsed for 1.64 s of CPU, `dharma_residuals`,
`VarCorr`; dev/lane-rules.md: such a NOTE measures load) and the HTML
manual's skipped math rendering (V8). Only core changed, so no
extension was checked. The tarball was built before three later edits:
two comments in `R/fit.R` (`mo_face_tol`'s measured weights and a
clause of `mo_search()`'s), and `test-optima.R`'s seed-54 assertion
(`> 0` to `>= 0`, which passes either way with the reference BLAS the
check used). The installed library carries the final source.

## Defects found, not fixed

- RTMB `pnorm(log.p = TRUE)`: derivative 1.2 percent off at 9.4e6 and
  NaN near 4.46e9 (upstream, `dev/optima-pnorm-scan.R`).
- FIXED in punch round 1 (m4): `escape_stationary()` kept the best of
  its runs but not the objective's state at it. The first version of
  this item called it unreached; the review found it on an ordinary
  `skew_normal()` with `(1 | g)`: 5 of the 14 escapes of seeds 1 to
  40 reported modes 0.002 to 0.143 away from the optimum's.
- The SE check misses an exactly flat simplex direction that is not
  along a coordinate: with four categories and the third unobserved
  (`dev/optima-segap.R`), the flat direction warned on 6 of 9 fits and
  was silent on 3 (tier 1's knife edge, backlog "Filed at 0.68.1",
  lane setier).
- Seed 194 of the interaction study: the search misses the other
  cone's optimum (0.0489 below; base found it).
- frmtmb.sample has no simplex prior (improper posterior in either
  coordinates).
- The cs() ordinal-mixture wall: no message says that a row's category
  closed; the convergence warnings are right but do not name the cause.

## Files

R: `R/objective.R` (`mo_simplex()`, `mo_sphere_basis()`, `mo_coords()`
new; `lp_eta_fixed()`), `R/predict.R` (`mo_col_values()`), `R/fit.R`
(`nlminb_trial_fn()`, `nlminb_best_par()` new, `run_optimizer()`,
`fit_assembled()`, `mo_search_terms()`, `mo_search()`, `mo_face_tol`
new, the `frm()` Rd section "Monotonic effects"), `R/families.R`
(`ord_cumulative_logpmf()`, `ord_cumulative_logpmf_cs()`,
`ord_touch_nan()` new, `ord_log_interior()` documentation only),
`R/links.R` (`probit`, `cloglog`, `softit` `logit_eta`), `R/priors.R`
(the `simo` message), `R/confint.R` (a comment), `man/frm.Rd`,
`NEWS.md`.

Punch round 1 added: `R/objective.R` (`mo_simplex()` and `mo_coords()`
take `chart`; `mo_chart_frame()`, `mo_frame_terms()` new),
`R/sampling-api.R` (the four exported, documented), `R/fit.R`
(`obj_settle()` new; `escape_stationary()`, `mo_search()`,
`fit_assembled()`, `frmtmb_control()` with `mo_search`, the
`check_se` and "Monotonic effects" documentation), `R/methods-fit.R`
(`summary_mo_frame()` new, `summary.frmtmb_fit()`,
`print.summary.frmtmb_fit()`), `R/se-check.R` (one warning clause and
three comments), `R/links.R` (comments), `vignettes/diagnostics.Rmd`,
`man/frm.Rd`, `man/frmtmb_control.Rd`, `man/frmtmb-sampling-api.Rd`,
`NAMESPACE`, `NEWS.md`; frmtmb.sample: `R/sample.R`
(`prior_augmented_obj()`, `ncp_objective()`, `frm_sample()`,
`default_prior_notes()`, the default-prior documentation;
`mo_sample_frame()`, `mo_sample_pars()`, `mo_simplex_nlp()`,
`mo_draws_to_chart()`, `draws_mo_cols()` new), `R/draws-brms.R`
(`draws_natural_cols()`), `man/frm_sample.Rd`, `NEWS.md`; tests
`test-optima.R` (four tests added), `test-brms-likelihood.R`,
`test-brms-agreement.R`, `test-brms-names.R` (comments), frmtmb.sample
`tests/testthat/test-mo-simplex-draws.R` (new); records
`dev/test-backlog.md`, `dev/upstream-bugs.md` (RTMB-3).
