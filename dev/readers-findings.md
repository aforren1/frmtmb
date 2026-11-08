# Lane readers: prediction readers that read the wrong quantity

Round opened 2026-11-06, base 3cd2a0c8 (frmtmb 0.79.1). Worktree
`frmtmb-wt-readers`, library `wt-readers-lib`, base library `rellib-r18`.
Status: DONE, uncommitted. This file is the record if the report does not
arrive.

## Shared functions this lane changes

- `R/me.R`: `me_row_loads()` (rewritten), `me_node_value()` (rewritten),
  `me_band_end()`, `me_check_new_levels()`; new `me_src_add()`,
  `me_reads_groups()`, `psd_root()`, `chol_rows()`, `me_gh_moments()`,
  `me_band_dir()`, `me_warn_not_finite()`, `me_gh_size`,
  `me_gh_slope_size`.
- `R/predict.R`: `frm_linpred()` (the body branch), `lp_basis_layer()`,
  `lp_extra_var()`, `predict_mean_se()` (the body branch),
  `prob_band_draw()`, `fitted_me_moments()`, `lp_basis_nl()`; new
  `nl_basis_extra()`, `nl_extra_var()`; roxygen of `fitted()`.
- `R/conditional-effects.R`: the `se__` lines of the mean and dpar
  displays, the category band, `ce_pool_scale()`, new
  `ce_cat_draw_band()`, roxygen of `conditional_effects()`.
- `R/compat.R`: three rows (below).
- `extensions/frmtmb.sample/R/methods-draws.R`:
  `draws_refuse_new_levels()` and its three predicting callers.

`lp_basis_nl()` holds `<<-` sites that lane superassign rewrites. This lane
deletes the refusal block inside `build()` and adds lines after it; it does
not touch a `<<-` line.

## Scripts and logs

Every script sources `dev/readers-setup.R`: `READERS_ARM=before` reads
`rellib-r18` alone, anything else puts `wt-readers-lib` first. Logs are in
`dev/readers-out/` (gitignored). Tests run through `dev/readers-run-tests.R`
and `dev/readers-run-par.sh`, copies of `dev/release/run-tests.R` and
`run-par.sh` on `rellib-r18` with `FRMTMB_RELLIB` set to the lane library
(the release copies name `rellib-r17`).

## The cause shared by items 1, 4 and 7

`fitted()` on a new row integrates the expected response over the draws a
new row reads (`fitted_me_moments()`). Through 0.79.1 that integral had one
dimension: `me_row_loads()` took the first predictor that read a `me()`
value as the lead, moved every other predictor by its regression on the
lead, and, when a latent residual sat on `mu` while the `me()` value sat on
another predictor, dropped the `me()` value. A new grouping level's effects
were not a draw at all: they sat at 0, the median.

The fix treats every draw of a new row as one normal source over the
predictors it moves: each `me()` call (one draw whichever predictors read
it), the product part of a `me()` interaction, a latent residual on `mu`,
and each column of a root of a new level's block covariance (one draw per
block and level, shared by every predictor that holds a component). Per row
the covariance over the predictors is factored with the lead first
(`chol_rows()`), and the mean and variance are a product Gauss-Hermite rule
over the factor's directions: 41 nodes for one direction, as before, then
21, 11 and 7. A node is one `frm_linpred()` call, about 1 ms
(`dev/readers-cost.R`: 0.95 ms per call at 3 and at 100 rows). The band
reads the draw's quantiles along the direction in which the value moves
most, widened to the whole variance (`me_band_dir()`); with one direction
this is the 0.79.1 band.

## Items

### 1. `fitted()` on a new random-effect level

Script `dev/readers-item1.R`, data seed 4101 (poisson `(1 | g)`), 4102
(correlated slope), 4103 (bernoulli); Monte Carlo seed in the script, 2e6
draws of the effect at the estimates. Logs `item1-before.txt`,
`item1-after.txt`.

| case, row | 0.79.1 | lane | reference |
|---|---|---|---|
| poisson, sd 0.999, x = -0.5 | 1.512779 | 2.491786 | closed form 2.49179, MC 2.49162 |
| same, ratio closed form / Estimate | 1.64716 | 1 | |
| same, Est.Error | 1.532480 | 3.288065 | MC sd over the effect alone 3.27867 |
| `(1 + x | g)`, x = 0.8 | 2.732690 | 4.667388 | closed form 4.66739 |
| bernoulli, sd 1.593, x = -0.5 | 0.4440574 | 0.4614946 | MC 0.461261 |
| same, Est.Error | 0.4004766 | 0.2836988 | MC sd over the effect 0.278832 |

0.79.1 read the median `exp(eta)`. The bernoulli `Est.Error` was too wide
(the delta method through a logit at a variance of 2.5), and the poisson
one too narrow.

Done: `fitted_me_moments()` and `me_row_loads()` take a new level as a
draw. `lp_extra_var()` drops the new level's variance when the fit carries
`me_no_extra` (`me_coef_only()`), so that the coefficient part of
`Est.Error` and the integral's variance are not counted twice.
`scale = "linear"`, `re_formula = NA` and `sample_new_levels =
"old_levels"` do not change (tested).

### 2. `conditional_effects()`'s band for a category under `me()`

Script `dev/readers-item2.R`, data seed 7971 (the design of
`dev/review-predread-r1-check1.R`), row `xo = 1`, `sx = 1.2`; the reference
is the draw's quantiles at the estimates, 2e4 draws, seed 11. Logs
`item2-before.txt`, `item2-after.txt`.

| category | `fitted()` Q, both builds | CE band 0.79.1 | CE band lane | draw at the estimate |
|---|---|---|---|---|
| cumulative 2 | 0.047 to 0.250 | 0.068 to 0.415 | 0.047 to 0.250 | 0.048 to 0.245 |
| cumulative 3 | 0.090 to 0.254 | 0.196 to 0.305 | 0.090 to 0.254 | 0.098 to 0.247 |
| categorical 2 | 0.170 to 0.356 | 0.209 to 0.426 | 0.170 to 0.356 | 0.190 to 0.345 |

Cause: the category display builds a Wald interval on the logit about the
value at the held `x` with `se__`, which holds the draw. A middle category
has a hump in the draw, which no interval about one point reads. Done: the
display takes `fitted()`'s `Q` columns on rows that hold a draw
(`ce_cat_draw_band()`); a row with none keeps the Wald band (tested).

### 3. `se__` on the link scale

Script `dev/readers-probe.R`, seed 7971, `x` at -1, 0, 1. Log
`probe-before.txt`. 0.79.1 at `x = 0`: poisson `se__` 0.04538 against
`fitted()`'s `Est.Error` 0.06473; lognormal 0.02666 against 0.04503; hurdle
0.04703 against 0.05425. brms's `se__` is the MAD of its draws of the
displayed value (`conditional_effects.R` of brms 2.23.0, lines 462 to 466:
`posterior_summary(out, robust = robust)`, robust `TRUE` by default). Done:
`se__` is the response-scale delta-method standard error on the mean and
the dpar displays; the band is unchanged. `ce_pool_scale()` now says
`se_on = "response"` and names the band's own link (the log on a lognormal
mean), so `conditional_effects()` of a `frm_multiple()` result pools on the
band's scale; through 0.79.1 a lognormal's log-scale `se__` was pooled as if
on the identity link.

Against brms 2.23.0, run (`dev/readers-item3-brms.R`, data seed 7971,
brms 4 chains of 1000 draws, seed 3; logs `item3-brms-before.txt`,
`item3-brms-after.txt`), `se__` over brms's `se__` at x = -1, 0, 1:

| fit | 0.79.1 | lane |
|---|---|---|
| poisson mean | 1.493, 0.694, 0.329 | 1.027, 1.010, 1.012 |
| lognormal mean, `sigma ~ x` | 1.302, 0.561, 0.223 | 1.049, 1.015, 0.972 |

### 4 and 7b. A latent on `mu` and a `me()` value on another predictor

Script `dev/readers-probe2.R`, data seeds 9201 and 9202, Monte Carlo seed
11, 2e6 draws of each source at the estimates. Logs `probe2-before.txt`,
`probe2-after.txt`. `hurdle_poisson()`, `me(zo, sz)` on `mu` and `hu`, a
latent on `mu`:

| row | 0.79.1 mean, draw sd | lane mean, draw sd | MC mean, sd |
|---|---|---|---|
| 1 | 1.36224, 0.365888 | 1.35583, 0.408586 | 1.35598, 0.408559 |
| 2 | 2.11044, 0.732977 | 2.14134, 0.925262 | 2.14175, 0.925734 |

Item 4 as filed, `me()` on `hu` alone beside a latent on `mu`
(`dev/readers-item4.R`, data seed 9204, Monte Carlo seed 11, 2e6 draws;
logs `item4-before.txt`, `item4-after.txt`):

| row | 0.79.1 mean, draw sd | lane mean, draw sd | MC mean, sd |
|---|---|---|---|
| 1 | 1.777451, 0.558723 | 1.739476, 0.574969 | 1.739612, 0.575456 |
| 2 | 0.468529, 0.179825 | 0.522019, 0.358066 | 0.521947, 0.358196 |

0.79.1 dropped the `me()` draw on `hu`: row 2 was 10 percent low, its sd
half.

The lognormal arms of both items read `sigma` through a draw, and their
mean is infinite: item 7a below.

### 7a. lognormal, `me()` on `mu` and `sigma`

The mean `exp(mu + sigma^2 / 2)` with `log(sigma)` normal has no finite
integral: `E[exp(c exp(2 b z))]` diverges for any `c > 0` and `b != 0`.
`dev/readers-probe2.R`, seed 9203, `sx = 0.8`, 2e6 Monte Carlo draws:

| row | 0.79.1 Estimate, Est.Error | MC mean | MC sd of four quarters of the draws |
|---|---|---|---|
| 1 | 1.700338e+162, Inf | 6.11522e+20 | 1.730e+24, 2.803e+16, 3.405e+19, 5.400e+19 |
| 2 | 1.117850e+218, Inf | 9.63198e+210 | Inf, 4.781e+149, Inf, Inf |

So the review's "207 against 3.56" compares two finite truncations of an
infinite integral. The first round returned `Inf` by a numeric tail test;
the review found that wrong (M1, M2) and the fix round replaced it: see
"Fix round" below, which supersedes this paragraph.

### 5. `conditional_effects()` on a body whose `nlpar` reads `me()`

Script `dev/readers-probe.R`, seed 5101. 0.79.1 stopped: "Its chain rule
through a nonlinear body has not been measured". Done: `lp_basis_nl()` gives
each linear part that holds a `me()` draw or a `gp()` kriging variance a
shift of its own on the tape, one input per row, reads `d eta / d eta_part`
off the diagonal of the Jacobian, and passes each source through at first
order, a `me()` call summed over the parts it reaches before the square
(`nl_extra_var()`). `frm_lp_basis(extra_cov = TRUE)` on such a body stops
by name: the covariance of one draw across rows through a body is not
written. Test: `a * exp(b z)`, `a ~ me(xo, sx)`, where the first-order
variance is exact, `(exp(b z) b_a sx)^2`, to 1e-8 relative.

### 6. frmtmb.sample `posterior_epred(newdata, allow_new_levels = TRUE)` at a new `me(gr = )` level

0.26.0 refused it in `draws_refuse_new_levels()`, which asked core's design
builder whether a level is new, and core refuses a new `me(gr = )` level
unless `allow_new_levels` reaches it. Done: under the flag the draws carry
`fit$me_new_levels`, which core's `me_check_new_levels()` reads, so the
`me()` level passes and each draw takes its latent value from
`normal(x, sdx)` (`me_newdata_latent()`, unchanged); a group-level term's new
level is still refused. Test `test-me-gr-new-level-draws.R` (seed 7511,
4000 draws at the estimate): the draws' mean is within 4 Monte Carlo
standard errors of `fitted()`'s integrated `Estimate`, their sd within
`4 / sqrt(2 n)` of the integral's.

### 8. A nonlinear body's `Est.Error` is `NA`

Script `dev/readers-probe.R`, seed 8101: 0.79.1 `NA` at both rows. Done:
`frm_linpred(se.fit = TRUE)` on a body reads `lp_basis_nl()`'s taped
Jacobian, and `predict_mean_se()` takes a body predictor's Jacobian as its
design, so a hurdle or lognormal mean with a body answers too. Test: the
delta method by hand to 1e-6, and the sd of the body over 20000 draws of
the coefficients from `N(estimate, V)` (seed 11) within `4 / sqrt(2 n)`
plus 2 percent of it.

## Seen to fail on 0.79.1

`tests/testthat/test-predict-draw-mean.R` on `rellib-r18`
(`newtest-before.txt`): pass 13, fail 23, error 2. Lane
(`newtest-after.txt`): pass 47, fail 0, error 0, warn 0. Every test fails
behaviorally on 0.79.1 except the unit test of `me_gh_moments()`, which is
the missing-symbol form and is pinned behaviorally by the others.
`extensions/frmtmb.sample/tests/testthat/test-me-gr-new-level-draws.R` on
`rellib-r18` (`sampletest-before.txt`): error 1 ("cannot predict a grouping
level the fit did not see"); lane (`sampletest-after.txt`): pass 5.

## Compat rows changed

- `nl` x `predict`: conditional to works, with the reason.
- `me()` x `fitted`: the sentence on several draws and the lognormal.
- `group:autocor` x `fitted`: the sentence on a second draw.

## What this lane did not do, and why

- The kriging variance of an exact `gp()` at a new position is not one of
  the draws `fitted()` integrates: it stays in `Est.Error` by the delta
  method, as in 0.79.1. brms's mean over the field's draws is the same
  kind of integral; it is filed.
- A new level of a Student-t block is integrated as a normal with the t
  variance, as `lp_extra_var()` reads it.
- A new level inside a nonlinear body: `frm_lp_basis()` still refuses
  `allow_new_levels = TRUE` on a body (`test-predict-lp-basis.R` pins
  it), and a nonlinear parameter's new level is not one of
  `me_row_loads()`'s draws. Not measured this round.
- More than four independent directions on one row: the first four are
  integrated, with a warning of class `frmtmb_draw_dims`.

## Test runs

All on `wt-readers-lib` ahead of `rellib-r18`, one file per process, three
at a time, `NOT_CRAN=true`; every log's `lib:` line names
`wt-readers-lib` for frmtmb.

- Core, ungated, all 309 files (`dev/readers-out/core1.log`, generated
  sums over its RESULT lines): pass 22558, fail 1, err 0, skip 239,
  warn 0. The fail was `test-method-residue.R`, "se.fit still refuses a
  nonlinear predictor's response mean": it pinned the refusal item 8
  removes. The test now asserts the answer (identity link equal to the
  link-scale se; a lognormal body `log(a) - b * x` equal to the linear
  `yp ~ x` to 1e-3). Skips: 222 are the brms and Stan gates, 13 drmTMB,
  one each scale and fuzz tier, two brms rows.
- Three other files pinned the same refusal and were changed with it, and
  pass on the lane: `test-nl.R` (now the delta method by hand),
  `test-nlf.R` (now equal to the linear `sigma ~ z`),
  `test-predict-lp-basis.R` (now equal to `frm_lp_basis()`'s A V A').
- `test-method-residue.R` rerun after the change: pass 67, fail 0.
- Extension suites, ungated, all 192 files (`dev/readers-out/ext1.log`;
  extensions read from `rellib-r18`, frmtmb.sample from the lane library):

  | package | files | pass | fail | err | skip | warn |
  |---|---|---|---|---|---|---|
  | frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
  | frmtmb.eam | 33 | 1848 | 0 | 0 | 5 | 0 |
  | frmtmb.latent | 12 | 683 | 0 | 0 | 2 | 0 |
  | frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
  | frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
  | frmtmb.sample | 86 | 3967 | 0 | 0 | 18 | 0 |
  | frmtmb.spatial | 8 | 478 | 0 | 0 | 6 | 0 |
  | frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |

## Cost

`dev/readers-timing.R`, minimum of 5 rounds per arm, each arm in its own
process (logs `timing-before.txt`, `timing-after.txt`). The control,
`frm_linpred()` at three rows, read 1.04 ms on 0.79.1 and 1.22 ms on the
lane, so the machine's load moved by about 17 percent between the arms.

| call | 0.79.1 | lane |
|---|---|---|
| `fitted()`, three rows at a seen level | 3.5 ms | 3.5 ms |
| `fitted()`, three rows at a new level | 6.0 ms | 189.0 ms |
| `fitted()`, hurdle, latent and `me()` on two predictors | 640.0 ms | 2793.3 ms |

The new-level call now integrates (one direction, 41 nodes plus the
slope's 82); the hurdle call has two directions, 441 nodes, and the band's
direction. The count of nodes does not grow with the rows.

## Not run: the OpenBLAS pass

`bash dev/ciharden-openblas.sh 0.3.26 lapack` was denied by the session's
permission classifier (it downloads an OpenBLAS release with `gh` and copies
the R installation). Per the round's rules it was not retried another way.
The numerical changes (the Gauss-Hermite product rule, the per-row
Cholesky, the body's Jacobian) have no OpenBLAS run; it needs the user.
- Gated tier, `FRMTMB_BRMS_FIT_TESTS=true`, the 46 files that skipped on
  the brms or Stan gate plus `test-me-newdata-brms.R` and
  `test-latres-newrow-brms.R` (`dev/readers-out/gated1.log`): 48 files,
  pass 4526, fail 0, err 0, skip 0, warn 0. Among them
  `test-latres-brms.R` 281, `test-me-newdata-brms.R` 53,
  `test-nl-brms.R` 91, `test-latres-newrow-brms.R` 2.
- One session per suite, three orders (`dev/readers-run-sessions.sh`, copy
  of the release driver on `rellib-r18`; `dev/readers-out/sess1.log`):
  frmtmb forward pass 22562, reverse 22712, preload 22805, fail 0, err 0,
  warn 0 in each; frmtmb.sample 3967 pass, fail 0 in each order. Every
  LEAK line is tagged `[with load of ...]`, an extension's own
  registration; none is untagged.
- `R CMD check --as-cran`, built and checked in `dev/readers-check/`
  (`R_LIBS` lane, rellib-r18, user library; pandoc and TinyTeX on PATH;
  `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`): frmtmb `Status: 1 NOTE`,
  frmtmb.sample `Status: 1 NOTE`. Both NOTEs are the HTML manual's
  "Skipping checking math rendering: package 'V8' unavailable".

## Bump and floor

- frmtmb: minor. Four behavior changes are marked BREAKING in NEWS.
- frmtmb.sample: patch or minor with the round; its `frmtmb (>= ...)` floor
  must rise to the frmtmb release that carries this lane, because
  `draws_refuse_new_levels()` sets `fit$me_new_levels` and only this
  core's `me_check_new_levels()` reads it. On 0.79.1 core the mark is
  ignored and the call is refused as before, so an older core fails loud,
  not silent. No other extension changes or needs a new floor.

## The survey and the ledger

The equivalence survey reads `frm()`, `logLik()` and gradients, not
`fitted()` or `conditional_effects()`, so no sampled pair reaches this
lane's change and none was remeasured. No row of
`dev/brmsport-verdicts.tsv` that I can identify reads `se__`, a new
level's fitted mean or a body's `Est.Error`; the ledger is regenerated at
consolidation.

## Defects found and not fixed

- The kriging variance of an exact `gp()` at a new position is not
  integrated into `fitted()`'s mean (see "What this lane did not do").
- A two-direction `fitted()` costs about 2.8 s per call
  (`dev/readers-timing.R`); a coarser grid would cut it if it matters.
- The other seven extensions in one session, three orders each
  (`dev/readers-out/sess2.log`): 21 of 21 sessions, fail 0, err 0, warn 0
  in all; pass counts equal the one-file runs (coupling 551, eam 1848,
  latent 683, learn 474, ode 549, survival 552, spatial 478). The 12 LEAK
  lines are all tagged `[with load of ...]`.

## Fix round (review `dev/reviews/2026-11-06-readers.md`, verdict FIX)

Everything below is on the lane build after the fix (`wt-readers-lib`,
reinstalled after the last source edit) unless a line says otherwise. The
reviewer's scripts were rerun from `dev/review-readers-*` with `REV_ARM`
unset; logs `dev/readers-out/fix-*.txt`.

### M1 and M2: divergence decided from the model

- `me_mean_diverges()`: the exact mean is infinite when, and only when, the
  family is a lognormal (`lognormal`, `shifted_lognormal`,
  `hurdle_lognormal`), the mean is read (no `dpar`), and a draw moves
  `sigma`: `E[exp(c exp(2 u))]` diverges at any sd of `u`. No numeric tail
  test is left; `me_gh_moments()` no longer returns `inf`.
- Where it diverges, `fitted()` returns the Gauss-Hermite mean and sd on
  `me_gh_div = 11` nodes a direction (reach 5.19 sds), with a warning of
  class `frmtmb_mean_not_finite` that says the exact mean is infinite and
  gives how far 9 nodes move the value. The `Q` columns are kept.
- The choice of rule. The coordinator asked for "the Gauss-Hermite value".
  The 41-node rule would add about 1e16 on the reviewer's case 2: its outer
  node is at 11.6 sds, where `sigma` is about 14.6, `exp(14.6^2 / 2)` is
  2e46 and the weight is 2.3e-30. So the value is a rule of stated reach.
  The alternatives the user may prefer: cut the integral at K sds and state
  K (this is that, with K = 5.19 set by the rule), or refuse by class.
- `dev/review-readers-inf.R` (`fix-inf.txt`):

  | case | review build | fix | reference |
  |---|---|---|---|
  | 2: `sigma ~ (1 \| g)`, sd 0.2745, new level | `Inf` | 2.089686, warns | 4e6 draws 2.089732; 20 x 4000 draws 2.0805 to 2.1031 |
  | 1: `sigma ~ me()`, link sd 0.005 to 0.149 | 2.6907 to 2.7825 | the same, and each warns | |
  | 1: link sd 0.249 | `Inf` | 3.00858, warns | |
  | 1: link sd 0.399 | `Inf` | 9247987, warns | 20 x 4000 draws 4.47 to 61.2: no stable mean |
  | 3: poisson `me()`, link sd 3 to 8 | `Inf` from 6.6; `Est.Error` 0.0047 of exact at 6.4 | ratio to `exp(eta + s^2 / 2)` 1, `Est.Error` ratio 1, at every sd | closed form |
  | 6: shifted lognormal, `sz` 0.8 | `Inf` | 1.953374, warns | |
  | 5: lognormal `mu` only, link sd 5 | 0.99994 of closed form | 1 | closed form |

  At link sd 0.4 the value is not useful, and neither is any mean of
  draws; the warning's 9-node figure says how far it moves.

### M3: a body's band under `me()`

`me_row_loads()` adds, for each body among the read dpars, the nonlinear
parameters it reads (`me_body_parts()`), so `conditional_effects()`, which
asks for `mu`, integrates the row as `fitted()` does. `ce_me_band()` now
sets `se__` to `fitted()`'s `Est.Error`, `sqrt((cf slope)^2 + var)`.
`dev/review-readers-body.R` part B (`fix-body.txt`, data seed 8202, MC
seed 41, 1e6 draws):

| link sd | ce `se__` review / fix | MC sd | ce band, fix | MC 2.5 / 97.5 % |
|---|---|---|---|---|
| 0.25 | 0.4548 / 0.4766 | 0.4757 | 1.1117 to 2.9666 | 1.1121 to 2.9656 |
| 0.50 | 0.9084 / 1.0971 | 1.0968 | 0.6814 to 4.8401 | 0.6810 to 4.8427 |
| 1.00 | 1.8163 / 3.9251 | 3.9286 | 0.2558 to 12.8929 | 0.2554 to 12.9135 |
| 1.50 | 2.7242 / 16.2972 | 16.2536 | 0.0960 to 34.3506 | 0.0958 to 34.4352 |

Pinned by "a body's band under me() is fitted()'s, and a draw's own".

### M4: the Cholesky cross term

New test "the per-row Cholesky factor is chol() row by row": six random
3 x 3 covariances (seed 3301) against `chol()`, to 1e-12 of the largest
entry. On the reviewer's mutant build (`rev-readers-lib`),
`REV_MUT=chol_cross` fails it 12 times (`mut-chol-cross.txt`); with
`REV_MUT` unset it passes (`mut-control.txt`: 0 fails in that block).

### Fix-brief item 4: closed forms and cost

- `me_log_linear()`: when the draws move one predictor and the value read
  is a constant times `exp(eta)` of it (a log link whose inverse is the
  mean, with no truncation, lag SAR or response-side latent; a lognormal's
  `mu`; a dpar read on a log link), `fitted_me_moments()` uses
  `exp(s2 / 2)` times the value at 0 for the mean, computed as
  `exp(log(v0) + s2 / 2)`, and `m^2 expm1(s2)` for the variance, with no
  nodes. Draws on two predictors have no closed form. A draw that moves
  nothing read is dropped first (`me_keep_dirs()`), so a Gamma or
  negbinomial with a draw on `shape` is closed form.
- `mean_only = TRUE` skips the second moment.
- Times (`dev/readers-timing.R`, minimum of 5; `timing-fix.txt` and
  `timing-fix-before.txt`, run back to back; the control `frm_linpred()`
  read 0.89 ms on 0.79.1 and 0.61 ms on the fix):

  | call | 0.79.1 | first round | fix |
  |---|---|---|---|
  | `fitted()`, 3 rows at a new level | 3.0 ms | 189.0 ms | 9.0 ms |
  | `fitted()`, hurdle, 2 directions | 283.3 ms | 2793.3 ms | 993.3 ms |
  | `fitted_me_moments(mean_only = TRUE)`, new level | | | 1.0 ms |
  | the same, hurdle, 2 directions | | | 583.3 ms |

  The first-round column ran under a heavier load (its control read
  1.22 ms).
- After the merge with mfxse (44 `get_predict()` calls per
  `avg_predictions()`, the review's count): about 44 x 1.0 ms, 0.05 s, on a
  one-direction log-link new level, and 44 x 0.58 s, about 26 s, on the
  two-direction hurdle with a latent. Estimated from per-call times, not
  run on a merge.
- The body Jacobian takes one shift input per part, not one per row
  (`dev/readers-bodycost.R`, `bodycost-fix.txt`, memory peak reset per
  size): 243 / 178 / 192 MB at 1000 / 4000 / 20000 rows, which is the
  session's own size. The review measured 488 MB at 4000 rows. The time at
  20000 rows is 16.6 s with the draw and 16.3 s with `sx = 1e-12`, so it is
  the body's se.fit itself, not the draw.

### Minors

- m1: NEWS and `?fitted` name brms's `sample_new_levels = "gaussian"`;
  `"uncertainty"` stays refused (the user's Q2).
- m2: `?fitted` states that `Est.Error` at a new level holds the sd as
  known (0.90 of brms's on the review's poisson). No delta term in theta
  was added.
- m3: `me_keep_dirs()` drops a direction whose unit step moves no cell and,
  past four, keeps the four that move the value most. The review's mixture
  (`dev/review-readers-dirs.R` B, `fix-dirs.txt`): 1.072292, Est.Error
  1.550123, against MC 1.073573 +- 0.001094 and sd 1.547351, with no
  `frmtmb_draw_dims` warning, since the two `sigma` directions move nothing.
  Part A (three directions) is unchanged: 1.55807 and 1.80258 against MC
  1.55777 (0.00080) and 1.80451 (0.00102).
- m4: past two directions, `me_rule_check()` compares the rule with two
  fewer nodes and warns (`frmtmb_draw_rule`) when the mean moves by more
  than 1e-3 or the sd by more than 2 percent. `dev/readers-rule.R` (the
  review's part 2c, seed 6103, `rule.txt`): sd error -2.8e-7, -5.9e-3,
  -8.1e-2 and -0.34 at mu link sd 1.01, 2.02, 2.53 and 3.03; the warning
  is silent at the first and fires at the other three. The sd is not made
  more exact there; one direction on a log link is now exact.
- m5: `?fitted` states the two-direction `Q` coverage, 0.926 to 0.957 on
  the review's four hurdle rows.
- m6: `?fitted` says `fitted()` integrates a `me()` value that a nonlinear
  parameter reads, and `frm_linpred(se.fit)` and `frm_lp_basis()` read it
  at first order.
- m7: `lp_basis_nl()` builds a row of seen levels under the flag; a new
  level that a body reads refuses `se.fit` with class
  `frmtmb_nl_new_level` and a reason, and `fitted()` takes `Est.Error` from
  the integral and its coefficient part (`fitted_point_se()` reads that
  first now). Review part A (`fix-body.txt`): `Est.Error` 0.01541 /
  0.01241 with no group term, 0.02653 / 0.01129 at a seen level (equal to
  the call without the flag), 0.2219 / 0.06746 at a new level; the review
  build read `NA` in all three.
- m8: above.
- m9: frmtmb.sample NEWS no longer names the internal mark.

### Merge note for consolidation

`frm_linpred()` (its body branch) and `predict_mean_se()` (its body branch)
are also edited by lane latperf, which makes a latent-residual fit's joint
covariance sparse. After the merge, run `frm_linpred(se.fit = TRUE)`,
`fitted()` and `conditional_effects()` on a nonlinear body beside
`ar(time = t, gr = s, cov = TRUE)`, on new rows and in sample, and compare
`Est.Error` with `frm_lp_basis()`'s `A V A'`. This lane has no such test:
the base refuses `a * x` with `ar(cov = TRUE)` (`test-latres.R`), and a
nonlinear `shape` or `zi` beside a linear `mu` with a latent is the case
that reaches both. With mfxse, time `avg_predictions()` on a two-direction
row (estimate above: about 26 s).

### Seen to fail, fix round

`test-predict-draw-mean.R`, 120 expectations on the fix:

- On `rellib-r18` (`newtest-fix-before.txt`): pass 29, fail 34, err 6;
  every block fails.
- On the first-round lane build (`rev-readers-lib`, `REV_MUT` unset,
  `mut-control.txt`, run before the rule-check block was added): pass 89,
  fail 28, err 1. The blocks that fail there are the fix round's: the
  lognormal quadrature mean 7, the log-linear sweep 6, the body band 10,
  the body under the flag 4, the mixture 2.

### Reruns, fix round

- The reviewer's 18 files, one per process (`dev/readers-out/fixrev.log`):
  0 fail, 0 err, 0 warn in each. `test-predict-draw-mean.R` 120,
  `test-nl.R` 92, `test-nlf.R` 78, `test-predict-lp-basis.R` 30,
  `test-method-residue.R` 67, `test-me-newdata.R` 178, `test-ce-aterms.R`
  19 (1 skip), `test-ce-bands.R` 170, `test-ce-facets.R` 92,
  `test-ce-levels.R` 101, `test-ce-options.R` 32, `test-ce-parity.R` 30,
  `test-ce-plot.R` 73, `test-ce-predict-median.R` 11 (1 skip),
  `test-conditional-smooths.R` 39, `test-conditions-census.R` 11,
  `test-conditions.R` 150, frmtmb.sample `test-me-gr-new-level-draws.R` 5.
- Core, ungated, all 309 files after the fix (`dev/readers-out/core2.log`,
  sums over its RESULT lines): pass 22636, fail 0, err 0, skip 239,
  warn 0; every `lib:` line names `wt-readers-lib`.
- `R CMD check --as-cran` after the fix: frmtmb `Status: 1 NOTE`,
  frmtmb.sample `Status: 1 NOTE`, each the V8 math-rendering NOTE of the
  HTML manual.
- Not rerun in the fix round: the seven other extension suites, the gated
  tier and the one-session runs. The fix touches only core readers that
  those suites passed through in the first round; consolidation runs them.
- Still not run: the OpenBLAS pass (the user's).

## Short-check fix (review `dev/reviews/2026-11-06-readers.md`, "Short check of the fix round", verdict FIX)

### S1: two layers

- The structural rule. `me_lognormal_sigmas()` names the `sigma` dpars of a
  family's lognormal parts: `sigma` of a lognormal, shifted or hurdle
  lognormal, and `sigma<k>` of a mixture whose component `k` is one (read
  from the family's name, `mixture(lognormal, lognormal)`).
  `me_mean_diverges()` follows those dpars into the bodies that feed them
  (`me_body_parts()`), so a draw on a nonlinear parameter that a `sigma`
  body reads counts.
- The net, for a mean the rule does not mark, at one or two directions:
  the 41- or 21-node value against two fewer nodes (`me_rule_change()`).
  Past a change of 1e-3 (`me_net_tol`), or where the value is not finite
  while the value at the draws' mean is, the rule grows by 10 and 20
  nodes. A mean that settles there is kept, with no warning; one that grows
  is read on 11 nodes with a warning of class `frmtmb_mean_unsettled` (it
  also carries `frmtmb_mean_not_finite`).
- `dev/review-readers-sc2.R` (`fix2-sc2.txt`):

  | model | fix-round build | fix | brms-style 20 x 4000 |
  |---|---|---|---|
  | mixture, `sigma1 ~ me()`, link sd 0.166 | 2.41468, no warning | 2.41468, warns | 2.414 to 2.415 |
  | the same, link sd 0.333 | 2.42274 / Est.Error 92.51, no warning | 2.42274 / 0.1554, warns | 2.422 to 2.424 |
  | the same, link sd 0.555 | 2.495e207 / NaN, no warning | 222685 / 2.478e8, warns | 2.455 to 2.666 |
  | lognormal, `nlf(sigma ~ exp(s1))`, `s1 ~ me()`, `sz` 0.3 | 4.67e82, no warning | 2.81272 / 0.1712, warns | |
  | the same, `sz` 0.8 | NaN, no warning | NaN, warns | |

  At `sz` 0.8 the body's `sigma = exp(exp(s1))` overflows inside 5.19 sds,
  so the row reads `NaN` under the warning. At the mixture's link sd 0.555
  the 11-node value is not useful; the warning gives the change.

### The warning once per call

`me_said()` marks the warning in the frame of the outermost function of
the package family, as `me_warn_missing()` does, so a call says it once
with the largest 9-node change over its rows. Pinned: one
`conditional_effects()` call over two condition sets says it once.

### The dropped-draw shortcut

`me_keep_dirs()` now steps each direction by 0.5, 1 and 1.7 sds each way.
The reviewer's `sin(a)^2 + c` at `a = 0` and a draw sd of pi
(`dev/readers-sin.R`, data seed 7208, `sin.txt`): `fitted()` 1.009284
against the exact 1.009284 (relative -5.7e-13). Its rule did not settle at
21 nodes (a change of 1.2 percent from 19): it settled at 41 under the net,
with no warning. The fix-round build read 0.5093.

### Cost of the net

`dev/readers-timing.R` (`timing-fix2.txt`, control 0.49 ms): `fitted()` at
a new level 9.0 ms (closed form, no net), the two-direction hurdle 1330 ms
(993 ms before the net), `mean_only` on the hurdle 987 ms (583 ms). The
mfxse estimate for that row rises to about 44 x 0.99 s, 43 s per
`avg_predictions()`.

### Seen to fail, short-check fix

`test-predict-draw-mean.R` on the fix-round build (`rev-readers-lib`,
`REV_MUT` unset, `sc-before.txt`): pass 123, fail 7. All three new blocks
fail: the lognormal component and body 4, the warning once per call 2, the
period 1. On the lane (`newtest-fix5.txt`): pass 130, fail 0, warn 0.

### Reruns

The reviewer's 18 files, one per process (`dev/readers-out/fixrev2.log`):
0 fail, 0 err, 0 warn. `test-predict-draw-mean.R` 130, `test-nl.R` 92,
`test-nlf.R` 78, `test-predict-lp-basis.R` 30, `test-method-residue.R` 67,
`test-me-newdata.R` 178, `test-ce-aterms.R` 19 (1 skip), `test-ce-bands.R`
170, `test-ce-facets.R` 92, `test-ce-levels.R` 101, `test-ce-options.R` 32,
`test-ce-parity.R` 30, `test-ce-plot.R` 73, `test-ce-predict-median.R` 11
(1 skip), `test-conditional-smooths.R` 39, `test-conditions-census.R` 11,
`test-conditions.R` 150, frmtmb.sample `test-me-gr-new-level-draws.R` 5.
Code changed only in `R/me.R` and `R/predict.R`; the help of `fitted()` and
NEWS changed with it, so core was checked again (below).

`R CMD check --as-cran` of core after this fix (its tests run the whole
suite in one session): `Status: 1 NOTE`, the V8 math-rendering NOTE of the
HTML manual. frmtmb.sample did not change and was not checked again.

## Confirmation minor: irrational steps

`me_keep_dirs()` steps each direction by sqrt(2) / 2, the golden ratio / 1.3
and sqrt(3) sds each way (`me_keep_steps`). The net now decides row by
row (`me_row_change()`, `me_rows_take()`), so a row that does not settle
no longer moves its neighbours to 11 nodes. `dev/readers-sin.R` (seed 7208,
`sin3.txt`), one call over both rows: at a draw sd of pi 1.009284 against
the exact 1.009284; at 10 pi the draw is kept, the rule cannot integrate
cos(20 pi z) and does not settle, and the row reads 0.8826 against 1.0093
under the `frmtmb_mean_unsettled` warning, not 0.509 in silence.
`test-predict-draw-mean.R`: lane pass 134, fail 0 (`newtest-fix7.txt`); the
fix-round build fails the period block 5 times (`period-before.txt`).
