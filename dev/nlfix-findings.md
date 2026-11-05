# Lane nlfix: nonlinear bodies against brms 2.23.0

Round opened 2026-10-29, base d5d0056d (frmtmb 0.78.0). Worktree
`frmtmb-wt-nlfix`, library `wt-nlfix-lib2`. The tables below are generated
by `dev/nlfix-report.R` from the local logs. Scripts are `dev/nlfix-*`.

## 0a. Fix round 2

The short check of round 1 (`dev/reviews/2026-11-02-nlfix.md`, "Short check of fix
round 1") ruled FIX: 3 majors and 5 minors. This section supersedes what section 0 says of
the stop, of `nl_restart_credible()`, of the return code in the flag and of the
first-body restarts. Each item, what was done, and its evidence:

- **R1, a second body stayed on its saddle.** `nl_saddle_find()` now returns every
  constant body (`$saddle`, `$local`), and one restart moves the coefficients of all of
  them. The kept fit is looked at again, and a body not yet moved gets a further pass.
  Sign flips are per body, with an earlier body held at the kept fit: four flips shared
  over two bodies never reached a negative slope on the second, 3.6 nats short on `mvbf()`.
  `dev/nlfix-r1.log`, seed 11001:

  | model | round 1 | now | every coefficient moved, or the separate fits |
  |---|---|---|---|
  | gaussian, `mu` and `sigma` bodies | -422.3639 | -405.5429 | -405.5429 |
  | mixture, `mu1` and `mu2` bodies | -470.6004 | -404.9296 | -404.9296, linear too |
  | `mvbf()`, two ordinal bodies | -750.5521, first form of this fix | -746.9868 | -746.9868 |

  The `mvbf()` test starts from a first fit with both bodies on the point
  (`iter.max = 0`), so nothing rests on its joint start. The two gaussian and mixture
  designs are tested, with the range of the second body's predictor.
- **R2, the stop.** `frm_stop()` became a warning of class `frmtmb_nl_saddle` in every
  case, naming each body, saying the data may simply have no effect there and that a
  `start` can test it. `nl_restart_credible()` is removed. `dev/nlfix-r2.log`: slope 0 and
  slope 2e-5 on `a + b1 * b2 * x` return the linear value with one warning; 400 null data
  sets, seed 11002: no error, 1 warning, where one stopped; sigma 1e-3 to 1e-8, seed 11004,
  the default start reaches the linear value at code 0 (2475.3633 at 1e-6). At 1e-9 both
  the default start (3663.80) and a started fit (3677.15) stop short of 3856.91 at code 1,
  with the convergence warning: an optimizer limit, not this rule.
- **R3, `beta_binomial()`'s density.** `beta_binomial_lpmf()` replaces
  `RTMBdist::dbetabinom()` in both lpdfs, built on `lgamma_count_shift()`. `osa` and
  `simref` responses still go to `dbetabinom()`, which refuses the one and simulates the
  other. Before and after, the review's script (`dev/nlfix-rc-c4-before.log`, on
  rellib-r16; `dev/nlfix-rc-c4-after.log`):

  | seed | before, linear | before, body | after, linear and body | binomial limit |
  |---|---|---|---|---|
  | 2 | -541.0861, code 1 | -541.2367, code 1 | -541.2368, code 0 | -541.2368 |
  | 4 | -543.7369, code 1 | -543.7369, code 1 | -543.7369, code 0 | -543.7369 |
  | 6 | -541.4800, code 1 | +94400, code 1 | -541.4800, code 0 | -541.4800 |
  | 9 | +94208, code 0 | -548.9821, code 1 | -548.9821, code 0 | -548.9821 |
  | 12 | -543.5460, code 1 | -543.5460, code 1 | -543.5460, code 0 | -543.5460 |

  The grid at seed 9 stays within 2e-11 of the finite product from log phi -3 to 100
  (`dev/nlfix-rc-c5bbgrid-after.log`). The zero-inflated lpdf, which the review did not run:
  `dev/nlfix-lgscan-after.log` holds -543.82211 at every log phi from 20 to 100, where
  rellib-r16 reads -80064 at 39, +196608 at 45 and 1.5e+30 at 100; all 5 seeds reach code 0.
  `test-beta-binomial-density.R` holds both against the finite product by ratio. Overdispersed
  fits keep their value: the cost table's logLik agree to 6 decimals, and the test holds
  trials = 200 within 1e-11 relative of `dbetabinom()` and of the finite product.

  **Cost, profiled.** The time is proportional to the tape (retape of the gradient, MakeADFun,
  the optimizer), not to the R code that writes it: a sparse sum over the rows, which saves
  the per-step copies, changes nothing. So the lever is the node count. The sum of `log(a + j)`
  is about 3 nodes a step against about 80 for `lgamma_shift_diff()`, counted from the code, and it is exact, so rows with
  a shift of at most 24 use it (`dev/nlfix-bbcost.log`: 6, 12, 24, 48 and 96 measured). Whole
  fits of the shipped families at n = 3e4, min of 4, `dev/nlfix-bbcost2-*.log`:

  | trials | before | after | ratio |
  |---|---|---|---|
  | 5 | 2.71 s | 1.92 s | 0.71 |
  | 10 | 3.34 s | 2.89 s | 0.87 |
  | 24 | 3.34 s | 5.27 s | 1.58 |
  | 30 | 2.97 s | 5.77 s | 1.94 |
  | 200 | 3.42 s | 10.58 s | 3.09 |

  The reviewer's form costs 3.0 to 3.7 at 10, 30 and 200 trials (`dev/nlfix-bbcost.log`, as a custom family; 2.89 in the review). Above 24 trials the price
  stays, because `size` is long on every row and nothing there can be cached: the data terms
  are already off the tape. Zero-inflated: 0.61, 0.81, 1.48, 1.73, 2.67.

  **The other users of these differences.**
  - The CDF, `beta_binomial_cdf_log()`: had the same raw form and is fixed, used by `cens()`
    and `trunc()` on both families. Calling the new mass at each support point would have
    paid for the whole sum at every point, so it reads a table of the shifts instead, one
    `log()` per step. Censored fits at n = 3000: 0.55 to 0.50 s at 10 trials, 6.52 to 6.97 s
    at 200 (`dev/nlfix-bbcost2-cens-*.log`), logLik equal to 6 decimals. Tested against the
    finite product at log phi 2 to 45, per-row trials, and the support's edges.
  - `simulate()`: `rbetabinom()` reads `rbeta()`, no `lgamma()`, and is sound at log phi 39
    and 60 (`dev/nlfix-lgscan-*.log`). Filed below: it stops at a mean of exactly 0 or 1.
  - `dirichlet_multinomial()`: already `lgamma_shift_diff()`. Flat to 8 decimals from log phi
    30 to 100. Unchanged.
  - `negbinomial()`, `nbinom1()`, `zero_inflated_negbinomial()`, `hurdle_negbinomial()`: go
    through `dnbinom_robust()`, flat to 8 digits from log precision 15 to 100 on Poisson data,
    5 seeds. Unchanged.
  - `Beta()` and its inflated forms (`RTMB::dbeta()`): raw, and rounding at log phi 30 and up
    (-lbeta off by 800 at 39). A plain fit does not go there, since that density falls with
    phi unless every y equals its mean. A latent residual does, and that is latres's cap
    (`dirichlet()`, `R/families-simplex.R`, is the same case). Unchanged, filed below.
  - `R/families.R:2340` is a unit deviance, off the tape.

  **The check.** `discrete_loglik_warning()`, class `frmtmb_discrete_loglik`, in
  `fit_end_checks()`: every response discrete and logLik above 0, past rounding of 1e-8 per
  row. Not under a prior, where `logLik()` holds the log prior, nor under REML. Tested with a
  custom family that drops a constant. On rellib-r16's +94208 it would have fired; it is
  one-sided, and that is what the density fix is for.

  **Marked BREAKING** in NEWS, since these fits change. **For latres:** rerun the R1 table
  (`dev/latres-r1.R`) after consolidation. Expect the `beta_binomial()` "no cap" row to lose
  its 6.5e+17 and read a largest logLik of -298, the capped arms' value, with fewer nonzero
  codes than 9, a largest log dispersion near 20 to 22 (where the density has gone flat and
  the gradient has no more to say), and the capped arms unmoved within about 1e-5 (at log phi
  15 the density is within 7e-6 of its limit on seed 9). So the box at 15 may no longer be
  needed for this family, and it does no harm. Expect the `Beta()` (1.9e+30) and
  `zero_one_inflated_beta()` (4.1e+07) rows not to change: this round does not touch those
  densities, so `latres_disp_upper()` stays for them.
- **Ruling (a).** A constant body the curvature exempts now gets a message of class
  `frmtmb_nl_local`: "... which is the thresholds-only model ... The point is a local maximum
  along its coefficients. A `start` elsewhere can test for a higher one". Seed 11005,
  `b^2 * x * (1 - b^2)`: -412.6343, against brms's -383.7955 (`dev/nlfix-rc-c1b-after.log`).
  A `start` of 1.5 returns to the same point there, so the sentence offers a test and not a
  cure.
- **Ruling (b).** The return-code arm is gone: only the flat direction in the standard errors
  flags a fit. The acat ridge, on three BLASes (`dev/nlfix-acat.log`,
  `dev/nlfix-acat-ob0.3.26.log`, `dev/nlfix-acat-ob0.3.32.log`): -372.2279 from the default
  start and under `se = FALSE`. Under `check_se = "ignore"`, all three stay at -378.0992. So
  rounding no longer decides, and that one setting turns the check off along with what it
  reports: a limit, stated here. `iter.max = 0` bypasses the wrapper: -413.5626, estimates
  0, no message. `iter.max = 3` no longer restarts a capped fit that has left the point.
  Review's script, `dev/nlfix-rc-c3b-after.log`, ratio to one fit: `a + b` RE 2.4 (5.0 in the
  review), `b0 + b1 * x` cumulative 1.8 (4.1), `iter.max = 3` RE 1.0 (4.0).
- **Two agreeing restarts stop the round, if their estimates agree too.** Found by measuring:
  equal values are not enough. `c1 * exp(c2 * z)` with a negative slope reads -393.0311 from the
  base start and the first flip, at c1 = -12077 and -14523, and -388.3165 only from the second
  flip (`dev/nlfix-cost3b.log`). `b1 * (1 - exp(-b2 * x))` under `cratio()` and `sratio()` reads
  -377.0010 twice, at b1 = 21293 and 22565, then -373.2498; under `bernoulli()`, -180.7943 twice,
  then -180.3394 (`dev/nlfix-cost4.log`). A first form of this rule that read values and the flat
  warning returned those three fits short, at code 0 on the ordinals, and the saddle table showed
  it. The rule is now values within 1e-7 and estimates within 1e-3. A pair at the saddle's own
  value does not stop it, and on a flat direction a first restart that agrees with the first fit
  does. Saddle table: at or above brms's best on 50 of 51 rows, the same as round 1, and the one
  short is bernoulli's true second maximum. **The price, stated plainly:** a ridge is the case the
  review timed, so on such data the third restart still runs. `b1 * exp(b2 * x)` on a linear truth,
  3e4 rows, seed 6401 (`dev/nlfix-cost3b.log`): restarts of 8.0, 8.3 and 2.6 s, 20.1 s in all, and
  the third is the one that reads highest, -38506.5831 against -38506.8830. So the review's "restart
  3 alone wastes 12 s" is saved only where the first two land on one point. An identified body
  stops at two. `a + b1 * b2 * x` stops where two runs land on one point of its flat. A limit stands
  for a restart from a flat direction: it stops at the first agreement with the first fit, as the
  ruling asks, so a flip it did not try is not tried.
- **The restart's message.** "reached a higher point, log-likelihood X, optimizer return code
  C", with the digits that tell the two values apart. No "better maximum". A model is named
  only from `nl_model_words()`: "thresholds-only" needs a `mu` body, every body of the response
  constant, and nothing else on it that varies; the message of a kept restart names none.
- **The help.** `?par_template` and `start` in `?frm`. Measured first (`dev/nlfix-bern2.log`):
  `frm_allfit()` does not show bernoulli's second maximum, either way: all four optimizers read
  -193.4967 from the estimates, and from the original start. A `start` of `b1 = 1` reads
  -189.999. So the sentence says a `start` can and `frm_allfit()` cannot, not the "or" the
  message suggested.
- **The `post$fit_check` warning on `nlf(disc ~ ...)`.** Filed first, then fixed on the
  coordinator's addition, from predread's trial merge. `ord_fit_check()` read `lp[["X"]]`, which a
  body does not have, so `colnames(NULL)` raised "'data' must be of a vector type", the fit-end
  guard turned that into the hook's warning, and the check never ran; my test allowed the
  warning. It now returns on a body, the allowance is gone, and a test holds both spellings. Why
  it returns and does not read the body: finding a free constant in an arbitrary body is what the
  standard errors already do, and the same test shows they report `d0 + d1 * z` as flat.
- **The abandoned session run (`dev/release/nlfixsess-sessions/`).** Not clean, and not a
  fault of the code: my own stop of the chain caused it. The first final chain ran on
  `wt-nlfix-lib`. To pick up the `ord_fit_check()` addition I stopped its background task at
  00:51:10, and that killed only the outer shell. The chain's inner bash and its sessions went
  on, and from that moment they could not start a child R process. The evidence is the order in
  which the files fell. Every failing file spawns a child (`system2`, or
  `helper-absent-package.R`). In forward and preload order only the files after the stop fail:
  `test-refit-data.R` and `test-registry-helper.R` (fail=4, err=1). `test-generic-collision.R`
  passes there. In reverse order only the files before it fail: `test-generic-collision.R`
  (75), `test-extension-home-mixed.R` (13), `test-guardrail1.R` (9), `test-cs-bar-brms.R` (7),
  `test-conditions.R` and `test-compat-extensions.R` (1 each), with an empty `out`.
  `test-refit-data.R` passes there. `frmtmb.sample` reverse: 37 in the same file. The `frmtmb.sample`
  preload session, which starts a PSOCK cluster, wrote its last line at 00:53:37 and hangs (PID
  7544): a worker that did not start is the likely cause, inferred and not confirmed. The
  mechanism itself (why a killed parent stops children from starting) is inferred too. The same
  sessions were clean in round 1, and they are clean on the final build below, which ran without
  a stop. PID 7544 was still alive at 05:31, holding `wt-nlfix-lib`'s DLLs. Stopping it was
  denied by the auto mode classifier, so it stays for the user. `wt-nlfix-lib` is therefore
  stale: the final build is `wt-nlfix-lib2`.
- **`frmtmb.learn` skip=13 and `frmtmb.spatial` skip=6** are gates, not load failures: 11 + 2 and
  6 reasons, all `FRMTMB_BRMS_FIT_TESTS` or `FRMTMB_SCALE_TESTS`, the same counts as round 1, and
  their passes are unchanged.
- **Final pass, `wt-nlfix-lib2`, every `lib:` line read.** Core and all eight extensions installed
  there, so no `lib:` line mixes builds. R sources dated before the core check's build (07:16:50).

  | tier | result |
  |---|---|
  | core, one file per process | 290 files, PASS 20980, FAIL 0, ERROR 0, SKIP 217, WARN 0 |
  | gated | 35 files, PASS 3040, FAIL 0, ERROR 0, SKIP 15 |
  | extensions | 183 files, PASS 8861, FAIL 0, ERROR 0, SKIP 47 |
  | frmtmb.sample gated | 7 files, PASS 731, FAIL 0, ERROR 0, SKIP 1 |
  | sessions, tag `nlfixsess3` | 27 of 27, fail 0, err 0; 50 LEAK lines, as in round 1, all "with load of" |
  | OpenBLAS 0.3.26 and 0.3.32 with LAPACK | 7 files each, FAIL 0, ERROR 0: nl-family 321, nl-stan-funs 122, nl-brms 91, beta-binomial-density 93, cdfs 234, cens-trunc 72, numerical-robustness 723 |
  | `R CMD check --as-cran` | frmtmb, sample, eam, survival: 1 NOTE (HTML version of manual); learn: OK |
  | a8 | 530 of 530 identical |
  | equiv | 35 pairs, 2 changed verdict: p049, p158 |
  | sweep | 41 nl honored, 7 differ from both, the same 7 families as before |

- **Filed, not fixed:** section 7.

The review (`dev/reviews/2026-11-02-nlfix.md`) ruled FIX: 2 blockers, 2 majors, 13
minors and the NITs. The `is_primary` fix held bitwise on 530 of 530 items, and the L1
identities stand. Each finding, what was done, and its evidence:

- **Blocker 1, a default start returned the thresholds-only model.** Cause: at
  b1 = b2 = 0 the body is constant, the thresholds start at their exact maximum, so
  every gradient of a body coefficient is 0 and the optimizer stops at once, code 0.
  My specs all carried a `start`, and L1 cannot see it, since the objectives agree at
  the saddle. `nl_fit_off_saddle()` wraps frm()'s fit of any model with a body. A fit
  whose body varies across the rows and whose optimizer converged is returned as it
  was. On a body constant at the stop that can vary at a moved point
  (`nl_saddle_find()`), the curvature along the body's coefficients decides: strictly
  convex negative log-likelihood is a true maximum and is left alone
  (`a + b^2 * x` with a negative slope in the data, tested). Otherwise frm() restarts
  from up to 4 deterministic moved starts (`nl_saddle_starts()`, no RNG, tested with two
  seeds), keeps a higher point and says so in a message, and the first fit's warnings
  are dropped with it. Round 2 (section 0a) replaced the stop with a warning, moved every
  body, removed the return code from the flag and removed the smoothness test on a restart.
  A default start whose standard errors report a flat direction also gets the
  restarts, because `acat()` runs up the ridge to the linear limit (-378.10): that case
  has no constant body. It keeps the first fit unless a restart is higher.

  From the default start, every body and family, with brms's `optimizing` from 6 starts
  (`dev/nlfix-saddle.R`, seeds 9202, 1 and 9205):

<!-- nlfix:saddle begin -->
51 of 52 cases have a row. frm() from the default start is at or above brms's best of 6 on 50; `restarted` is 32. No row: poisson_id_sat.

| family | body | default | started | no_body | brms_best | brms_minus_frm | restarted | conv | err |
|---|---|---|---|---|---|---|---|---|---|
| acat | b1 * exp(b2 * x) | -372.2279 | -372.2279 | -413.5626 | -372.2279 | -7.1e-08 | TRUE | 0 |  |
| acat | b1 / (1 + exp(b2 * x)) | -378.0987 | -378.0987 | -413.5626 | -378.0989 | -1.6e-04 | FALSE | 0 |  |
| acat | a + b1 * b2 * x | -378.0987 | -378.0987 | -413.5626 | -378.0987 | -7.9e-08 | TRUE | 0 |  |
| acat | b1 * (1 - exp(-b2 * x)) | -372.2279 | -378.1001 | -413.5626 | -372.2279 | -2.8e-07 | TRUE | 0 |  |
| bernoulli | b1 * exp(b2 * x) | -193.4967 | -189.9990 | -204.9945 | -189.9990 | 3.5e+00 | FALSE | 0 |  |
| bernoulli | b1 / (1 + exp(b2 * x)) | -187.5073 | -187.5073 | -204.9945 | -187.5073 | -2.7e-09 | FALSE | 0 |  |
| bernoulli | a + b1 * b2 * x | -177.1772 | -177.1772 | -204.9945 | -177.1772 | 9.0e-12 | TRUE | 0 |  |
| bernoulli | b1 * (1 - exp(-b2 * x)) | -180.3394 | -180.7943 | -204.9945 | -180.3394 | -3.7e-11 | TRUE | 0 |  |
| cratio | b1 * exp(b2 * x) | -373.2498 | -373.2498 | -413.5626 | -373.2498 | -2.9e-08 | FALSE | 0 |  |
| cratio | b1 / (1 + exp(b2 * x)) | -377.0006 | -377.0006 | -413.5626 | -377.0008 | -2.0e-04 | FALSE | 0 |  |
| cratio | a + b1 * b2 * x | -377.0006 | -377.0006 | -413.5626 | -377.0006 | -1.1e-08 | TRUE | 0 |  |
| cratio | b1 * (1 - exp(-b2 * x)) | -373.2498 | -377.0020 | -413.5626 | -373.2498 | -2.0e-08 | TRUE | 0 |  |
| cumulative | b1 * exp(b2 * x) | -371.7584 | -371.7584 | -413.5626 | -371.7584 | -9.9e-09 | TRUE | 0 |  |
| cumulative | b1 / (1 + exp(b2 * x)) | -377.2793 | -377.2793 | -413.5626 | -377.2794 | -1.2e-04 | TRUE | 0 |  |
| cumulative | a + b1 * b2 * x | -377.2793 | -377.2793 | -413.5626 | -377.2793 | -5.7e-08 | TRUE | 0 |  |
| cumulative | b1 * (1 - exp(-b2 * x)) | -371.7584 | -377.2799 | -413.5626 | -371.7584 | -3.7e-08 | TRUE | 0 |  |
| gamma_log | b1 * exp(b2 * x) | -350.1928 | -350.1928 | -448.4053 | -350.1928 | -8.5e-11 | FALSE | 0 |  |
| gamma_log | b1 / (1 + exp(b2 * x)) | -354.8181 | -354.8181 | -448.4053 | -354.8181 | -6.6e-10 | FALSE | 0 |  |
| gamma_log | a + b1 * b2 * x | -317.4805 | -317.4805 | -448.4053 | -317.4805 | 1.7e-13 | TRUE | 0 |  |
| gamma_log | b1 * (1 - exp(-b2 * x)) | -357.8319 | -389.9904 | -448.4053 | -357.8319 | -1.6e-07 | TRUE | 0 |  |
| gaussian | b1 * exp(b2 * x) | -194.2316 | -194.2316 | -487.0357 | -194.2316 | -3.7e-13 | FALSE | 0 |  |
| gaussian | b1 / (1 + exp(b2 * x)) | -313.9613 | -313.9613 | -487.0357 | -313.9613 | -2.2e-09 | FALSE | 0 |  |
| gaussian | a + b1 * b2 * x | -286.3017 | -286.3017 | -487.0357 | -286.3017 | -2.0e-10 | TRUE | 0 |  |
| gaussian | b1 * (1 - exp(-b2 * x)) | -615.0010 | -615.0007 | -487.0357 | -614.9997 | 1.3e-03 | TRUE | 1 |  |
| hurdle_cumulative | b1 * exp(b2 * x) | -447.0257 | -447.0257 | -481.6299 | -447.0257 | -1.5e-07 | TRUE | 0 |  |
| hurdle_cumulative | b1 / (1 + exp(b2 * x)) | -451.4982 | -451.4982 | -481.6299 | -451.4983 | -1.7e-04 | TRUE | 0 |  |
| hurdle_cumulative | a + b1 * b2 * x | -451.4982 | -451.4982 | -481.6299 | -451.4982 | -5.8e-08 | TRUE | 0 |  |
| hurdle_cumulative | b1 * (1 - exp(-b2 * x)) | -447.0257 | -451.4994 | -481.6299 | -447.0257 | -4.0e-07 | TRUE | 0 |  |
| hurdle_poisson | b1 * exp(b2 * x) | -517.0674 | -517.0674 | -564.7344 | -517.0674 | -4.3e-11 | FALSE | 0 |  |
| hurdle_poisson | b1 / (1 + exp(b2 * x)) | -517.3206 | -517.3206 | -564.7344 | -517.3206 | -2.1e-10 | FALSE | 0 |  |
| hurdle_poisson | a + b1 * b2 * x | -513.8280 | -513.8280 | -564.7344 | -513.8280 | -4.0e-09 | TRUE | 0 |  |
| hurdle_poisson | b1 * (1 - exp(-b2 * x)) | -583.0238 | -587.8831 | -564.7344 | -587.8831 | -4.9e+00 | TRUE | 0 |  |
| poisson | b1 * exp(b2 * x) | -488.5901 | -488.5901 | -562.9226 | -488.5901 | -1.6e-11 | FALSE | 0 |  |
| poisson_id | b1 * exp(b2 * x) | -473.6457 | -473.6457 | -562.9226 | -473.6457 | -7.6e-12 | FALSE | 0 |  |
| poisson_id | b1 / (1 + exp(b2 * x)) | -483.3459 | -483.3459 | -562.9226 | -483.3459 | 5.7e-12 | FALSE | 0 |  |
| poisson_id | a + b1 * b2 * x | -486.3386 | -486.3386 | -562.9226 | -493.7689 | -7.4e+00 | TRUE | 0 |  |
| poisson | b1 / (1 + exp(b2 * x)) | -487.2919 | -487.2919 | -562.9226 | -487.2919 | -1.3e-11 | FALSE | 0 |  |
| poisson | a + b1 * b2 * x | -473.6457 | -473.6457 | -562.9226 | -473.6457 | -2.4e-10 | TRUE | 0 |  |
| poisson | b1 * (1 - exp(-b2 * x)) | -503.8690 | -507.4595 | -562.9226 | -507.4519 | -3.6e+00 | TRUE | 0 |  |
| probit | b1 * exp(b2 * x) | -372.1003 | -372.1003 | -413.5626 | -372.1003 | -9.7e-10 | TRUE | 0 |  |
| probit | b1 / (1 + exp(b2 * x)) | -377.5830 | -377.5830 | -413.5626 | -377.5831 | -1.5e-04 | TRUE | 0 |  |
| probit | a + b1 * b2 * x | -377.5829 | -377.5829 | -413.5626 | -377.5829 | -3.5e-09 | TRUE | 0 |  |
| probit | b1 * (1 - exp(-b2 * x)) | -372.1003 | -377.5833 | -413.5626 | -372.1003 | -1.1e-07 | TRUE | 0 |  |
| sratio | b1 * exp(b2 * x) | -373.2498 | -373.2498 | -413.5626 | -373.2498 | -2.9e-08 | FALSE | 0 |  |
| sratio | b1 / (1 + exp(b2 * x)) | -377.0006 | -377.0006 | -413.5626 | -377.0008 | -2.0e-04 | FALSE | 0 |  |
| sratio | a + b1 * b2 * x | -377.0006 | -377.0006 | -413.5626 | -377.0006 | -1.1e-08 | TRUE | 0 |  |
| sratio | b1 * (1 - exp(-b2 * x)) | -373.2498 | -377.0020 | -413.5626 | -373.2498 | -2.0e-08 | TRUE | 0 |  |
| zip | b1 * exp(b2 * x) | -439.4212 | -439.4212 | -477.4951 | -439.4212 | -1.3e-09 | FALSE | 0 |  |
| zip | b1 / (1 + exp(b2 * x)) | -443.6930 | -443.6930 | -477.4951 | -443.6930 | -1.4e-12 | FALSE | 0 |  |
| zip | a + b1 * b2 * x | -432.6954 | -432.6954 | -477.4951 | -432.6954 | -5.6e-09 | TRUE | 0 |  |
| zip | b1 * (1 - exp(-b2 * x)) | -449.5135 | -450.8689 | -477.4951 | -450.8690 | -1.4e+00 | TRUE | 0 |  |
<!-- nlfix:saddle end -->

  Reading it. `default` at or above `brms_best` on every row but two. `bernoulli`
  `b1 * exp(b2 * x)` is at -193.4967, a true local maximum and no saddle: gradient
  6.5e-5, reached from 4 of 6 starts of mine at b1 -0.36, b2 -0.87, where a start
  with b1 > 0 reaches brms's -189.999 at b1 0.0017, b2 5.01 (`dev/nlfix-bern.R`). Two
  maxima of a body the data do not fit; a restart rule cannot see it, and frm()
  reports code 0 honestly. `gaussian` `b1 * (1 - exp(-b2 * x))` has no level and is
  equal to brms within 1e-3 of a nat at its ridge, code 1 with its warning. Where
  `default` is above `started`, the fixed `start` of 0.3 stopped on a ridge or a
  worse maximum, as it does in brms. `poisson('identity')` `b1 * (1 - exp(-b2 * x))`
  has no row: frm() stops with "failed from its default starting values", correctly,
  since that body is negative on some rows at every coefficient tried, and the
  script's process then dies in brms's `rstan::optimizing()` arm, not in frm()
  (`dev/nlfix-segv2.R`).

  The restart's cost and that nothing else moved, digit for digit:

<!-- nlfix:a8 begin -->
```
keys only before:  
keys only after:  
identical: 530 of 530 
cases: 53 
```
<!-- nlfix:a8 end -->

<!-- nlfix:cost2 begin -->
```
plain 2e3                n=2000    fit    0.420 s  guard  0.0001 s  domain  0.0005 s  saddle look  0.0003 s  (all)/fit = 0.002
plain 2e5                n=200000  fit    2.750 s  guard  0.0102 s  domain  0.0306 s  saddle look  0.0114 s  (all)/fit = 0.019
RE 2e4                   n=20000   fit    0.670 s  guard  0.0091 s  domain  0.0032 s  saddle look  0.0013 s  (all)/fit = 0.020
no start, linear body    n=200000  fit    1.690 s  guard  0.0103 s  domain  0.0216 s  saddle look  0.0024 s  (all)/fit = 0.020
ar cov 4e3               n=4000    fit    0.400 s  guard  0.0015 s  domain  0.0067 s  saddle look  0.0004 s  (all)/fit = 0.021
cumulative 3e4           n=30000   fit    0.590 s  guard  0.0065 s  domain  0.0000 s  saddle look  0.0003 s  (all)/fit = 0.012
cumulative 3e4, saddle   n=30000   fit   20.410 s  guard  0.0066 s  domain  0.0000 s  saddle look  0.0017 s  (all)/fit = 0.000
the restart's cost, where it runs: a saddle fit against a started fit, data of
saddle fit 0.048 s, started fit 0.014 s, ratio 3.45
```
<!-- nlfix:cost2 end -->

- **Blocker 2, `simulate(newdata =)` under `cov = FALSE` counted the AR part twice.**
  `sim_newdata_dpars()` now reads `nl_body_value()`. Rerun of the review's script:

<!-- nlfix:a4c begin -->
```
slope of nl's row-mean excess on the observed shift: 0.003219543 
slope of lin's row-mean excess on the observed shift: 0.003219868 
```
<!-- nlfix:a4c end -->

  The other callers my change to `frm_linpred()` reaches, read one by one. The linear
  branch and the body branch now take the same route in each, so none differs from a
  linear `mu`: `dpars_natural()` (R/predict.R:1665), `predict_dpar_values()`
  (R/predict-brms.R:1169, 1183), `predict_cross_rows()` (:402), `fitted_point()`
  (R/predict.R:3582, 3592), `ord_thres_linpred()` (:3695), the summary's se call (:3810),
  `conditional_effects()` (R/conditional-effects.R:875, 879, 2650, 2693, 2713, 2740, which
  run on a fit with its `cov = FALSE` terms removed), `curve_*` (R/curve-cov.R:187, 311;
  R/curve.R:410; on a body the curve carries no shift, filed by the review for
  consolidation), `get_predicted` (R/insight.R:191), `predictions` (R/interop.R:936),
  R/spectral.R:629, R/families-extreme.R:491, 800, frmtmb.coupling (R/coupling.R:124, 212),
  frmtmb.eam (R/ddm-shared.R:661), frmtmb.survival (R/rp-check.R:72) and frmtmb.sample
  (R/methods-draws.R:946, 1071, 1361, 1378, 3694, R/conditional-effects-draws.R:190).
  `sim_newdata_dpars()` was the only place with a split, its linear branch on the design
  and its body branch on `frm_linpred()`. `lp_eta_or_body()` reads `nl_body_value()` and
  so carries no shift.

- **Major 1, the guard refused correct models.** It now fires only when the body's value
  moved by more than 1e-4 of the response's sd and the objective did not, on two probes,
  with a separate irrational step per coefficient. It reads the taped gradient first,
  which needs a fixed-effects tape, and asks the objective only for what is left.
  `(2 * a - b) * x`, `2 * a - b` with priors (on `frm_sample()` too), `a - 2 * b + c`, `a *
  x` with x about 1e-17, `b * x` with x about 1e-10, `step(x2 - c)` and `exp(-k * xh)` fit;
  `mu * 0` and `categorical(refcat = NA)` still refuse; a density reading 2 rows passes
  (test-nl-family.R). Rerun of the review's script:

<!-- nlfix:a2 begin -->
```
step(x2 - c), one coef             ll=-233.42633 conv=0 beta=0 WARN[Standard errors are not available for 1 of 2 parameters. c_(Intercept): the likelihood is flat along it, or too 
step(x2 - c), start 15             ll=-41.8046711 conv=0 beta=15 WARN[Standard errors are not available for 1 of 2 parameters. c_(Intercept): the likelihood is flat along it, or to
a * step(x2 - c)                   ll=-173.593842 conv=0 beta=0.47835,0.00000 WARN[Standard errors are not available for 1 of 3 parameters. c_(Intercept): the likelihood is flat al
pmax(x2 - c, 0)                    ERROR The nonlinear formula body could not be evaluated: Comparison is generally unsafe for AD types
pmin(x2, c), flat at start         ERROR The nonlinear formula body could not be evaluated: Comparison is generally unsafe for AD types
ifelse in body                     ll=-41.7083339 conv=0 beta=0.98671 
exp(-k xh) start k=1               ll=137.84281 conv=0 beta=1 WARN[Standard errors are not available for 1 of 2 parameters. k_(Intercept): the likelihood is flat along it, or too n
a tiny multiplier a*1e-300         ll=-300.711233 conv=0 beta=8.7215e-299 WARN[Standard errors are not available for 1 of 2 parameters. a_(Intercept): a step along it raises the lo
tiny x, linear                     ll=-275.780407 conv=0 beta=5.0728e+19 
tiny x, a * xt                     ll=-300.711233 conv=0 beta=8.7215e-19 
2 * a - b                          ll=-300.687086 conv=0 beta=-0.0067647, 0.0033823 WARN[The coefficients a_Intercept, b_Intercept of the nonlinear parameters 'a', 'b' are not iden
2 * a - b + 0 * x, prior           ll=-302.524992 conv=0 beta=-0.0067567, 0.0033783 
(2 * a - b) * x                    ll=-275.780407 conv=0 beta= 0.20291,-0.10146 WARN[The coefficients a_Intercept, b_Intercept of the nonlinear parameters 'a', 'b' are not identifi
a - 2 * b + c                      ll=-300.687086 conv=0 beta=-0.0028182, 0.0056365,-0.0028182 WARN[The coefficients a_Intercept, b_Intercept, c_Intercept of the nonlinear paramete
a*b - b*a (really absent)          ll=-300.711233 conv=0 beta=0,0 WARN[Standard errors are not available for 2 of 3 parameters. a_(Intercept), b_(Intercept): the likelihood is flat
a * x0                             ll=-300.711233 conv=0 beta=0 WARN[Standard errors are not available for 1 of 2 parameters. a_(Intercept): the likelihood is flat along it, or too
a * x0, linear                     ll=-300.711233 conv=0 beta=0 WARN[Standard errors are not available for 1 of 2 parameters. x0: the likelihood is flat along it, or too nearly fla
b + a * x0                         ll=-300.687086 conv=0 beta= 0.00000,-0.01691 WARN[Standard errors are not available for 1 of 3 parameters. a_(Intercept): the likelihood is flat 
a * b * x                          ll=-275.780407 conv=0 beta=0.74405,0.68178 WARN[Standard errors are not available for 2 of 3 parameters. a_(Intercept), b_(Intercept): the likeli
constant prior on a                ERROR Unsupported prior distribution 'constant' (supported: normal, student_t, cauchy, exponential, logistic, gamma, inv_gamma, beta, uniform, lo
constant prior, all                ERROR Unsupported prior distribution 'constant' (supported: normal, student_t, cauchy, exponential, logistic, gamma, inv_gamma, beta, uniform, lo
a ~ 0 + offset                     ll=-275.786235 conv=0 beta= 
bound coefficient                  ll=-275.780407 conv=0 beta=0.50728 
custom, mu at y > 2 only           ll=-275.699329 conv=0 beta= 2.487501,-0.062172 
custom, mu * 0                     ERROR The nonlinear body of 'mu' does not reach the likelihood of family 'rvbranch0': moving the coefficients of 'b0', 'b1' moves the body and le
custom, reads mean(mu) only        ll=-300.711233 conv=0 beta=4.0303e-16 WARN[Standard errors are not available for 1 of 2 parameters. b1_(Intercept): the likelihood is flat along 
nlf(sigma) body                    ERROR Cannot interpret bf() argument 1: expected a dpar formula, a named numeric constant, or a dpar equated to another by name, as sigma1 = "sig
hurdle all-zero y                  ERROR The nonlinear body of 'mu' does not reach the likelihood of family 'hurdle_poisson': moving the coefficients of 'b0', 'b1' moves the body a
hurdle all-zero y linear           ll=-4.64395601e-296 conv=1 beta=0,0 WARN[Optimizer did not report convergence: function evaluation limit reached without convergence (9). The cur
frm_sample body the family skips: ERROR The nonlinear body of 'mu' does not reach the likelihood of family 'rvbranch0': moving the coefficients of 'b0', 'b1' moves the body and lea
log(a) default start               ll=-Inf conv=0 beta=0 WARN[The gradient at the reported optimum is not finite (beta, betad), so convergence could not be checked: the objective's
custom skip + log(a)               ll=-300.711233 conv=0 beta=0 WARN[Standard errors are not available for 1 of 2 parameters. a_(Intercept): the likelihood is flat along it, or too
custom skip, gamma-like            ll=-300.711233 conv=0 beta=0 WARN[Standard errors are not available for 1 of 2 parameters. b0_(Intercept): the likelihood is flat along it, or to
control silent                     ERROR `check_se` must be one of "warning", "ignore", "stop", not logical "FALSE"
```
<!-- nlfix:a2 end -->

- **Major 2, `conditional_effects()`.** `lp_eta_or_body()` at both sites, and
  `ord_prob_se()` reads the body's Jacobian from `frm_lp_basis()`, so an ordinal body now
  has a Wald band on the mean and category displays, equal to the linear spelling's to
  1e-3 of the largest se (test-nl-family.R). `band = "boot"` works. The refusal stays for a
  display that runs through several predictors (a zero-inflated mean), where boot,
  `method = "predict"` and `dpar = "mu"` all work, tested; it omits `method` on a category
  display, where it has no meaning.
- **Minor 1, cost.** Above: the guard stops at the first difference and reads the gradient
  first, `nl_domain_start()` reuses the objective built for the fit and skips a latent or
  category predictor, whose domain is the real line.
- **Minor 2, extension rows.** Section 5. **Minor 3.** Section 2, and NEWS and the row.
  **Minor 4.** Withdrawn in section 2. **Minor 5.** NEWS and `?frmtmb-autocor`, and
  `dev/upstream-bugs.md`.
- **Minor 6.** `bf(y ~ 1, autocor = ~ ar()) + nlf(mu ~ ...)` and `bf(y ~ 1) + nlf(mu ~ ...) +
  acformula(~ ar())` keep the term in the residual, equal to the `nl = TRUE` spelling; a
  right-hand side with anything else still stops, and its advice is now true.
- **Minor 7.** A helper `cs()` that the formula's environment finds is called again, with
  brms's `cs()` on the search path not counting as one.
- **Minor 8.** Value and derivative tests for all 20 names, against written closed forms,
  ungated, on the tape and on doubles (test-nl-stan-funs.R). `d_cloglog`, `d_lse`, `d_lde`,
  `d_hypot` and `d_Phi` are now killed ungated by construction: each compares `$jacobian`
  with the analytic form. **`ac_refit` is not a refit bug.** The real build keeps the term
  through `update(formula. =)` and a later `update(data =)` (-200.6437 in the review's
  log). -207.27 is the mutant losing it from the stored call, and also what a new formula
  without the term asks for, the same on a linear formula. A test now pins it.
- **Minor 9.** `vignettes/habit.Rmd` says what `nl = TRUE` does there.
- **Minor 10.** Section 5.
- **Minor 11.** `nl_pow()`: `x^b` and `pow(x, b)` with a data zero have Stan's value at every
  b and a gradient of 0 there, through `sign()`, and everything else is the same operation.
  L1 on the review's `r_powzero` body in the table below (`powzero`, `powzero_fun`).
- **Minor 12.** Section 7, with brms's evidence. **Minor 13.** Section 2.
- **NITs.** Search path: the help and NEWS say it. A column and a call of one name: resolved,
  tested. `hypot()` no longer overflows or underflows, `inv_square(0)` has -Inf, and
  `inv_cloglog` has 5e-63 at x = 5, each tested. The absolute tolerance in test-nl-family.R is
  a ratio. `acformula` twice: kept as brms has it, with a line in the help. `log(a)` in NEWS:
  reworded. `mixture(Gamma(), Gamma())` moves its start now. The mixture prior table: filed.

Default-start L1 specs (`sad_*`) are in section 2's table, and the rulings are recorded where
they apply: D7 upheld (section 2), the two `predictor_ac()` entries (section 2 and
`dev/upstream-bugs.md`), the segfault (section 2), p166 and tweedie (section 5), and
`files=279/289` is not a defect: the 10 are the `test-brms-suite-*.R` files, which source to
nothing in the ungated tier.

## 1. What changed

- **D1, the cause.** `R/frame.R`: the ordinal rule that thresholds replace
  the intercept tested `is_primary`, and `primary_dpars` also holds the
  nonlinear parameters a location's body reaches (for REML). So `b1 ~ 1`
  lost its one column, `b1` was 0 and the body was 0 on every row.
  `b1 ~ 1 + z` lost the intercept only, which no drop of the whole body
  would show. The rule now reads the family's own location. Seen on
  rellib-r16: `dev/nlfix-probe1.R`, `beta` empty, logLik -267.1211.
- **D1, what brms does, run and not recalled.**
  `dev/nlfix-brms-ordinal.R` (stancode, seed 49) and the L1 specs below
  (cached compile, log_prob): brms writes the body as `mu[n]`, keeps the
  thresholds as `Intercept`, and keeps the nlpar's intercept, under
  cumulative, sratio, cratio, acat and hurdle_cumulative, with and without
  `thres()` and `disc`. `cs(z)` in an nlpar formula is refused: "Category
  specific effects are only supported for the main parameter 'mu'". `cs(z)`
  in the body is written into Stan as a call to `cs()`, and the compile
  fails (pair p014). So the body is honored, and both `cs()` spellings are
  refused. frm() already refused the first in brms's words. The second has
  a new sentence.
- **D1, the post-fit paths.** `fitted()`, `predict()` and the rest were
  refused on an nl ordinal fit (nothing reached them, since the fit dropped
  the body). `ord_probs()`, `ord_mix_probs()`, `ord_dpars_more()` and
  `cat_probs()` now read the body through `lp_eta_or_body()`. The prior
  table lists the thresholds (`add_thres_rows()`), and the sampling table
  gives them `student_t(3, 0, 2.5)`, as `brms::default_prior()` does
  (`dev/nlfix-priors.R`).
- **D1, the audit.** The objective reads every body the same way:
  `dparv[[resp]][[dpar]]`. Whether a family honors it depends on two things:
  whether its design keeps the nlpar's coefficients (the bug) and whether its
  density reads that dpar. `dev/equiv-nl-sweep.R` measures the first on
  every core family (section 3), `dev/nlfix-extfam.R` and `dev/nlfix-asym.R`
  on cox, royston_parmar and wiener. The multi-mu families refused: see
  below. No other family drops the body.
- **The guard.** `nl_check_body_read()` (R/par-template.R), called from
  `fit_assembled()` once the tape exists, refuses (a) a parameter a body reads
  whose formula has terms and whose design has no column and (b) a body whose
  value moved while the objective did not (section 0 has the rule of fix round
  1, and the first form's false refusals). Present and absent cases are in
  test-nl-family.R: a custom family that skips mu refuses, the same family
  reading mu fits at the linear maximum; a hand-hollowed design refuses, the
  real frame passes; `b ~ 0` and `b ~ 0 + offset(o)` fit (brms builds
  `b ~ 0`). It fails open on a non-finite value or an objective that cannot
  run on doubles. The sweep below is the false-alarm run: no core family
  fires.
- **Multi-mu families, beyond the drop.** nl = TRUE on categorical,
  multinomial, dirichlet, dirichlet2, logistic_normal, dirichlet_multinomial
  and mixture was refused where brms builds and compiles them. brms writes
  the one body into each location with one set of nlpars
  (`dev/nlfix-brms-multimu.R`: `mub[n]` and `muc[n]` both read `nlp_b0`). The
  parser now copies the body the same way, and L1 is an identity. A family
  with several locations and no brms form (frmtmb.latent's `hmm()`) keeps the
  refusal, tested.
- **D7.** `bf(nl = TRUE, autocor =)` and `+ acformula()` keep the term beside
  the formula (slot `autocor`), since a body is evaluated and not read as
  terms, and the parser hands it to `pull_autocor()`. The objective already
  read the mean from the body. See section 2 for the measurement. brms
  refuses `ar(p = 2, cov = TRUE)` ("only possible for effects of maximal
  order one"), and frm() fits it, as on a linear formula (a decision of
  0.77.0, not this lane).
- **Default starts.** `nl_domain_start()`: only where the objective at the
  start is not finite and the body itself is, one intercept the body reads
  moves to where the body's mean equals `link(mean(y))`, by one linear solve
  per parameter, tried in the order the body names them; the first that makes
  the objective finite stays, with a message. No search of the objective (the
  reason `nl_start_message()` gives stands). Undefined bodies (`log(a)` at 0)
  and bodies with no start in the domain (`b1 * x` under an inverse link with
  x of both signs) keep the error. `start` and a located prior stand. On a
  multivariate model only the body that needs it moves.
- **Stan functions.** `nl_stan_funs` (R/ad-env.R). brms's vignettes and Rd
  examples call only `inv_logit()` in a body beyond base math:
<!-- nlfix:stanfuns begin -->
```
files read: 297 
calls found in nl bodies:

      exp inv_logit     nlfun      sqrt       log 
       21         5         2         2         1 

by name:
exp            vignettes/brms_multilevel.ltx: ult * (1 - exp(-(dev/theta)^omega)) ;; vignettes/brms_nonlinear.Rmd: b1 * exp(b2 * x) ;; man/brmsformula-helpers.Rd: a * exp(b * x)
inv_logit      vignettes/brms_nonlinear.Rmd: 0.33 + 0.67 * inv_logit(eta) ;; vignettes/brms_nonlinear.Rmd: guess + (1 - guess) * inv_logit(eta) ;; tests/testthat/tests.stancode.R: 
nlfun          tests/testthat/tests.stancode.R: nlfun(a, b, c, X) ;; tests/testthat/tests.stancode.R: nlfun(a, b, c, X2)
sqrt           man/brmsformula.Rd: tau * sqrt(eta)
log            tests/testthat/tests.stancode.R: log(eta) + a

brms numeric-helpers.R defines:
logit, inv_logit, cloglog, inv_cloglog, Phi, incgamma, square, cbrt, exp2, pow, inv, inv_sqrt, inv_square, hypot, log1m, step, int_step, logm1, expp1, logit_scaled, inv_logit_scale
```
<!-- nlfix:stanfuns end -->
  The wider set brms resolves on both of its sides is
  `brms:::names_stan_functions()`, 38 names, of which Stan 2.32 parses these
  with the argument types a row-wise body passes (`dev/nlfix-stanc.R`):
<!-- nlfix:stanc begin -->
```
in Stan: logit, inv_logit, inv_cloglog, Phi, square, cbrt, exp2, pow, inv, inv_sqrt, inv_square, hypot, log1m, step, log1p_exp, log1m_exp, log_diff_exp, log_sum_exp, log_inv_logit,
not in Stan: cloglog, incgamma, logm1, expp1, logit_scaled, inv_logit_scaled, multiply_log, log_mean_exp, log_expm1, scale_unit, fabs, inv_odds, softit, inv_softit, log_inv_softit,
```
<!-- nlfix:stanc end -->
  `softmax()` and `log_softmax()` take a vector only, so no row-wise brms body
  can call them; they are left out. The other 20 resolve in a body. A gated
  test evaluates each in Stan's generated quantities and asserts values within
  1e-12 relative (test-nl-stan-funs.R, "the values are Stan's"; the file's
  count with the gate on is in section 6). A function of the same name in the formula's environment
  is found first; `stats::step()` does not count as one.

## 2. Measurements against brms

L1 through the instrument: const is `log_prob - ours` at the estimate, extra
at the two stepped points, grad the max |gradient| (bound 1e-3).

<!-- nlfix:l1 begin -->
66 specs (`dev/nlfix-l1-specs.R`, through `dev/nlfix-l1.R`, the equiv instrument unchanged).

| spec | family | const | grad | extra | l1_ok | verdict | l1_error |
|---|---|---|---|---|---|---|---|
| ac_ar | gaussian()  | -1.14e-13 | 3.81e-05 | -1.14e-13 -2.84e-14 | TRUE | defect |  |
| ac_ar_cov | gaussian()  | 5.12e-13 | 2.18e-07 | 1.99e-13 -7.11e-13 | TRUE | equivalent |  |
| ac_arma | gaussian()  | 2.84e-14 | 4.72e-05 | 0.00 5.68e-14 | TRUE | defect |  |
| ac_arma_cov | gaussian()  | 3.98e-13 | 0.000864 | -1.08e-10 6.57e-10 | TRUE | equivalent |  |
| ac_cosy2 | gaussian()  | -5.68e-14 | 2.81e-05 | -3.13e-13 3.41e-13 | TRUE | equivalent |  |
| ac_fcor | gaussian()  | -5.68e-14 | 1.88e-06 | -1.14e-13 2.84e-14 | TRUE | equivalent |  |
| ac_ma | gaussian()  | 0 | 0.000178 | 8.53e-14 5.68e-14 | TRUE | defect |  |
| ac_ma_cov | gaussian()  | -1.71e-13 | 8.88e-05 | -5.97e-13 -4.55e-13 | TRUE | equivalent |  |
| ac_student | student()  | 2.56e-13 | 2.09e-05 | 3.13e-13 -4.26e-13 | TRUE | equivalent |  |
| ac_unstr | gaussian()  | 1.78e-14 | 0.000165 | 1.42e-14 -1.78e-14 | TRUE | cannot measure |  |
| fam_asym_laplace | asym_laplace()  | 0 | 2.72 | 0.00 -2.84e-14 | FALSE | defect |  |
| fun_guess | bernoulli('identity')  | -2.84e-14 | 2.13e-06 | 0 0 | TRUE | equivalent |  |
| fun_guess33 | bernoulli('identity')  | 0 | 3.14e-07 | 0.00 -1.42e-14 | TRUE | equivalent |  |
| mm_categorical | categorical()  | 0 | 0.000368 | 0 0 | TRUE | equivalent |  |
| mm_dirichlet | dirichlet()  | 5.68e-14 | 0.000443 | 2.84e-13 1.99e-13 | TRUE | equivalent |  |
| mm_dirichlet2b | brmsfamily('dirichlet2')  | -2.27e-13 | 4e-05 | -1.71e-13 -1.14e-13 | TRUE | equivalent |  |
| mm_dirmult | dirichlet_multinomial()  | 9.09e-13 | 4.51e-07 | 3.08e-11 -5.12e-12 | TRUE | equivalent |  |
| mm_logistic_normal | logistic_normal()  | -1.71e-13 | 0.000262 | -1.14e-13 -1.71e-13 | TRUE | equivalent |  |
| mm_mixture | mixture(gaussian(), gaussian())  | 0 | 0.000178 | 0.00 5.68e-14 | TRUE | equivalent |  |
| mm_multinomial | multinomial()  | -1.14e-12 | 0.000437 | -7.96e-13 -1.02e-12 | TRUE | equivalent |  |
| mv_nl | gaussian()  | 0 | 0.000735 | 1.14e-13 5.68e-14 | TRUE | equivalent |  |
| ord_acat_disc | acat()  | 5.68e-14 | 3.21e-05 | 0 0 | TRUE | equivalent |  |
| ord_acat_exp | acat()  | 0 | 6.72e-05 | 0 0 | TRUE | equivalent |  |
| ord_acat_nonlin | acat()  | -1.08e-12 | 0.000197 | 0.00 -7.28e-12 | TRUE | equivalent |  |
| ord_acat_plain | acat()  | 5.68e-14 | 1.21e-05 | 0 0 | TRUE | equivalent |  |
| ord_acat_re | acat()  | 0 | 7.82e-14 | 0 0 | TRUE | equivalent |  |
| ord_acat_thres | acat()  | 0 | 8.22e-06 | 0 0 | TRUE | equivalent |  |
| ord_cratio_disc | cratio()  | 0 | 0.000325 | 0 0 | TRUE | equivalent |  |
| ord_cratio_exp | cratio()  | 0 | 0.000124 | 0 0 | TRUE | equivalent |  |
| ord_cratio_nonlin | cratio()  | 0 | 9.78e-05 | 0 0 | TRUE | equivalent |  |
| ord_cratio_plain | cratio()  | 0 | 0.000343 | 0 0 | TRUE | equivalent |  |
| ord_cratio_re | cratio()  | 0 | 2.31e-12 | 0 0 | TRUE | equivalent |  |
| ord_cratio_thres | cratio()  | 0 | 0.00043 | 0 0 | TRUE | equivalent |  |
| ord_cumulative_disc | cumulative()  | -2.84e-14 | 0.000551 | 0.00 -5.68e-14 | TRUE | equivalent |  |
| ord_cumulative_exp | cumulative()  | 0 | 6.63e-05 | 0 0 | TRUE | equivalent |  |
| ord_cumulative_nonlin | cumulative()  | 0 | 0.0181 | 0.00 -1.09e-11 | FALSE | defect |  |
| ord_cumulative_plain | cumulative()  | 0 | 0.000522 | 0 0 | TRUE | equivalent |  |
| ord_cumulative_probit | cumulative(link = 'probit')  | 0 | 7.09e-05 | 0.00 5.68e-14 | TRUE | equivalent |  |
| ord_cumulative_re | cumulative()  | 0 | 8.83e-20 | 0 0 | TRUE | equivalent |  |
| ord_cumulative_thres | cumulative()  | 0 | 0.000293 | 0 0 | TRUE | equivalent |  |
| ord_equidistant | cumulative(threshold = 'equidistant')  | -5.68e-14 | 0.000285 | 0 0 | TRUE | equivalent |  |
| ord_grouped | cumulative()  | 0 | 0.000472 | 0 0 | TRUE | equivalent |  |
| ord_hurdle_hu | hurdle_cumulative()  | 0 | 0.000633 | 0 0 | TRUE | equivalent |  |
| ord_hurdle_plain | hurdle_cumulative()  | 0 | 0.000564 | 0 0 | TRUE | equivalent |  |
| ord_mixture | mixture(cumulative(), cumulative())  |  |  |  |  | equivalent | no brms design matrix X_mu1 for linear predictor mu1 |
| ord_mv | gaussian()  |  |  |  |  | equivalent | no brms design matrix X_y for linear predictor y |
| ord_sratio_disc | sratio()  | 0 | 8.94e-05 | 0 0 | TRUE | equivalent |  |
| ord_sratio_exp | sratio()  | 0 | 0.000124 | 0 0 | TRUE | equivalent |  |
| ord_sratio_nonlin | sratio()  | 0 | 0.0617 | -Inf -Inf | FALSE | defect |  |
| ord_sratio_plain | sratio()  | 0 | 2.33e-05 | -5.68e-14 -5.68e-14 | TRUE | equivalent |  |
| ord_sratio_re | sratio()  | 0 | 7.64e-10 | -5.68e-14 0.00 | TRUE | equivalent |  |
| ord_sratio_thres | sratio()  | 0 | 2.73e-05 | 0 0 | TRUE | equivalent |  |
| ord_sum_to_zero | cumulative(threshold = 'sum_to_zero')  | 0 | 0.000151 | 0 0 | TRUE | equivalent |  |
| powzero | gaussian()  | -2.13e-14 | 7.58e-07 | -3.55e-14 -5.68e-14 | TRUE | equivalent |  |
| powzero_fun | gaussian()  | -2.13e-14 | 7.58e-07 | -3.55e-14 -5.68e-14 | TRUE | equivalent |  |
| sad_acat | acat()  | 1.14e-13 | 1.91e-05 | 2.33e-10 1.46e-10 | TRUE | equivalent |  |
| sad_cratio | cratio()  | 0 | 2.62e-05 | 0 0 | TRUE | equivalent |  |
| sad_cumulative | cumulative()  | -1.14e-13 | 5.17e-05 | 0 0 | TRUE | equivalent |  |
| sad_gaussian | gaussian()  | 0 | 0.00021 | -5.68e-14 -5.68e-14 | TRUE | equivalent |  |
| sad_hurdle_cumulative | hurdle_cumulative()  | 0 | 3.71e-05 | 0 0 | TRUE | equivalent |  |
| sad_poisson | poisson()  | -5.68e-14 | 5.54e-06 | -3.41e-13 -1.71e-13 | TRUE | equivalent |  |
| sad_probit | cumulative('probit')  | 0 | 8.87e-05 | 0.00 5.68e-14 | TRUE | equivalent |  |
| sad_sat | cumulative()  | 0 | 1.04e-06 | -1.14e-13 5.68e-14 | TRUE | equivalent |  |
| sad_sratio | sratio()  | 0 | 2.62e-05 | 0 0 | TRUE | equivalent |  |
| start_gamma | Gamma()  | 2.84e-13 | 5.12e-06 | 1.42e-13 2.84e-14 | TRUE | equivalent |  |
| start_ig | inverse.gaussian()  | 8.53e-14 | 1.13e-05 | -1.14e-13 NaN | TRUE | equivalent |  |
<!-- nlfix:l1 end -->

Reading the table:

- The `ac_*` rows come with their methods run, so their verdict column carries the
  method half: `ac_ar`, `ac_ma` and `ac_arma` read "defect" on brms's post-processing
  (below), and `ac_unstr` "cannot measure" because its method half reached no level.
  `const`, `grad` and `l1_ok` are the model. `ar(p = 2, cov = TRUE)` has no row, as
  brms refuses it.
- `ord_cumulative_nonlin` and `ord_sratio_nonlin` read l1_ok FALSE on the
  gradient alone (const 0). On those data b1 runs off to about 3400 and b2 to
  1e-4, the linear limit of `b1 * exp(b2 * x)`, a ridge, and frm() warns of
  false convergence and of a flat direction. Corrected by the review (minor 13):
  "a restart stays" holds only near the ridge. On `ord_sratio_nonlin` the best of
  30 frm() starts is -249.3396, at b1 -0.13 and b2 -1.07, which brms's `optimizing`
  also finds (1 of 6), and that spec's `start` stops at -252.3176 with its warnings
  (`dev/review-nlfix-a9.R`). The `_exp` specs
  use a latent `1.5 exp(0.6 x)`, where the body is identified: const 0, gradient
  below the bound, for all four families.
- `fam_asym_laplace`: const 0 and extra at most 3e-14, gradient 2.72 at a kink
  of the density, where a gradient means nothing. Hence "conditional".
- `ord_mixture` and `ord_mv` are translator gaps in
  `stan_pars_from_fit()` (`no brms design matrix X_mu1`, `X_y`), and the
  instrument falls to L4 (`rstan::optimizing`, 4 and 6 starts): maxima differ
  by 2.25e-08 and 3.23e-10 nats. Not L1.
- `ord_hurdle_disc`: **cannot measure on brms's own program, which is a brms
  2.23.0 bug (ruling of the review).** `hurdle_cumulative_logit_lpmf` tests
  `y == nthres + 2` for the top category, which is `nthres + 1`, so every top row
  reads `thres[4]` and Stan stops, as `brms-1` of `dev/upstream-bugs.md` records. The
  instrument's process dies there on the linear spelling too (`ctl_hurdle_disc`).
  brms takes that function whenever `disc` is estimated. frm() is right: linear and
  nl logLik both -308.6792 (`dev/nlfix-segv.R`), and with the test patched to
  `nthres + 1` the log density equals frm()'s to 5.7e-14 on both spellings
  (`dev/review-nlfix-segv.R`). `ord_hurdle_hu` without disc is an identity. NEWS and
  the row say this now; my first text claimed `disc` for all five families.
- `start_ig` extra has a NaN: the stepped point left the link's domain on both
  sides.
- `ac_cosy` (no longer in the table, data with no correlation): rho was -0.033
  and brms's bound [0, 1] cannot hold it. `ac_cosy2`, with a real group effect,
  is the identity.

Methods at the same parameter values (L2; the first figure is the rule's
verdict, the second its metric: relative difference for fitted, the z ratio to
its bound for simulate and predict):

<!-- nlfix:l2 begin -->
35 specs. Each cell is ok (metric).

| spec | fitted | simulate | predict |
|---|---|---|---|
| ac_ar_cov_m | TRUE (0) | TRUE (0.639) | TRUE (0.627) |
| ac_ar_m | FALSE (0.0872) | TRUE (0.711) | FALSE (2.53) |
| ac_ar_m | FALSE (0.0872) | TRUE (0.711) | FALSE (2.53) |
| ac_arma_cov_m | TRUE (0) | TRUE (0.571) | TRUE (0.713) |
| ac_arma_m | FALSE (0.0337) | TRUE (0.673) | FALSE (1.23) |
| ac_arma_m | FALSE (0.0337) | TRUE (0.673) | FALSE (1.23) |
| ac_cosy2_m | TRUE (0) | TRUE (0.858) | TRUE (0.762) |
| ac_fcor_m | TRUE (0) | TRUE (0.768) | TRUE (0.671) |
| ac_ma_cov_m | TRUE (0) | TRUE (0.707) | TRUE (0.673) |
| ac_ma_m | FALSE (0.0634) | TRUE (0.845) | FALSE (2) |
| ac_ma_m | FALSE (0.0634) | TRUE (0.845) | FALSE (2) |
| ac_student_m | TRUE (0) | TRUE (0.741) | TRUE (0.643) |
| ac_unstr_m | TRUE (0) | NA (NA) | NA (NA) |
| ctl_dirichlet_m | TRUE (1.11e-16) | FALSE (NA) |  |
| ctl_hurdle_m | TRUE (1.67e-16) | FALSE (7.67) |  |
| ctl_multinomial_m | TRUE (4.87e-16) | FALSE (NA) |  |
| ctl_mv_m | TRUE (1.99e-16) | FALSE (NA) |  |
| fun_guess_m | TRUE (0) | TRUE (0.565) | TRUE (0.557) |
| fun_guess33_m | TRUE (0) | TRUE (0.621) | TRUE (0.534) |
| mm_categorical_m | TRUE (1.11e-16) | TRUE (0.679) | TRUE (0.742) |
| mm_dirichlet_m | TRUE (1.11e-16) | FALSE (NA) | TRUE (0.788) |
| mm_dirmult_m | TRUE (4.22e-16) | FALSE (NA) | TRUE (0.807) |
| mm_logistic_normal_fit | TRUE (NA) |  |  |
| mm_mixture_m | TRUE (0) | TRUE (0.568) | TRUE (0.635) |
| mm_multinomial_m | TRUE (4.32e-16) | FALSE (NA) | TRUE (0.667) |
| mv_nl_m | TRUE (2.21e-16) | FALSE (NA) | TRUE (0.63) |
| ord_acat_plain_m | TRUE (2.22e-16) | TRUE (0.557) | TRUE (0.616) |
| ord_cratio_plain_m | TRUE (1.11e-16) | TRUE (0.723) | TRUE (0.572) |
| ord_cumulative_disc_m | TRUE (2.08e-16) | TRUE (0.764) | TRUE (0.692) |
| ord_cumulative_plain_m | TRUE (2.22e-16) | TRUE (0.53) | TRUE (0.683) |
| ord_cumulative_re_m | TRUE (1.94e-16) | TRUE (0.674) | TRUE (0.797) |
| ord_hurdle_plain_m | TRUE (1.67e-16) | FALSE (7.67) | TRUE (0.671) |
| ord_sratio_plain_m | TRUE (1.67e-16) | TRUE (0.732) | TRUE (0.713) |
| start_gamma_m | TRUE (0) | TRUE (0.63) | TRUE (0.724) |
| start_ig_m | TRUE (0) | TRUE (0.581) | TRUE (0.661) |
<!-- nlfix:l2 end -->

Reading it:

- **`ac_ar_m`, `ac_ma_m`, `ac_arma_m`, `fitted` and `predict` FALSE: this is a
  brms 2.23.0 defect, not frm()'s.** `predictor.bprepnl`
  (`brms/R/predictor.R`) has no `predictor_ac()`, so on an nl model with
  `cov = FALSE` brms's R side drops the ARMA shift its own Stan program has.
  Measured on one model written both ways (`dev/nlfix-arnl2.R`, seed 1582):
<!-- nlfix:arnl begin -->
```
nl fitted - body, sd: 0.03465509 
simulate nl: slope of (mean - body) on the shift: -0.01165141 
simulate linear: slope of (mean - linear body) on the shift: -0.01165142 
predict nl: slope: 0.9984374 
predict linear: slope: 0.9946748 
brms nl epred - body, max: 0 
brms linear epred - frm nl fitted, max: 2.602378e-06 
brms linear epred - frm linear fitted, max: 0 
brms log_lik sum, nl and linear: -149.1356 -148.6737 
frm logLik: -148.6737 
```
<!-- nlfix:arnl end -->
  brms's nl `posterior_epred` is the body alone (difference 0), and its nl
  `log_lik` sums to -149.1356 where its own linear spelling and its Stan
  `log_prob` give -148.6737. frm() on the nl spelling equals frm() and brms on
  the linear spelling, shift included. The project tiebreaker applies: brms is
  clearly wrong against its own program, so frm() does not copy this. **Ruling of
  the round: upheld.** The review's own model has phi 0.531: brms's nl `log_lik` sums
  to -185.87332294, the density at phi = 0, against -150.99326554 from its Stan
  `log_prob`, its linear program, its linear `log_lik` and frm(), 34.9 nats. It is a
  divergence by decision under section 2.6, with the reason in NEWS and in
  `?frmtmb-autocor`. There are two upstream entries for `predictor_ac()`: latres's
  brms-45, where the latent arm misses the AR sum, and this one, where
  `predictor.bprepnl` calls neither arm. Filed in `dev/upstream-bugs.md`, number to
  be settled at consolidation. Nothing is posted upstream.
- **`simulate` NA (`mm_multinomial_m`, `mm_dirichlet_m`, `mm_dirmult_m`,
  `mv_nl_m`) and 7.67 (`ord_hurdle_plain_m`):** instrument. The same on linear
  formulas: `ctl_multinomial_m`, `ctl_dirichlet_m`, `ctl_mv_m` NA (draw shapes
  the instrument does not pair) and `ctl_hurdle_m` 7.67, max |z| 36.13: frm()
  simulates a factor whose codes run from 1, brms codes from 0.
- `mm_logistic_normal_m` did not finish its `simulate` and `predict` in 30
  minutes and was stopped; `mm_logistic_normal_fit` reads fitted alone.
- `simulate()` in sample of a `cov = FALSE` model: slope of the row means on the
  observed shift -0.012 on both spellings, while `predict()` carries it (0.998 and
  0.995). **This is correct and was wrongly filed in the first version of these
  findings.** A replicate has no observed past, so its row mean is the body, and it
  runs the recursion over its own draws (review, minor 4: row sd 0.6025 against the
  marginal 0.6075, lag-1 correlation 0.508 against 0.520 in the data, mean AR over 40
  refits on replicates 0.514). Withdrawn. The defect was on `newdata` (blocker 2).

Each route's maxima against the linear fit:

<!-- nlfix:d1check begin -->
```
cumulative         linear -255.5164 nl -255.5164 (diff 0.00e+00) thresholds-only -267.1211; b0+b1x -255.5164; nonlinear -255.4415; disc nl -255.5051, linear -255.5051
sratio             linear -255.9867 nl -255.9867 (diff 0.00e+00) thresholds-only -267.1211; b0+b1x -255.9867; nonlinear -255.9304; disc nl -255.9552, linear -255.9552
cratio             linear -255.9867 nl -255.9867 (diff 0.00e+00) thresholds-only -267.1211; b0+b1x -255.9867; nonlinear -255.9304; disc nl -255.9552, linear -255.9552
acat               linear -255.7561 nl -255.7561 (diff 0.00e+00) thresholds-only -267.1211; b0+b1x -255.7561; nonlinear -255.6704; disc nl -255.0796, linear -255.0796
hurdle_cumulative  linear -299.7922 nl -299.7922 (diff 0.00e+00) thresholds-only -311.6032; b0+b1x -299.7922; nonlinear -299.2921; disc nl -299.7208, linear -299.7208
cumulative         thres(3): linear -255.5164 nl -255.5164
sratio             thres(3): linear -255.9867 nl -255.9867
cratio             thres(3): linear -255.9867 nl -255.9867
acat               thres(3): linear -255.7561 nl -255.7561
```
<!-- nlfix:d1check end -->

Every core family's default link, and the bounded non-default links, from the
default starts (`dev/nlfix-startcheck.R`, seed 20261029):

<!-- nlfix:starts begin -->
```
Gamma()                                linear  -296.7449 nl  -296.7449 diff  0.00e+00 
Gamma(link = 'inverse')                linear  -296.7449 nl  -296.7449 diff  0.00e+00 
inverse.gaussian()                     linear  -311.0008 nl  -311.0008 diff  0.00e+00 
Gamma(link = 'identity')               linear  -300.0357 nl  -300.0357 diff  0.00e+00 
Gamma(link = 'sqrt')                   linear         NA nl         NA diff        NA ERR 'sqrt' is not a supported link for family 'Gamma'. Supported links are: 'log', 'identity',
Gamma(link = 'log')                    linear  -293.9139 nl  -293.9139 diff -1.52e-09 
inverse.gaussian(link = 'log')         linear  -301.4470 nl  -301.4470 diff -2.05e-11 
inverse.gaussian(link = 'inverse')     linear  -302.9087 nl  -302.9087 diff  0.00e+00 
poisson(link = 'identity')             linear  -487.8291 nl  -487.8291 diff  4.72e-12 
poisson(link = 'sqrt')                 linear  -486.3169 nl  -486.3169 diff -1.78e-11 
lognormal()                            linear  -290.6830 nl  -290.6830 diff -7.39e-12 
weibull()                              linear  -311.6559 nl  -311.6559 diff -1.22e-11 
exponential()                          linear  -450.4048 nl  -450.4048 diff -1.52e-09 
frechet()                              linear  -314.6481 nl  -314.6481 diff -6.10e-11 
hurdle_gamma()                         linear  -293.9139 nl  -293.9139 diff  8.65e-09 
-296.7449 msg:  
-297.6888 msg: Nonlinear starting values placed at the prior locations: b0_(Intercept) = 0.4. Give `start` to choose your own.
-739.2238  msg: Nonlinear starting values moved to where the body equals the response mean, because the link cannot evaluate it at 0: y_b0_(Intercept) = 0.5793. Give `start` to cho
sum of the two linear fits: -739.2238 
```
<!-- nlfix:starts end -->

Extension families (`dev/nlfix-extfam.R`, `dev/nlfix-asym.R`, seed 7101) and
`sar()` (`dev/nlfix-sar.R`, seed 9101). `b1 * x` under cox has no level, as in
brms, so its row in the first script is another model; the `b0 + b1 * x` lines
are the identity:

<!-- nlfix:ext begin -->
```
cox                    linear  -298.9452 nl  -406.5825 intercept  -323.6844 
royston_parmar         linear  -299.4877 nl  -299.9709 intercept  -323.9121 
wiener                 linear  -228.5402 nl  -228.5402 intercept  -242.6380 
custom reads mu        linear  -454.0095 nl  -454.0095 intercept  -481.2498 
custom skips mu        linear  -482.6647 nl         NA intercept  -482.6647 nl ERR: The nonlinear body of 'mu' does not reach the likelihood of family 'nlfixskip': moving the coeff
asym_laplace: linear -451.4893 nl default start -451.4871 nl from linear -451.4900
cox linear -298.9452 nl with b0 -298.9452 
rp linear -299.4877 nl with b0 -299.4877 
logLik linear, nl: -91.84050596 -91.84050596 
max |fitted difference|: 1.032745e-06 
brms stancode: 2626 
```
<!-- nlfix:ext end -->

### Sampler route (item 6)

frm_sample() builds from core's frame and objective, so D1 and D7 were there too
and the core fix reaches them; no code in frmtmb.sample changes.
`dev/nlfix-sample.R`, seeds 49 and 158, 1 chain, 500 iterations.
Before (rellib-r16, `dev/nlfix-sample-before.log`): variables
`b_Intercept[1:3], lp__`, no b1, `posterior_epred` refused. After: `b_b1_Intercept`
posterior mean 0.675, epred 250 x 200 x 4,
the prior summary has brms's `student_t(3, 0, 2.5)`; fcor and `ar(cov = TRUE)`
with a body sample (`ar[1]` in the draws), and categorical with one body reads
epred 200 x 150 x 3. `test-nl-body-draws.R` in frmtmb.sample holds these, and the
cov = FALSE shift by hand against the draws.

## 3. The nl x family sweep, rerun

`FRMTMB_EQUIV_LIB=<lane> Rscript dev/equiv-nl-sweep.R dev/nlfix-equiv-out`, not
recorded. The record (`dev/equiv-data/nl-sweep.tsv`, 0.77.0) had nl honored 35,
differs 1, dropped 5, refuses 6.

<!-- nlfix:sweep begin -->
| outcome | families |
|---|---|
| 1 | 7 |
| 2 | 41 |

| family | status_now | frm_linear | frm_nl | frm_intercept | outcome | proposal | default_linear | default_nl_default_start |
|---|---|---|---|---|---|---|---|---|
| gaussian | works | -422.597 | -422.597 | -451.955 | nl honored |  |  |  |
| student | works | -425.3746 | -425.3746 | -454.8267 | nl honored |  |  |  |
| lognormal | works | -290.683 | -290.683 | -352.5927 | nl honored |  |  |  |
| shifted_lognormal | works | -290.3482 | -290.3482 | -352.2866 | nl honored |  |  |  |
| skew_normal | works | -417.174 | -417.174 | -449.2481 | nl honored |  |  |  |
| exgaussian | works | -300.8931 | -300.8931 | -369.6141 | nl honored |  |  |  |
| asym_laplace | conditional | -446.1533 | -446.9826 | -469.3713 | nl differs from both |  |  |  |
| Gamma | works | -293.9139 | -293.9139 | -357.3398 | nl honored |  | -296.7449 | -296.7449 |
| weibull | works | -311.6559 | -311.6559 | -373.8665 | nl honored |  |  |  |
| exponential | works | -450.4048 | -450.4048 | -463.7779 | nl honored |  |  |  |
| inverse.gaussian | works | -301.447 | -301.447 | -354.1034 | nl honored |  | -311.0008 | -311.0008 |
| beta | works | 235.5179 | 235.5179 | 130.4089 | nl honored |  |  |  |
| tweedie | works | -425.1249 | -425.1249 | -443.1407 | nl honored |  |  |  |
| poisson | works | -442.6096 | -442.6096 | -488.0239 | nl honored |  |  |  |
| negbinomial | works | -442.6096 | -442.6096 | -483.564 | nl honored |  |  |  |
| negbinomial2 | works | -442.6096 | -442.6096 | -483.564 | nl honored |  |  |  |
| nbinom1 | works | -442.6096 | -442.6096 | -483.564 | nl honored |  |  |  |
| geometric | works | -490.2892 | -490.2892 | -507.8078 | nl honored |  |  |  |
| compois | works | -442.4324 | -442.4324 | -483.6511 | nl honored |  |  |  |
| binomial | works | -538.5875 | -538.5875 | -621.0102 | nl honored |  |  |  |
| bernoulli | works | -191.7266 | -191.7266 | -201.9035 | nl honored |  |  |  |
| beta_binomial | works | -546.2155 | -546.2155 | -605.3154 | nl honored |  |  |  |
| zero_inflated_beta_binomial | works | -583.4681 | -583.4681 | -630.0495 | nl honored |  |  |  |
| multinomial | works | -1023.5413 | -1083.3188 | -1088.7264 | nl differs from both |  |  |  |
| dirichlet_multinomial | works | -1041.0492 | -1092.0275 | -1103.3302 | nl differs from both |  |  |  |
| zero_inflated_poisson | works | -416.6716 | -416.6716 | -452.4958 | nl honored |  |  |  |
| zero_inflated_negbinomial | works | -416.6716 | -416.6716 | -450.2649 | nl honored |  |  |  |
| zero_inflated_binomial | works | -573.8408 | -573.8408 | -628.6198 | nl honored |  |  |  |
| zero_inflated_beta | works | 35.7166 | 35.7166 | -60.4096 | nl honored |  |  |  |
| zero_one_inflated_beta | works | -4.0226 | -4.0226 | -78.4173 | nl honored |  |  |  |
| xbeta | works | 248.0806 | 248.0806 | 138.2045 | nl honored |  |  |  |
| hurdle_poisson | works | -506.908 | -506.908 | -541.1679 | nl honored |  |  |  |
| hurdle_negbinomial | works | -506.908 | -506.908 | -540.4004 | nl honored |  |  |  |
| hurdle_gamma | works | -394.4969 | -394.4969 | -433.8501 | nl honored |  |  |  |
| hurdle_lognormal | works | -393.0446 | -393.0446 | -432.6877 | nl honored |  |  |  |
| cumulative | works | -395.9068 | -395.9068 | -411.3239 | nl honored |  |  |  |
| hurdle_cumulative | works | -462.6041 | -462.6041 | -477.4804 | nl honored |  |  |  |
| sratio | works | -397.7524 | -397.7524 | -411.3239 | nl honored |  |  |  |
| cratio | works | -397.7524 | -397.7524 | -411.3239 | nl honored |  |  |  |
| acat | works | -396.2385 | -396.2385 | -411.3239 | nl honored |  |  |  |
| categorical | works | -328.7267 | -329.3203 | -328.964 | nl differs from both |  |  |  |
| von_mises | works | -228.3185 | -228.3185 | -334.9075 | nl honored |  |  |  |
| dirichlet | works | 275.5537 | 267.8864 | 240.1086 | nl differs from both |  |  |  |
| dirichlet2 | works | 275.2745 | 236.2029 | 245.9485 | nl differs from both |  |  |  |
| logistic_normal | works | 287.1736 | 269.6627 | 259.8158 | nl differs from both |  |  |  |
| frechet | works | -314.6481 | -314.6481 | -384.7677 | nl honored |  |  |  |
| gen_extreme_value | works | -431.1333 | -431.1333 | -467.9903 | nl honored |  |  |  |
| discrete_weibull | works | -445.9729 | -445.9729 | -483.6242 | nl honored |  |  |  |
<!-- nlfix:sweep end -->

Reading it: `dev/equiv-nl-sweep.R` compares a body with the linear fit, which
fits only the single-location families. On categorical, multinomial,
dirichlet_multinomial, dirichlet, dirichlet2 and logistic_normal "nl differs from
both" is the model: brms's one body with shared parameters is not a linear
predictor per location, so its maximum sits below theirs by construction, and the
evidence for those rows is L1 in section 2. `asym_laplace` is the kink. No family
reads "nl body DROPPED", and no core family stops at the guard.

## 4. Equiv pairs rerun against the baseline

The 34 pairs `dev/nlfix-equiv-ids.R` says the fix can reach (a body, an ordinal
family or `autocor =`), remeasured; the other 138 keep the baseline's `<id>.tsv`.

<!-- nlfix:verdicts begin -->
35 pairs rerun, 2 changed verdict.

Changed:

| id | pair | before | before_level | after | after_level | l1_const | l1_grad | frm_error |
|---|---|---|---|---|---|---|---|---|
| p049 | thres() x nl | defect | L4 | equivalent | L1 | 0 | 0.000198 |  |
| p158 | fcor() x nl | defect |  | compat row wrong | L1 | -5.68e-14 | 1.88e-06 |  |

All rerun pairs:

| id | pair | before | before_level | after | after_level | l1_const | l1_grad | frm_error |
|---|---|---|---|---|---|---|---|---|
| p001 | prior x nl | equivalent | L1 | equivalent | L1 | -2.84e-14 | 1.87e-05 |  |
| p010 | cens() x disc | compat row wrong |  | compat row wrong |  | NA | NA | cens()/trunc() need a family with a CDF, and 'cumulative' su |
| p011 | cens() x equidistant | compat row wrong |  | compat row wrong |  | NA | NA | cens()/trunc() need a family with a CDF, and 'cumulative' su |
| p013 | cs_pred() x frm_curve | defect | L1+L2 | defect | L1+L2 | 0 | 0.000224 |  |
| p014 | cs_pred() x nl | compat row wrong |  | compat row wrong |  | NA | NA | cs() is not supported in a nonlinear formula body: b1 * x +  |
| p015 | cs_pred() x prior | equivalent | L1 | equivalent | L1 | 0 | 0.00056 |  |
| p016 | cs_pred() x simulate | equivalent | L1+L2 | equivalent | L1+L2 | 0 | 0.000497 |  |
| p019 | equidistant x frm_curve | equivalent | L1+L2 | equivalent | L1+L2 | 0 | 0.000372 |  |
| p027 | importance x equidistant | equivalent | L1+L4 | equivalent | L1+L4 | 0 | 1.78e-14 |  |
| p031 | nl x hypothesis_profile | equivalent | L1+L2+L4 | equivalent | L1+L2+L4 | -5.68e-14 | 1.34e-07 |  |
| p032 | nl x simulate | equivalent | L1+L2 | equivalent | L1+L2 | 0 | 7.29e-05 |  |
| p042 | s() x nl | equivalent | L4 | equivalent | L4 | NA | NA |  |
| p047 | t2() x nl | equivalent | L4 | equivalent | L4 | NA | NA |  |
| p049 | thres() x nl | defect | L4 | equivalent | L1 | 0 | 0.000198 |  |
| p059 | vint() x sum_to_zero | compat row wrong |  | compat row wrong |  | NA | NA | cumulative: the addition term `vint()` is not one this famil |
| p064 | cumulative x frm_curve | equivalent | L1+L2 | equivalent | L1+L2 | 0 | 0.000258 |  |
| p075 | cumulative x hypothesis_profile | equivalent | L1+L2+L4 | equivalent | L1+L2+L4 | 0 | 0.000196 |  |
| p076 | cumulative x smooth | equivalent | L1 | equivalent | L1 | 0 | 1.27e-19 |  |
| p093 | hurdle_cumulative x residuals | equivalent | L1+L2 | equivalent | L1+L2 | 0 | 0.000353 |  |
| p107 | sratio x frm_curve | equivalent | L1+L2 | equivalent | L1+L2 | 0 | 0.000809 |  |
| p111 | sum_to_zero x mmc() | equivalent | L1 | equivalent | L1 | -6.82e-11 | 4.34e-11 |  |
| p119 | gaussian x nl | equivalent | L1 | equivalent | L1 | 0 | 6.71e-07 |  |
| p123 | beta_binomial x verbose | equivalent | L1 | equivalent | L1 | -1.59e-12 | 0.00019 |  |
| p127 | cens() x nl | equivalent | L1 | equivalent | L1 | -1.14e-13 | 2.45e-06 |  |
| p129 | cratio x mvbf | equivalent | L4 | equivalent | L4 | NA | NA |  |
| p131 | cs_pred() x mixture | defect | L4 | defect | L4 | NA | NA |  |
| p136 | gp x nl | equivalent | L4 | equivalent | L4 | -152 | 341 |  |
| p143 | smooth x nl | equivalent | L4 | equivalent | L4 | NA | NA |  |
| p144 | sratio x prior | equivalent | L1 | equivalent | L1 | 0 | 0.000479 |  |
| p145 | thres() x prior | equivalent | L1 | equivalent | L1 | -0.107 | 0.000146 |  |
| p153 | acat x rescor | refusal correct |  | refusal correct |  | NA | NA | rescor = TRUE requires all responses to be gaussian, or all  |
| p155 | cumulative x bar_crossing | refusal correct |  | refusal correct |  | NA | NA | A random-effect term cannot be crossed with '*': x * (1 \| g) |
| p158 | fcor() x nl | defect |  | compat row wrong | L1 | -5.68e-14 | 1.88e-06 |  |
| p166 | nbinom1 x nl | compat row wrong |  | compat row wrong |  | NA | NA |  |
| p171 | zero_inflated_beta x nl | compat row wrong | L1 | compat row wrong | L1 | 0 | 1.96e-05 |  |
<!-- nlfix:verdicts end -->

Reading it:

- p049 is the D1 pair: defect at L4 to equivalent at L1, const 0, gradient 2e-4.
- p158 is the D7 pair: L1 const -5.68e-14, gradient 1.88e-6. It reads "compat row
  wrong" and not "equivalent" because the flow compares with the status frozen in
  `dev/equiv-data/sample.tsv` ("refused"), which the lane may not edit; the table
  now says works. The same holds for p166 and p171 (measured L1 const 0 on p171),
  frozen at "refused", and for p014, frozen at "untested", where both sides refuse
  and the table now says refused. Their rows are fixed and their verdicts follow at
  consolidation.
- `summary.md` of this outdir does not move the headline, because
  `dev/equiv-data/adjudication.tsv` sets p049, p158 and p166 by hand with the 0.77.0
  reasons, and the summary reads the adjudication over the measurement. **Those
  three adjudication rows must go at consolidation** (the record is rerun once then).
  The per-pair `.tsv` above are the measurement.
- p013 and p131 are other lanes' or filed defects, and stay.
- p136's const -152 is the exact gp() gap, scored at L4 as before (lane gpgap).

## 5. Compat rows changed

Every rule is on its own line in R/compat.R.

- Removed: `nl x kind:family`, refused (the blanket rule).
- `nl x` student, lognormal, shifted_lognormal, skew_normal, exgaussian, weibull,
  exponential, beta, poisson, negbinomial, negbinomial2, geometric, compois,
  binomial, bernoulli, beta_binomial, zero_inflated_beta_binomial,
  zero_inflated_poisson, zero_inflated_negbinomial, zero_inflated_binomial,
  zero_inflated_beta (p171), zero_one_inflated_beta, xbeta, hurdle_poisson,
  hurdle_negbinomial, hurdle_gamma, hurdle_lognormal, von_mises, frechet,
  gen_extreme_value, discrete_weibull: refused to **works** (the sweep's proposals).
- `nl x` tweedie, nbinom1 (p166): refused to **works**, and each note opens "Works,
  no brms reference", as the survey proposed. Section 2.6 makes them cannot-measure at
  consolidation (review, p166).
- `nl x` Gamma, inverse.gaussian: refused to **works**. The proposal was
  "conditional (default link needs starting values)"; the start fix removes the
  condition.
- `nl x` cumulative, hurdle_cumulative, sratio, cratio, acat: refused to **works**
  (proposal was "broken", D1; the fix removes it).
- `nl x` categorical, multinomial, dirichlet_multinomial, dirichlet, dirichlet2,
  logistic_normal: refused to **works**; not among the proposals, decided by the
  measurement above.
- `nl x asym_laplace`: refused to **conditional**; the proposal was "untested".
- `thres() x nl` (p049): untested to **works**.
- `nl x` disc, equidistant, sum_to_zero, mixture, mvbf: untested to **works**.
- `cs_pred() x nl` (p014): untested to **refused**.
- `group:autocor x nl`: refused to **conditional** (fix round 1; **works** in the
  first version): gaussian and student work, a family whose term would be a latent
  residual stays refused until latres's code says otherwise. Six rows, p158 among
  them. On its own line, since it is the one text conflict of the trial merges.
- Extension rows, each on its own line in that package's rules (fix round 1, minor
  2): `nl x` cox, royston_parmar, wiener and wiener_gng **works**, measured; `nl x` rdm
  and lba **refused**, naming `nlf(v2 ~ ...)` as the route; `nl x` the frmtmb.learn
  families **refused**, untested before, naming `nlf()`. gddm reads "untested" through
  the table's default: not measured. `sar() x nl` stays untested in frmtmb.spatial's
  rules, measured honored above.

## 6. Tests

New files: `tests/testthat/test-nl-family.R` (ungated),
`tests/testthat/test-nl-stan-funs.R` (ungated, one gated test),
`tests/testthat/test-nl-brms.R` (gated),
`extensions/frmtmb.sample/tests/testthat/test-nl-body-draws.R`. One existing test
changed: `tests/testthat/test-conditions-census.R` and its generator
`dev/conditions-census-gen.R` gain three exemptions, the rethrows of a held warning or
error in `nl_fit_off_saddle()`. `tests/testthat/helper-brms.R` moved four lines of the
translator so an ordinal body, which has no design for mu, reaches its threshold rule.
Runner: `dev/nlfix-testrun.R`, `dev/nlfix-suite.sh`, `dev/nlfix-final3.sh`.

Seen to fail on rellib-r16, behaviorally: a body under each ordinal family reads the
intercept-only maximum, and `expect_lt(abs(nl - lin) / (lin - int), 1e-3)` fails;
`b_b1_Intercept` is absent from the sampler's draws; inv_logit() stops with "could
not find function"; the frm_compat() rows read refused. The `nl_stan_funs` tests'
missing symbol is the weak form, and their fit tests are the behavioral one. The
test-nl-brms errors on r16 are the translator's domain error on the missing b1
(D1), the refused term (D7), the failed start and the missing function.

<!-- nlfix:suite_r16 begin -->
```
Before, the new test files on rellib-r16 (0.78.0):
FILE test-nl-family.R PASS=102 FAIL=51 ERROR=26 SKIP=0 WARN=2
FILE test-nl-stan-funs.R PASS=2 FAIL=0 ERROR=8 SKIP=1 WARN=0
FILE test-nl-brms.R PASS=2 FAIL=1 ERROR=9 SKIP=0 WARN=0
FILE test-nl-body-draws.R PASS=1 FAIL=1 ERROR=3 SKIP=0 WARN=0
Behavioral failures among them, first lines:
[the nl body reaches the likelihood of cumulative] expectation_failure: Expected `abs(nl - lin)/(lin - int)` < 0.001.
[the nl body reaches the likelihood of cumulative] expectation_failure: Expected `st$status[st$feature_b == f]` to be identical to "works".
Error: The nonlinear formula body could not be evaluated: could not find function "inv_logit"
[a body can call inv_logit(), the vignette's item response model] expectation_error: <frmtmb_error/error/condition>
Error: The nonlinear formula body could not be evaluated: could not find function "inv_logit"
[an nl body under an ordinal family is sampled] expectation_failure: Expected `"b_b1_Intercept" %in% variables(s)` to be TRUE.
[the nl body reaches the likelihood of cumulative] expectation_failure: Expected `abs(nl - lin)/(lin - int)` < 0.001.
Actual comparison: 1.000 >= 0.001
Difference: 0.999 >= 0
[the nl body reaches the likelihood of cumulative] expectation_failure: Expected `st$status[st$feature_b == f]` to be identical to "works".
```
<!-- nlfix:suite_r16 end -->

<!-- nlfix:suite_core begin -->
```
Core, ungated, one file per process, lane library:
TOTAL files=290 PASS=20980 FAIL=0 ERROR=0 SKIP=217 WARN=0
FILE test-nl-body-vars.R PASS=23 FAIL=0 ERROR=0 SKIP=0 WARN=0
FILE test-nl-brms.R PASS=0 FAIL=0 ERROR=0 SKIP=10 WARN=0
FILE test-nl-family.R PASS=321 FAIL=0 ERROR=0 SKIP=0 WARN=0
FILE test-nl-lexical.R PASS=31 FAIL=0 ERROR=0 SKIP=0 WARN=0
FILE test-nl-rtmb-scope.R PASS=105 FAIL=0 ERROR=0 SKIP=0 WARN=0
FILE test-nl-stan-funs.R PASS=122 FAIL=0 ERROR=0 SKIP=1 WARN=0
FILE test-nl.R PASS=92 FAIL=0 ERROR=0 SKIP=0 WARN=0
FILE test-nlf.R PASS=77 FAIL=0 ERROR=0 SKIP=0 WARN=0
```
<!-- nlfix:suite_core end -->

<!-- nlfix:suite_ext begin -->
```
Extensions, ungated, one file per process, lane core, base extensions:
TOTAL files=183 PASS=8861 FAIL=0 ERROR=0 SKIP=47 WARN=0
frmtmb.coupling files=11 PASS=551 FAIL=0 ERROR=0 SKIP=5
frmtmb.eam files=32 PASS=1838 FAIL=0 ERROR=0 SKIP=5
frmtmb.latent files=12 PASS=681 FAIL=0 ERROR=0 SKIP=2
frmtmb.learn files=15 PASS=430 FAIL=0 ERROR=0 SKIP=13
frmtmb.ode files=11 PASS=549 FAIL=0 ERROR=0 SKIP=1
frmtmb.sample files=79 PASS=3782 FAIL=0 ERROR=0 SKIP=15
frmtmb.spatial files=8 PASS=478 FAIL=0 ERROR=0 SKIP=6
frmtmb.survival files=15 PASS=552 FAIL=0 ERROR=0 SKIP=0
FILE test-nl-body-draws.R PASS=15 FAIL=0 ERROR=0 SKIP=0 WARN=0
```
<!-- nlfix:suite_ext end -->

<!-- nlfix:suite_gated begin -->
```
Gated (FRMTMB_BRMS_FIT_TESTS=true), the files that skipped ungated:
TOTAL files=35 PASS=3040 FAIL=0 ERROR=0 SKIP=15 WARN=0
FILE test-curve-scale.R PASS=0 FAIL=0 ERROR=0 SKIP=1 WARN=0
FILE test-drmtmb-agreement.R PASS=0 FAIL=0 ERROR=0 SKIP=13 WARN=0
FILE test-fuzz.R PASS=0 FAIL=0 ERROR=0 SKIP=1 WARN=0
FILE test-nl-brms.R PASS=91 FAIL=0 ERROR=0 SKIP=0 WARN=0
FILE test-nl-stan-funs.R PASS=162 FAIL=0 ERROR=0 SKIP=0 WARN=0
frmtmb.sample, gated:
TOTAL files=7 PASS=731 FAIL=0 ERROR=0 SKIP=1 WARN=0
OpenBLAS 0.3.26 and 0.3.32 with LAPACK, new test files (gated for test-nl-brms):
ob0.3.26-test-beta-binomial-density.txt:FILE test-beta-binomial-density.R PASS=93 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.26-test-cdfs.txt:FILE test-cdfs.R PASS=234 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.26-test-cens-trunc.txt:FILE test-cens-trunc.R PASS=72 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.26-test-nl-brms.txt:FILE test-nl-brms.R PASS=91 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.26-test-nl-family.txt:FILE test-nl-family.R PASS=321 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.26-test-nl-stan-funs.txt:FILE test-nl-stan-funs.R PASS=122 FAIL=0 ERROR=0 SKIP=1 WARN=0
ob0.3.26-test-numerical-robustness.txt:FILE test-numerical-robustness.R PASS=723 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.32-test-beta-binomial-density.txt:FILE test-beta-binomial-density.R PASS=93 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.32-test-cdfs.txt:FILE test-cdfs.R PASS=234 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.32-test-cens-trunc.txt:FILE test-cens-trunc.R PASS=72 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.32-test-nl-brms.txt:FILE test-nl-brms.R PASS=91 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.32-test-nl-family.txt:FILE test-nl-family.R PASS=321 FAIL=0 ERROR=0 SKIP=0 WARN=0
ob0.3.32-test-nl-stan-funs.txt:FILE test-nl-stan-funs.R PASS=122 FAIL=0 ERROR=0 SKIP=1 WARN=0
ob0.3.32-test-numerical-robustness.txt:FILE test-numerical-robustness.R PASS=723 FAIL=0 ERROR=0 SKIP=0 WARN=0
```
<!-- nlfix:suite_gated end -->

`test-drmtmb-agreement.R`, `test-curve-scale.R` and `test-fuzz.R` have gates of
their own and are the core files that skip under the brms gate. They do not reach
this change beyond what the ungated tier covers.

<!-- nlfix:suite_sess begin -->
```
One session, forward, reverse and preload:
SESSION frmtmb forward files=280/290 pass=20980 fail=0 err=0 skip=217 warn=0
LEAK test-arg-refusal.R compat features: 1 -> 2 entries [with load of frmtmb.sample]
LEAK test-arg-refusal.R compat rules: 0 -> 1 entries [with load of frmtmb.sample]
LEAK test-arg-refusal.R prior defaults: 0 -> 1 entries [with load of frmtmb.sample]
LEAK test-cox-home.R compat features: 2 -> 4 entries [with load of frmtmb.survival]
LEAK test-cox-home.R compat rules: 1 -> 2 entries [with load of frmtmb.survival]
LEAK test-cox-home.R aterm registry: bhaz [with load of frmtmb.survival]
LEAK test-families-extreme.R compat features: 4 -> 6 entries [with load of frmtmb.eam]
LEAK test-families-extreme.R compat rules: 2 -> 4 entries [with load of frmtmb.eam]
LEAK test-families-extreme.R aterm registry: dec ndt_group [with load of frmtmb.eam]
LEAK test-families-extreme.R frame checks: 0 -> 2 entries [with load of frmtmb.eam]
LEAK test-nl-family.R compat features: 6 -> 9 entries [with load of frmtmb.latent]
LEAK test-nl-family.R compat rules: 4 -> 7 entries [with load of frmtmb.latent]
SESSION frmtmb reverse files=280/290 pass=21130 fail=0 err=0 skip=217 warn=0
LEAK test-nl-family.R compat features: 1 -> 4 entries [with load of frmtmb.latent]
LEAK test-nl-family.R compat rules: 0 -> 3 entries [with load of frmtmb.latent]
LEAK test-families-extreme.R compat features: 4 -> 6 entries [with load of frmtmb.eam]
LEAK test-families-extreme.R compat rules: 3 -> 5 entries [with load of frmtmb.eam]
LEAK test-families-extreme.R aterm registry: dec ndt_group [with load of frmtmb.eam]
LEAK test-families-extreme.R frame checks: 0 -> 2 entries [with load of frmtmb.eam]
LEAK test-cox-home.R compat features: 6 -> 8 entries [with load of frmtmb.survival]
LEAK test-cox-home.R compat rules: 5 -> 6 entries [with load of frmtmb.survival]
LEAK test-cox-home.R aterm registry: bhaz [with load of frmtmb.survival]
LEAK test-arg-refusal.R compat features: 8 -> 9 entries [with load of frmtmb.sample]
LEAK test-arg-refusal.R compat rules: 6 -> 7 entries [with load of frmtmb.sample]
LEAK test-arg-refusal.R prior defaults: 0 -> 1 entries [with load of frmtmb.sample]
SESSION frmtmb preload files=280/290 pass=21223 fail=0 err=0 skip=219 warn=0
SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam forward files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
LEAK test-newdata-dec.R compat features: 2 -> 3 entries [with load of frmtmb.sample]
LEAK test-newdata-dec.R compat rules: 2 -> 3 entries [with load of frmtmb.sample]
LEAK test-newdata-dec.R prior defaults: 0 -> 1 entries [with load of frmtmb.sample]
SESSION frmtmb.eam reverse files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
LEAK test-sampling.R compat features: 2 -> 3 entries [with load of frmtmb.sample]
LEAK test-sampling.R compat rules: 2 -> 3 entries [with load of frmtmb.sample]
LEAK test-sampling.R prior defaults: 0 -> 1 entries [with load of frmtmb.sample]
SESSION frmtmb.eam preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.latent forward files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent reverse files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent preload files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn forward files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn reverse files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
LEAK test-records.R envvars: rxTempDir [with load of backports cachem checkmate data.table dparser dplyr farver fastmap ggplot2 gtable lotri memoise pillar pkgconfig Prec
SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
LEAK test-records.R envvars: rxTempDir [with load of backports cachem checkmate codetools data.table dparser dplyr farver fastmap ggplot2 glue gtable lotri MASS memoise p
SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
LEAK test-records.R envvars: rxTempDir [with load of backports cachem checkmate data.table dparser dplyr farver fastmap ggplot2 gtable lotri memoise pillar pkgconfig Prec
SESSION frmtmb.sample forward files=74/79 pass=3782 fail=0 err=0 skip=15 warn=0
LEAK test-compat-preflight.R compat features: 1 -> 7 entries [with load of frmtmb.eam frmtmb.latent frmtmb.ode]
LEAK test-compat-preflight.R compat rules: 1 -> 7 entries [with load of frmtmb.eam frmtmb.latent frmtmb.ode]
LEAK test-compat-preflight.R aterm registry: 0 -> 2 entries [with load of frmtmb.eam frmtmb.latent frmtmb.ode]
LEAK test-compat-preflight.R frame checks: 0 -> 3 entries [with load of frmtmb.eam frmtmb.latent frmtmb.ode]
LEAK test-cox-draws.R compat features: 7 -> 9 entries [with load of frmtmb.survival]
LEAK test-cox-draws.R compat rules: 7 -> 8 entries [with load of frmtmb.survival]
LEAK test-cox-draws.R aterm registry: bhaz [with load of frmtmb.survival]
SESSION frmtmb.sample reverse files=74/79 pass=3782 fail=0 err=0 skip=15 warn=0
LEAK test-simulators.R compat features: 1 -> 6 entries [with load of frmtmb.latent frmtmb.survival]
LEAK test-simulators.R compat rules: 1 -> 5 entries [with load of frmtmb.latent frmtmb.survival]
LEAK test-simulators.R aterm registry: 0 -> 1 entries [with load of frmtmb.latent frmtmb.survival]
LEAK test-compat-preflight.R compat features: 6 -> 9 entries [with load of frmtmb.eam frmtmb.ode]
LEAK test-compat-preflight.R compat rules: 5 -> 8 entries [with load of frmtmb.eam frmtmb.ode]
LEAK test-compat-preflight.R aterm registry: dec ndt_group [with load of frmtmb.eam frmtmb.ode]
LEAK test-compat-preflight.R frame checks: 0 -> 3 entries [with load of frmtmb.eam frmtmb.ode]
SESSION frmtmb.sample preload files=74/79 pass=3782 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.spatial forward files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
LEAK test-spde.R envvars: GDAL_DATA [with load of class classInt DBI e1071 fmesher KernSmooth proxy sf units]
SESSION frmtmb.spatial reverse files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
LEAK test-spde.R envvars: GDAL_DATA [with load of class classInt DBI e1071 fmesher glue KernSmooth pillar proxy rbibutils Rdpack reformulas sf units vctrs waldo]
SESSION frmtmb.spatial preload files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
LEAK test-spde.R envvars: GDAL_DATA [with load of class classInt DBI e1071 fmesher KernSmooth proxy sf units]
NLFIXSESS3 ran 27 of 27
```
<!-- nlfix:suite_sess end -->

<!-- nlfix:suite_check begin -->
```
frmtmb: Status: 1 NOTE
* checking HTML version of manual ... [25s] NOTE
frmtmb.sample: Status: 1 NOTE
* checking HTML version of manual ... NOTE
frmtmb.eam: Status: 1 NOTE
* checking HTML version of manual ... NOTE
frmtmb.learn: Status: OK
frmtmb.survival: Status: 1 NOTE
* checking HTML version of manual ... NOTE
```
<!-- nlfix:suite_check end -->

## 7. Not done, and filed

- Round 2, filed. `a * x` with x about 1e-20 stays at the start's scale: -299.33 against
  -272.50, the same on rellib-r16 (`dev/review-nlfix-c9.R`, seed 6201).
- Round 2, met on the way. `simulate()` on `beta_binomial()` stops with "shape1 and shape2
  must be positive" if a fitted mean rounds to exactly 0 or 1, since `sim` forms
  `1 - mu` where the lpdf uses the log-odds complement; the same on rellib-r16
  (`dev/nlfix-lgscan-*.log`). `frm_allfit(start_from_mle = FALSE)` reports bobyqa "failed:
  0 < ctrl$rhoend is not TRUE" on a body's fit from zeros (`dev/nlfix-bern2.log`; not compared
  with rellib-r16). `Beta()`, `dirichlet()` and the inflated Beta forms keep raw `lgamma()`
  differences: a plain fit does not reach the large phi where they round, and a latent residual
  has latres's cap.
- Reruns of the review's scripts: `dev/review-nlfix-c1.R` section H and `dev/review-nlfix-c3.R`'s
  decomposition call `nl_restart_credible()` and the old `nl_saddle_starts()`, so they stop
  there on this build. `dev/nlfix-r2.R` and `dev/nlfix-cost3.R` measure the same.

- `dev/brmsport-verdicts.tsv` is not edited. `brmsfit-methods:955` flips: the row
  was "defect" on `update(fit2, bf(count ~ a + b, nl = TRUE))` stopping at its zero
  start on Gamma("identity"), and it now refits (`dev/nlfix-u955.R`: logLik -116.1874,
  code 0, with the flat-direction warning `a + b` deserves). The ported suite's other
  rows did not change a result in the gated run.
- `dev/equiv-data/adjudication.tsv` rows p049, p158 and p166, and `sample.tsv`'s
  frozen statuses, are for the consolidation (section 4).
- brms 2.23.0 defect (above): `predictor.bprepnl` drops `predictor_ac()`.
- The `hurdle_cumulative` + disc failure is brms's (`brms-1`). Nothing to fix here.
- `stan_pars_from_fit()` has no rule for an ordinal mixture's or an mv ordinal's
  thresholds under a body (`ord_mixture`, `ord_mv`, L4 only), nor for the simplex
  `theta` outside the instrument's own rule, which is why test-nl-brms.R writes
  `theta1 ~ 1`.
- `nl x` the frmtmb.learn families and gddm, lba, rdm, wiener_gng: not measured.
  The guard stands behind them.
- `conditional_effects()` on a mean or category display of a body keeps its loud
  refusal naming `band = "boot"` and `method = "predict"`.
- `nlf(muB ~ ...)` with a flag and other locations: brms copies the main right-hand
  side ("1") into the others as a body (`muc[n] = (1)`), while frm() gives them an
  intercept. An edge case neither documents; not changed.
- Two brms spellings still refused, both older (review, minor 12). `bf(y ~ 1, nlf(sigma
  ~ ...), a + b ~ 1)`, brms's own: "Cannot interpret bf() argument 1", the same on
  rellib-r16, and it is what a per-category `nlf()` inside `bf()` hits; brms builds it.
  `bf(y ~ a + b * x, mu2 ~ x, a + b ~ 1, nl = TRUE)` under a mixture: brms builds it, and
  frm() says mu2 "cannot also be a nonlinear parameter". The `+ nlf()` route works and is
  an identity (`dev/review-nlfix-specs.R`, L1 per location -2.8e-14 to -1.0e-12).
- A linear mixture's prior table has no `theta` row (brms 9 rows, frm() 8), on a body
  too. Older (review, NIT).
- `beta_binomial()` with no overdispersion: filed in round 1, fixed in round 2 (section 0a).
- `log(a) * x` under gaussian returns logLik -Inf, conv = 0, with a warning, as on
  rellib-r16. Older (review, NIT). NEWS now cites the case measured, under Gamma().
- Filed by the review, not this lane's: a large response from the default start,
  `a + b * x` with y about 1e6, -3045.62 against -1658.76 linear, conv = 0, the same on
  rellib-r16; `frm_curve()` on a linear `cov = FALSE` stops with "assembled covariance
  ... disagrees" on rellib-r16 too, and on a body it returns the curve without the
  shift; student on near-normal data, L1 const 2.6e-4 on a body and 1.3e-4 on the linear
  control, 8.5e-14 on real t data.
- `bf(nl) + acformula(a) + acformula(b)`: the last wins, as in brms (the review verified,
  cosy kept), where a linear formula refuses two. Kept, and the help says so.
