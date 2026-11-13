# Lane boundfam: fits with no maximum, and boundaries of family parameters

Round of 2026-11-11, base e8abed0a (frmtmb 0.81.0), branch
`wt-r82-boundfam`. Brief: `round-briefs/r082-boundfam.md`. This file is
the record; logs named below are local (`dev/boundfam-out/`,
`dev/boundfam-suite-*/`, gitignored).

## Items and the rule that fixes each

| item | what | rule | status |
|---|---|---|---|
| 3333 | no maximum: all-zero poisson, all right-censored gaussian | A | fixed |
| 3098 | zero-width correlation interval | B | fixed |
| 2591 = 2769 | eam mixing weight at 6e-11 called "flat" | A | fixed |
| 2710 | mixture precision with an infinite MLE, silent | A | fixed |
| 3482 | dispersion at its limit with `(1 \| obs)`, sum(w) = n | A | fixed |
| 3022 | frechet nu = 1: no message naming the edge | not A | not done |
| 3037 | GEV -1 < xi < -0.5: no message | not A | not done |

- Rule A, `se_edge_limit()` in `R/se-check.R` (one function, as the
  brief asks, so that lane quadseam can call it): a coefficient of a
  linear predictor of any distributional parameter that has run to the
  end of its range, and the parameters this leaves undetermined.
- Rule B, `varcorr_trans_rows()` in `R/confint.R`: no zero standard
  error, and the boundary of a wider `us()` block's correlation space
  (`varcorr_us_row_edge()`).

## Rule A

### Before: rellib-r20 (`Rscript dev/boundfam-repro.R base <cases>`)

Logs `dev/boundfam-out/repro-base-{1,2,3}.txt`. Estimate, SE, lost
reason; W = warning, M = message.

- pois0 (seed 2, n = 50, y = 0, `y ~ x`, poisson): code 1 (nlminb 9),
  W "Optimizer did not report convergence ... (9). A Newton step from
  the stop gains 6.84e-276 ... so the estimates are at the maximum".
  (Intercept) -636.807 SE 2.703e+137; x -0.0552112 SE 2.414e+137.
  `glm()`: no warning, (Intercept) -26.30 SE 44231.32, x -4.2e-15 SE
  39475.91.
- cens (`y | cens(cen) ~ x`, every row "right", n = 60), seeds 1 to 6:
  - seed 1: code 0, NO warning; SEs 16080, 24610, 1979.
  - seed 2: code 9 W + "flat" W on (Intercept), x; sigma SE 5.923e+126.
  - seed 3: code 0, NO warning; SEs 6.351e8, 1.303e8, 2.112e8.
  - seed 4: as seed 2 (sigma SE 4.809e+128). Seed 5: "flat" W on
    (Intercept), x; sigma SE 2.654e9. Seed 6: code 0, NO warning; SEs
    9.867e7, 4.466e7, 4.746e7.
  - `survreg(Surv(y, 0) ~ x, dist = "gaussian")`: seed 2 no warning,
    SEs 348.3, 1281.6, 10235.0; seed 3 "Ran out of iterations and did
    not converge", Log(scale) -Inf, SEs 0.
- maling (test-bcm-latent-mixtures.R's Malingering_2): code 0, no
  condition; phi2_(Intercept) 23.3575 SE 1.255e+04; the others 0.2348,
  0.5028, 0.5170, 0.4286.
- eam (frmtmb.eam test-defects.R:69's mixture): W "flat" on mu1, ndt1,
  theta1 (theta1 at -23.4422); bs1_(Intercept) kept SE 5.575e+03.
- negbinomial on poisson data (`rpois(200, exp(1 + 0.4 x))`): seed 2
  log shape 15.872 SE 1218, seed 3 17.542 SE 4771, no condition;
  `MASS::glm.nb()` warns "iteration limit reached" twice (theta 9237 SE
  262989; 16488 SE 360085). Seed 1 has a maximum: log shape 4.037 SE
  1.768 (glm.nb theta 56.66 SE 100.1), no condition.
- wobs (3482; `y | weights(w) ~ x + (1 | obs)`, w / mean(w), n = 60):
  gaussian seeds 1 and 3 log sigma -10.06 and -9.99 "flat" (the olre
  warning explains it, so no SE warning); Gamma seed 1 code 8, log
  shape 15.58 "nocurv" under the code-8 warning; Gamma seed 3 code 0,
  log shape 15.58 "flat" W.

### The rule

`se_edge_limit(fit, H, skip)`, called by `cov_from_hessian()` beside
`se_edge_sd()` and `se_edge_cor()` (not under quadrature or importance
sampling, which give no boundary verdict either), its result joining
`edge`, so tier 3 removes those rows:

1. Candidates: the outer positions of every linear predictor's
   coefficients (`se_lp_outer_pos()`), not a nonlinear parameter's (left
   to `nl_flat_message()`), not held by a bound or an earlier edge rule,
   with curvature `|H_jj| u_j^2 <= se_edge_screen` (0.1, as
   `se_edge_sd()`).
2. `se_limit_dir()`: with the rest held, the objective changes by at
   most `se_edge_tol` (1e-6) at `se_edge_step` (2) natural units one way,
   or falls; does not rise by more than 1e-6 from there to 12 units;
   and rises by more than 1e-6 the other way (at 2 units, or, when both
   sides are level at 2, at zero or 10 units back, whichever is
   further; then the end is the side away from zero). The gain along
   the ray is kept: above `grad_tol` the coefficient is "rising".
3. Undetermined: with every end coefficient moved 2 units further, a
   parameter whose curvature falls below `se_limit_decay` (0.5) of its
   value at the stop (AD Hessian diagonal without random effects, else
   central differences of the gradient).
4. Reason "limit" (`se_relabel()`, before the boundary relabel; a
   proved separation still wins). Report (`se_report()`):
   - every limit parameter a coefficient of a family (non-primary)
     dpar, none rising: the boundary message, class
     `frmtmb_boundary_fit` (`se_limit_message()`), as a variance at zero;
   - otherwise the "Standard errors are not available" warning with the
     clauses `limit_end` ("the maximum likelihood estimate does not
     exist"), `limit_rise` ("the fit stopped on the way up") and
     `limit_undefined`.
   - `se_boundary_stop()` accepts "limit" as it accepts "boundary", so
     a nonzero code at such a stop gives this report in place of the
     optimizer's warning, whose text claimed the maximum.
5. `latres_shift_tk()` (R/autocor-latent.R) reads the propagation
   covariance when every lost parameter is "limit": a lognormal sigma
   at exp(-4.43) under `ar(cov = TRUE)` (test-latres.R) otherwise
   turned the latres check's "no curvature" warning on (shift Inf; with
   the propagation covariance, 5.97e-16, `dev/boundfam-latres.R`).

Probe of the internals (`dev/boundfam-probe.R <lib> <case>`, logs
`dev/boundfam-out/probe-*.txt`, lane build before the "rise" change):

- maling: phi2 objective change at -12, -2, 2, 12 units: 0.001033,
  4.055e-08, -5.485e-09, -6.343e-09; back to zero 11.15; dir +1. The
  four others' curvature ratio (q over p) is 1 to eight digits.
- eam: theta1 at -23.44: -7.658e-09, -6.622e-09, 4.893e-08, 0.001246;
  back to zero 118.4; dir -1. mu1, bs1, ndt1 ratio 0.1353353 (exp(-2));
  mu2 and sigma2 1.0000000.
- wobs Gamma seed 3: shape at 15.58: 0.2245, 5.254e-06, -9.426e-07,
  -0.04665. The ray gains 0.047 over 12 units with a gradient of
  5.5e-6 at the stop: a stop on the way up, so "rising" (Gamma seed 1:
  0.195).

### After: the lane build

Logs `dev/boundfam-out/repro-final-{1,2,3}.txt` (final build).

- pois0: one W, "Standard errors are not available for 2 of 2
  parameters. (Intercept) (at -636.8, toward -Inf): ... the maximum
  likelihood estimate does not exist ... x: the data do not determine
  it at that limit ..."; no optimizer warning (`se_boundary_stop()`);
  both SEs NaN, reason "limit".
- cens seeds 1 to 6: each one W naming (Intercept) toward +Inf and
  sigma_(Intercept) toward -Inf at the end of their ranges, x
  undetermined; all three SEs NaN. Seeds 2 and 4 (code 9) lose the
  optimizer's warning, whose text said the stop was at the maximum.
- maling: M (class `frmtmb_boundary_fit`) "Boundary fit:
  phi2_(Intercept) is at the end of its range (23.36 on its link
  scale, moving toward +Inf)"; phi2 NaN; the other four SEs as before
  (0.2348, 0.5028, 0.5170, 0.4286: phi2 is decoupled there).
- eam: one W, theta1_(Intercept) "(at -23.44, toward -Inf)" at the end
  of its range, and mu1, bs1, ndt1 undetermined there; bs1 loses its
  5575.
- negbinomial on poisson data seeds 2 and 3: M "shape_(Intercept) is
  at the end of its range (15.87 / 17.54 ...)"; shape NaN; the
  coefficients keep their SEs. test-se-limit.R checks that these equal
  `glm(family = poisson)`'s to a ratio within 1e-3.
- wobs gaussian seeds 1 and 3: sigma "limit" (NaN, as before), still
  explained by the olre warning, so no second condition. Gamma seeds 1
  and 3: W "shape_(Intercept) (at 15.58, toward +Inf, where the
  log-likelihood rises by a further 0.2 / 0.047): the likelihood still
  rises ..., so the fit stopped on the way up"; seed 1 (code 8) loses
  the optimizer's warning.
- latres (test-latres.R's lognormal `ar(cov = TRUE)`,
  `dev/boundfam-latres.R`): M on sigma at -4.432 toward -Inf. Kept
  SEs: (Intercept) and x unchanged (0.1533, 0.1488), theta_1 0.0638 ->
  0.0500 and theta_2 0.0759 -> 0.0757: the base inverse carried the
  flat sigma direction, the lane's is the limit model's (sigma held), as
  a variance at zero is handled.
- learn (frmtmb.learn test-engine.R's one-subject toy, 25 trials,
  `dev/boundfam-learn.R`, logs `learn-{base,lane}.txt`): base SEs
  24583 and 14643, silent. The profile of log tau (alpha refitted)
  reads -logLik 1.4686097 at 1.97, 1.3864144 at 2.97 and 1.3862943611
  (= log 4) at 3.97, 5.97, 8.97, 13.97 and 23.97: no maximum. Lane: W,
  tau toward +Inf, alpha undetermined, both NaN. The test now allows
  that warning.

### Guard-absent cases

Same logs. Each prints the same estimate and SE table on rellib-r20
and the lane, with no condition on either: pois1 (one count at the
median x; SEs 1.0030 and 0.8928, glm's 1.00255 and 0.89278),
negbinomial on poisson data seed 1 (log shape 4.037, SE 1.768), a
negbinomial with shape 30 (log shape 3.176, SE 0.5708), a
two-gaussian mixture with weight 0.08 (theta1 2.956, SE 0.2439), a
lognormal with sigma 0.02, a rare bernoulli (intercept -3.97), and the
GEV at xi0 = -0.3. A beta_binomial fitted to binomial data (log phi
19.94) is NOT an absent case: it is the binomial limit, and gets the
boundary message.

Units. The first version stepped in `par_units`, or 1 without them.
Under `autoscale = FALSE` a poisson `y ~ 0 + xs` with xs at scale 1e-6
(test-predfix.R) stopped with xs at 686670; a raw step of 2 moved the
predictor by 2e-6 and read level both ways, and zero read a loss, so
xs was called at its end: a false alarm, two escaped warnings in that
file. `se_limit_units()` now gives such a fit natural units (the
predictor's unit over the column's spread), and the false alarm is
gone (`dev/boundfam-suite2-lane-t2/frmtmb--test-predfix.txt`: pass=92
warn=0). Its first version also broke the SE computation of a
frmtmb.learn fit ("subscript out of bounds" on a missing predictor
unit); fixed, and `cov_from_hessian()` now keeps the other rules'
verdicts if the limit rule ever errors.

## Rule B: correlation intervals

- A zero se on any `varcorr_trans_rows()` row is NA: the delta method
  gives zero only where the jacobian row is flat or the quadratic form
  of an indefinite inverse rounds to zero or below (`pmax(..., 0)`).
  Shown by mocking `sdr_of()` to zero the theta block (test-se-limit.R,
  "a zero standard error gives no interval"): rellib-r20 fails 3
  expectations there (se_t 0, interval of zero width, no warning).
- A `us()` block of dimension 3 or more: the correlations of a row
  whose multiple correlation with the earlier rows, `sqrt(S / (1 + S))`
  with `S` the squared length of the unit-diagonal Cholesky row's
  off-diagonal part, is at least 0.9999 (the clamp; in dimension 2 the
  existing test exactly) get NA.
- The review's seed-7 script (`rev072-seed7.R`) is not in the tree or
  in history. `dev/boundfam-us3.R` builds a rank-2 three-term block
  (`u2 = 0.4 u0 - 0.6 u1`, 30 groups of 10). On rellib-r20
  (`dev/boundfam-out/us3-base-1.txt`) seeds 1 to 12: 10 stop with the
  third row past the clamp and lose theta_5 and theta_6 as "nocurv"
  (NA); seed 6 stops at code 1 with the row at (1770.5, -2270.7),
  nothing lost, and cor(z,(Intercept)) se_t 0.0009159 (a Fisher-z
  interval 0.0036 wide). Lane (`us3-lane-1.txt`): seed 6 NA for both
  row-3 correlations; seed 10 (code 0, row (4.66, -3.85)) keeps every
  interval, se_t 0.217 and 0.2096. No seed gave se_t exactly 0 on this
  design; the zero guard is shown by the mock above.

## SEs on the test suites, before and after

Every standard-error analysis (`cov_from_hessian()`) and covariance
interval table (`varcorr_trans_rows()`) of every test file was traced,
one file per process (`dev/boundfam-run-traced2.R` through
`dev/boundfam-par2.sh`; base core by `dev/boundfam-run-traced.R`, the
same record without the LIM detail). Base: rellib-r20, all 317 core
files (`dev/boundfam-suite-base/`) and the 21 extension files where the
rule fired (`dev/boundfam-suite2-base/`). Lane: the final build, all 510
files (`dev/boundfam-suite2-lane/`).

`Rscript dev/boundfam-cmp.R dev/boundfam-suite-base
dev/boundfam-suite2-lane` (`dev/boundfam-out/cmp-final.txt`), core:

    files compared 317
    COV bitwise=5353 ulps=3 moved=1 nan=54 length=0 count=0
    VCT bitwise=144 ulps=0 moved=0 nan=0 length=0 count=0

- "moved" 1: test-sugar.R record 2, `y ~ x` on `rnorm(50)` with no seed
  (different data each run).
- "nan" 54: each a record where the limit rule fired (below), in 25
  files. An earlier lane run also showed one in test-id-kron.R with no
  firing: three base runs of that file in fresh processes gave theta_2
  NaN in two and SE 0.9613 in one (hybrid-CPU rounding,
  dev/rtmb-pitfalls.md item 21), and the lane gave both too.

Every firing (`Rscript dev/boundfam-classify.R dev/boundfam-suite2-lane
dev/boundfam-suite-baseall`, `dev/boundfam-out/classify-final.txt`):

    TOTAL fired=89 already_nan=16 newly_nan=73 no_base=0

- 89 firings over core and extension suites; 16 only relabel parameters
  that were already NaN; 73 take finite base SEs to NaN. Those base SEs
  (86 values) range from 0.0494 to 6.4e12, median 6485.
- Every firing's estimate is far out on its link (|est| from 8.68 to
  4.96e4: student nu, beta_binomial and dirichlet_multinomial phi,
  negbinomial shape, xbeta kappa, wiener lambda and ndt, zero-inflation
  and mixing logits, gaussian sigma at exp(-8.7) to exp(-23.1) on
  perfect or obs-level fits, an HMM transition at -3043, a categorical
  category with no rows), except four, each checked: frmtmb.learn's
  one-subject tau at 3.97 (profile level from 3.97 up, above); its
  kalman sigma0 at -5.45 (rises 2.6 within 12 units: a stop short of
  the maximum, which 0.81.0 called "concave"); an HMM transition slope
  at -5.55 and the lognormal sigma at -4.43 under `ar(cov = TRUE)`
  (above). The smallest base SEs taken: 0.0494 on frmtmb.coupling's
  coherence fit whose logit ran to 18.5 with 647 more log-likelihood
  12 units out (a stop short; there tier 3 then also lost both powers,
  which it does in the next record on rellib-r20 too), and 4.42 on a
  gp() fit whose sigma ran to exp(-18) with 0.0128 still to gain.
- Kept SEs: of the 83 records with parameters kept, 74 keep them within
  0.1 percent; 9 move, all on fits with no maximum: test-ce-bands.R's
  perfect fit (coefficient SEs 1e-11, ratio 0.47 to 0.99), the latres
  lognormal (theta_1 ratio 0.784), test-ordinal-mixture.R's components
  run off in x (ratios down to 0.0127 on parameters whose base SE was
  10 to 8e5), the coupling coherence (0.998), frmtmb.eam's lba `A` at
  exp(-10.1) and exp(-10.7) (0.906, 0.955). What moved is the
  covariance of the other parameters with the limit one held, the limit
  model's, as for a variance at zero.

## Tests

- New `tests/testthat/test-se-limit.R` (9 tests in the first round).
  On rellib-r20: pass=41 fail=132 err=3. The log this line first cited
  was overwritten by a later run; the reviewer's rerun
  (`dev/reviews/boundfam-review-selimit-base.txt`) reads the same
  counts. Only one test errors on the missing `varcorr_us_row_edge()`;
  the other two errors (negbinomial, Malingering_2) follow a behavioral
  failure (no message), and the three-term fit test computes the row
  test inline and fails behaviorally. Lane: pass=186 fail=0. The fix
  round's counts are under "Fix round".
- frmtmb.eam `test-defects.R`: the mixture fit now requires the
  warning to name `theta1_(Intercept) (at `, `toward -Inf` and
  `bs1_(Intercept)`; rellib-r20 fails those 3.
- Changed to the new verdicts: `test-diagnostics-ux.R` (student `nu`
  run off: SE NaN, was "greater than the estimate"),
  `test-gp-start.R` (sigma e^-10 under a gp(): the flat warning or the
  boundary message), `test-cens-trunc.R` (the trunc window fixture:
  rows spread evenly over [-1, 2] put the truncated normal's maximum at
  sigma = Inf, the uniform; now rows concentrated in the window).

## Not done, with the reason

- A joint runaway (filed in the backlog). pois1max of
  `dev/boundfam-repro.R` (poisson, seed 2, n = 50, the one count at the
  largest x): rellib-r20 and the lane alike give a code-7 warning that
  says the stop is at the maximum, and the "flat" warning on
  (Intercept) -436.5 and x 208.8 (`repro-base-1.txt`,
  `repro-final-1.txt`); glm() warns "algorithm did not converge" and
  "fitted rates numerically 0 occurred". Each coefficient alone costs
  likelihood at that stop, so a one-coordinate probe cannot read the
  runaway. The verdict also depends on the coding (review, m2):
  `y ~ g` with an all-zero group is "flat", `y ~ 0 + g` is "limit".
- 3022 (frechet near nu = 1). Not the boundary rule: the rule does not
  mark nu_(Intercept) at -19.14, and rellib-r20 and the lane give the
  same "flat" warning on (Intercept) (17.38) and nu_(Intercept)
  (`repro-base-2.txt`, `repro-final-2.txt`). By brms's
  parameterization the mean intercept and nu run together (the mean is
  infinite at nu = 1), a joint ridge that a one-coordinate probe does
  not read; not measured further. The fix the item asks for is a
  family fit-end check like `gev_fit_check()`.
- 3037 (GEV xi in (-1, -0.5)). Not a boundary: the maximum is interior
  and non-regular, so the fix is a message in `gev_fit_check()`. Seen
  on rellib-r20: xi0 -0.63, -0.7, -0.89 (seed 4, n = 300) fit at xi
  exp(-1.087) - 1 = -0.663, -0.736 and -0.936 with code 0 and no
  condition (`repro-base-2.txt`).

## Cost

`dev/boundfam-cost.R` and `dev/boundfam-cost2.R` (logs `cost.txt`,
`cost2.txt`): `se_edge_limit()` traced with its time and its objective
and gradient calls, three replicates each, n = 2e4.

- No coefficient screened (a gaussian `y ~ x + (1 | g)`, negbinomial
  fits whose shape has a maximum, with `(1 | g)`): 0 evaluations, 0.00
  to 0.02 s of 1.0 to 17.8 s.
- Screened, not fired (negbinomial `y ~ x`, log shape 7.03 at a
  maximum): 3 objective evaluations, 0.02 to 0.03 s of 1.35 to 1.77 s.
- Fired with `(1 | g)` (negbinomial on poisson data, log shape 15.39):
  4 objective evaluations and 12 gradients (the curvature test, four
  per other outer parameter), 3.6 to 4.1 s of 23.6 to 25.3 s. The
  gradients are the cost; a fit pays it only where the rule fires.

## OpenBLAS

`dev/ciharden-openblas.sh 0.3.26 lapack` and `0.3.32 lapack`, test
files through `dev/boundfam-ob.sh` (logs `dev/boundfam-out/ob-*/`,
`lib:` lines show the lane's frmtmb). 0.3.26: test-se-limit (186),
test-se-tier (112), test-latres (305), test-cens-trunc (99),
test-diagnostics-ux (130), test-gp-start (238), test-review-v29 (186),
test-predfix (92), frmtmb.eam test-defects (63), frmtmb.learn
test-engine (20): every one fail=0 err=0 warn=0. 0.3.32: test-se-limit,
test-se-tier, test-latres, test-diagnostics-ux, test-review-v29,
frmtmb.eam test-defects: the same counts, fail=0 err=0 warn=0.
`dev/ciharden-out/` deleted afterwards.

Process count: the first 0.3.26 launch used detached subshells, which
the job limiter could not see, and 7 R processes ran for a few seconds
before I stopped 4 of them (my own, matched on the OpenBLAS build's
path in my worktree) and reran them through `dev/boundfam-ob.sh`. The
build script's own R check also ran beside the lane's 3 suite
processes for about 15 s.

## Files touched

Package code: `R/se-check.R` (se_edge_limit(), se_limit_dir(),
se_limit_units(), se_limit_info(), se_limit_candidates(),
se_limit_family_par(), se_limit_curv(), se_limit_decay,
se_limit_quiet(), se_limit_message(), se_limit_where(); hunks in
cov_from_hessian(), se_report(), se_relabel(), se_lost_clauses(),
se_boundary_stop()), `R/confint.R` (varcorr_trans_rows(),
varcorr_us_row_edge(), confint_varcorr() docs and warning text),
`R/autocor-latent.R` (latres_shift_tk(), lane arspat's file: one hunk),
`man/confint_varcorr.Rd`, `NEWS.md`.

Tests: new `tests/testthat/test-se-limit.R`; changed
`tests/testthat/test-diagnostics-ux.R`, `test-gp-start.R`,
`test-cens-trunc.R`, `extensions/frmtmb.eam/tests/testthat/test-defects.R`,
`extensions/frmtmb.learn/tests/testthat/test-engine.R`.

Records: `dev/test-backlog.md` (new section "Filed in the round of
2026-11-11", CLOSED lines under 2591, 2710, 2769, 3098, 3333, 3482),
this file, and the `dev/boundfam-*` scripts.

## Fix round (review `dev/reviews/2026-11-11-boundfam.md`: FIX)

Each finding was reproduced on the first round's build before its fix
(logs `dev/boundfam-out/fix/*-pre.txt`) and rerun after it
(`*-post.txt`), with the reviewer's scripts.

### B1: NaN at a real maximum when a smooth gives `par_units`

Reproduced (`dev/reviews/boundfam-review-units.R lane`,
`fix/units-pre.txt`): started at the autoscaled fit's maximum, with
`autoscale = FALSE`, `y ~ xs + s(z)` gave xs NaN and the "does not
exist" warning for poisson at scale 1e-6 (xs 325886) and bernoulli at
1e-5 and 1e-6 (12904, 129035). Poisson at 1e-5 did not fire.

Fix: `se_edge_limit()` now steps every candidate coefficient in
`se_limit_units()` (the predictor's unit over the column's spread),
whatever `par_units` holds; other parameters keep `par_units`.

After (`fix/units-post.txt`): no condition, and the base SEs 43402.7,
13007.9 and 130078 (poisson 1e-5: 4340.27, as before). The stop short
of the maximum from the default start (the review's "units" case,
`boundfam-review-fp.R lane units`, `fix/fp-units-post.txt`) no longer
fires: objective 502.41121629, xs 0.00123, SE 41370.7 and no condition,
as on rellib-r20, so the false "rises by a further 0.0022" is gone.
That fit stops 29.2 short of the maximum (473.18771786) with code 0 and
nothing says so, on both builds. It is an explicit `autoscale = FALSE`
on a column at scale 1e-6, outside this lane's rule; not filed
separately (test-predfix.R covers the default, which autoscales). The new test "a raw-unit column
beside a smooth keeps its SE at the maximum" passes on rellib-r20 (base
keeps the SEs) and on the lane. The first round's build failed it, as
`fix/units-pre.txt` shows.

### M1: a separating zi or hu slope was a quiet family limit

Reproduced (`boundfam-review-fn.R lane zisep`,
`boundfam-review-split.R lane`; `fix/fn-pre.txt`, `fix/split-pre.txt`):
`bf(y ~ 1, zi ~ x)`, seeds 1 and 2, gave the message only, code 1,
zi_(Intercept) SE 1063.36 and 1562.44. The hu form gave SEs 1436.65
and 2172.91.

Fix: `se_limit_separating()`. Take an end coefficient of a probability
parameter (dpar `zi`, `hu`, `zoi`, `coi` or `theta`, numeric suffix
allowed) whose design column has rows of both signs. It separates
that parameter's rows. It keeps reason "limit", with the clause "the
data separate the rows of its probability ... the maximum likelihood
estimate does not exist". Every other coefficient of that predictor
gets NaN ("sepother"). `se_limit_quiet()` is FALSE, so the fit warns.

After (`fix/fn-post.txt`, `fix/split-post.txt`): zi and hu, seeds 1
and 2, give one W each; both coefficients of the probability are NaN;
the mean intercept keeps its SE (0.0369, 0.0358, 0.0624, 0.0602, as
before). Quiet controls:

- `zi ~ g` with group b poisson has a maximum: SEs 0.186 and 1.13, no
  condition.
- `shape ~ g` seed 2: shape_gb at 16.64 gets the boundary message, and
  shape_(Intercept) keeps SE 0.150.
- A zi intercept at its limit (`dev/boundfam-ziseed.R`): seeds 2 and 6
  of 12 reach -19.6 and -21.0 and get the boundary message; the other
  ten have a maximum at -3.4 to -5.7 and no condition.

New tests: "a zero-inflation or hurdle slope that separates warns" and
"a probability intercept or one group at its limit stays quiet" (seed
2).

### M2: the predicate for lane ppdraw

`se_limit_draw_held(fit, lost)` in R/se-check.R, with `lost` =
`sdr$se_lost`. It returns a logical vector named by `names(lost)`:
TRUE exactly for a parameter with reason "limit" and kind "end" when
`se_limit_quiet(fit, lost)` is TRUE (every limit parameter a
coefficient of a family parameter, none on a rising ray, no
separation). It is FALSE for every "undefined" or "sepother"
parameter, and for every limit parameter of a fit that warns. Wiring
for the merge:

- `boundary_only` in `fit_draw_space()`:
  `all(lost == "boundary" | (lost == "limit" & se_limit_quiet(fit, lost)))`.
- The raw-row restore (`draw_limit_restore()`): only the positions
  where `se_limit_draw_held(fit, lost)` is TRUE.
- A warning-class limit fit keeps the plug-in fallback.

### m1: a limit missed behind a convergence warning

Reproduced (`boundfam-review-fn.R lane nb17`): log shape 23.2154,
code 8, shape "nocurv", with the false-convergence and large-gradient
warnings. A non-finite Hessian diagonal is now admitted
(`se_edge_limit()` asks such a candidate of the objective), but that
is not what misses here:

- On this fit the AD Hessian diagonal is finite (600, 803.7,
  -1.67e-5), and `se_edge_limit()` marks the shape "end" when called
  (`dev/boundfam-nb17.R`).
- The SE analysis is never applied. The gradient of 1.35e-3 exceeds
  `grad_tol`, so the convergence check warns and `se_explained` is
  "convergence". `se_boundary_stop()` requires a gradient within
  `grad_tol`, and `sdr_rescue_hessian()` keeps sdreport()'s inverse at
  a stop short, by design. Both arms behave the same.

Not changed: holding a limit parameter out of the convergence gradient
is a change to `check_convergence()` in R/fit.R, which four other lanes
edit. Filed in the backlog.

### m2: the coding dependence

Reproduced (`fix/fn-pre.txt`, grp0):

- `y ~ g` (group a all zero): "flat" on (Intercept) -22.31 and gb
  23.53.
- `y ~ 0 + g`: ga "limit" at -22.25, gb SE 0.0769 kept.
- zip and negbinomial `y ~ g`: "flat".

Unchanged after the fix (`fix/fn-post.txt`). Added to the filed
joint-runaway item.

### m3: stops short of an interior point

The rising-ray case (gain above `grad_tol`) was already reported as a
stop on the way up and never as "does not exist". The review found two
such stops: frmtmb.coupling's coherence (gain 647, a re-optimized
maximum about 6 units out) and frmtmb.learn's kalman sigma0 (gain
2.6).

Change: at a rising ray the undetermined test is skipped
(`se_edge_limit()`), because the curvature further out says nothing
about the stop. The coherence's two powers had been lost by tier 3
itself, not by the undetermined test (the first round's LIM line named
coh_(Intercept) only), and rellib-r20 loses them in the next record of
that file too.

A re-optimized check at +2 units was not built:

- The gain test already separates every case the review measured. The
  two stops short gain 647 and 2.6. The review's 10 random firings,
  all level or monotone under a re-optimized profile, gain at most
  7.76e-5.
- One refit at +2 units cannot tell a maximum further out from a rise
  without end; that needs the profile itself.

### m4 to m6

- m4: NEWS now says that the rule does not run under quadrature or
  importance sampling.
- m5: the cost where the rule fires (four gradients per other outer
  parameter) is stated in `se_edge_limit()`'s documentation and in
  NEWS. It is not counted in the fit-time budget.
- m6: the boundary message now names both ends (at its lower end
  student() becomes the Cauchy).
- m6, backlog: the section moved to just above "## Filed at the 0.81.0
  release (2026-11-10)", with a blank line before that heading. The
  CLOSED lines are indented two spaces under their items.
- m6, citations: the pois1max citation points at "Not done" above, and
  the test-se-limit.R base citation points at the reviewer's log.

### Runs (final build of the fix round)

- Changed test file on rellib-r20: test-se-limit.R pass=60 fail=141
  err=4 (`dev/boundfam-suite2-base/frmtmb--test-se-limit.txt`). Lane:
  pass=221 fail=0 err=0 warn=0 (`dev/boundfam-suite2-lane/`).
  test-se-limit.R ran on its own each time, because the full-suite job
  lists of both rounds were built before the file existed.
- One file per process, traced (`dev/boundfam-suite3-lane/`; every
  `lib:` line names the lane's core):
  - frmtmb: 317 files, pass=23487 fail=0 err=0 skip=239 warn=0 (plus
    test-se-limit.R above).
  - frmtmb.eam: 33 files, pass=1851 fail=0 err=0 skip=5 warn=0.
  - frmtmb.learn: 16 files, pass=474 fail=0 err=0 skip=13 warn=0.
  - frmtmb.spatial test-sar.R: pass=106 fail=0 err=0 skip=4 warn=0.
- Gated core (`FRMTMB_BRMS_FIT_TESTS=true`), the 28 test-brms-*.R
  files (`dev/boundfam-suite3-gated/`): pass=4737 fail=0 err=0 skip=0
  warn=0.
- SE comparison against rellib-r20, core (`dev/boundfam-out/cmp-fix.txt`):
  `COV bitwise=5353 ulps=3 moved=1 nan=54`, VCT `bitwise=144`. These
  are the same 26 files as the first round's final comparison, line for
  line. I removed that final run's traces by mistake during this round;
  `cmp-final.txt` and `classify-final.txt` keep their results.
- Against the first round's earlier full trace
  (`dev/boundfam-suite2-lane-old/`, `cmp-fix-vs-prev.txt`), the only
  firings that changed are the three removed in test-predfix.R (xs
  twice and alpha_xs, the raw-unit false alarm). The other differences:
  - test-sugar.R's unseeded data;
  - test-id-kron.R's known cross-process flip (dev/rtmb-pitfalls.md
    item 21);
  - one gp-start record that equals rellib-r20.
- `classify-fix.txt` (core, eam, learn, sar): fired=84 already_nan=15
  newly_nan=69 no_base=0. No fit in the suites is a separating
  probability.
- OpenBLAS 0.3.26 with LAPACK (`dev/boundfam-fixrun2.sh`, logs
  `dev/boundfam-out/ob-0.3.26/`, every `lib:` line the lane's core):
  test-se-limit 221, test-se-tier 112, test-latres 305, test-cens-trunc
  99, test-diagnostics-ux 130, test-gp-start 238, test-review-v29 186,
  test-predfix 92, frmtmb.eam test-defects 63, frmtmb.learn test-engine
  20. All of them fail=0 err=0 skip=0 warn=0. `dev/ciharden-out/`
  deleted.
- `R CMD check --as-cran`, core (`dev/boundfam-check.sh .`, log
  `dev/boundfam-check/frmtmb/check.log`): Status: 1 NOTE, the V8 math
  rendering note of the HTML manual, which is environmental. Tests ran
  27 min. The first round's pending checks, as the coordinator read
  them: core 1 NOTE (V8), frmtmb.eam 1 NOTE (V8), frmtmb.learn OK; the
  eam and learn sessions were clean.

### Held limit directions in predictions (ppdraw's short check)

Found by lane ppdraw's reviewer on a merged ppdraw + boundfam build
(`frmtmb-wt-ppdraw/dev/reviews/ppdraw-short-fitted*`). On backlog item
3269 (a zero_inflated_binomial whose zi intercept ran to its limit, a
quiet limit), fitted() gave NaN Est.Error on all 200 rows with the
"predictions move along a direction the fit does not determine"
warning. On 0.81.0 it equals the binomial limit model's to 2.8e-10,
and the rows move by at most 2.7e-9 per logit unit along zi. The null
basis kept zi's unit direction, and jc_nonest() reads its cosine in the
unit-diagonal metric, where zi's tiny curvature makes it large.

Fix: `se_null_drop_held()`, the last step of `sdr_rescue_hessian()`,
drops from `sdr$se_null` the columns that load only on parameters
`se_limit_draw_held()` holds (a family parameter at its end in a quiet
fit). The propagation covariance already has a zero row there, so the
parameter is read at its estimate with no variance, as a variance
component at zero is. Undetermined parameters, separations and the
limits of a warning fit keep their direction, so a quantity that reads
them stays NaN.

Measured (`dev/boundfam-held.R <lib> <dir>`, logs
`dev/boundfam-out/fix/held-{base,lane}.log`):

- 3269 fitted(), rellib-r20: no warning, no NaN, Est.Error 0.173144681621,
  0.245344298847, 0.133920325672, ...; max |zib - binomial| 2.764e-10.
- 3269 fitted(), lane: no warning, no NaN, Est.Error 0.173144681421,
  0.24534429879, 0.133920325599, ... (within 2.0e-10 of rellib-r20);
  max |zib - binomial| 2.263e-10.
- All-zero poisson, fitted(newdata = 3 rows): rellib-r20 Est.Error
  1.05e-139, 7.41e-140, 9.41e-140 in silence; lane NaN on all 3 with
  the warning (x is undetermined).
- 3269 predict(ndraws = 400), seed 11: rellib-r20 Est.Error 3.698,
  2.992, 2.443, ...; lane 1.254, 1.851, 1.406, ... with one warning. On
  this lane alone the draws take the plug-in fallback ("the covariance
  of the estimates could not be recovered"), because `boundary_only`
  in `fit_draw_space()` (R/brms-shapes.R, ppdraw's) does not accept
  "limit": the review's M2, wired at the merge through
  `se_limit_draw_held()`. 0.81.0's 3.7 comes from drawing zi's logit
  with its huge standard error, which puts many draws at zi near 1; the
  merged build holds zi and gives brms's mean (ppdraw's
  `ppdraw-short-merged-cases.txt`: brms mean 5.7022). A test that
  predict() equals 0.81.0 within 1e-8 would pin that defect, so I did
  not write it. The fitted() test is the one asked for.

New test "a family parameter held at its limit keeps fitted() SEs": the
3269 fitted() Est.Error over the binomial model's is 1 within a ratio of
1e-6, with no warning; the all-zero poisson's fitted values are NaN
with the warning. On rellib-r20 it fails 3 expectations: no "limit"
reason, and the poisson's finite SEs with no warning
(`dev/boundfam-suite2-base/frmtmb--test-se-limit.txt`: test-se-limit.R
pass=62 fail=144 err=4).

Runs on this build:

- Core, all 318 files, one per process (`dev/boundfam-suite4-lane/`,
  every `lib:` line the lane's core): pass=23713 fail=0 err=0 skip=239
  warn=0. test-se-limit.R alone: pass=226.
- Gated core, the 28 test-brms-*.R files (`dev/boundfam-suite4-gated/`):
  pass=4737 fail=0 err=0 skip=0 warn=0.
- `R CMD check --as-cran`, core: Status: 1 NOTE (V8), tests 20 min.
  `dev/boundfam-check/` deleted afterwards.

## Consolidation fixes (2026-11-12)

The consolidation (`frmtmb-wt-rel082`, rellib-r21; `dev/round-20261111.md`,
conflicts 2 and 4) found two defects that trace to the limit rule.
Both are reproduced on this lane's build alone and fixed here.

Files changed in this step, for the re-merge:

- `R/se-check.R`: `se_null_drop_held()` (marks limit directions in
  the attribute `limit`), `jc_nonest()` (the move test for those
  directions), `se_report()` and `se_boundary_act()` (the subclass
  `frmtmb_limit_fit`).
- `R/par-template.R` (lane nlstop's file, two hunks): `nl_flat_warning()`
  and `nl_fit_off_saddle()`.
- `R/predict.R` (lane ppdraw's file, one hunk):
  `joint_cov_repair_meta()` carries the attribute to the joint
  positions.
- Tests: `tests/testthat/test-se-limit.R` (two new tests),
  `extensions/frmtmb.eam/tests/testthat/test-defects.R` (residuals now
  finite), `tests/testthat/test-conditions-census.R` (regenerated, one
  exemption).
- `dev/conditions-census-gen.R`: the exemption and its reason.
- `NEWS.md`: two bullets.

### 1. fitted() NaN on a zi or hu slope at its limit

Reproduced with the consolidator's `dev/rel082-fitted-held.R` (copy in
`dev/boundfam-out/cons/`), run from the merged tree on this lane's
library before the fix (its output was later overwritten by the post-fix
run, so the numbers are as printed then). With `zi ~ x` (ppdraw's case,
zi at its end) and `hu ~ x` (no y is 0), fitted() gave NaN Est.Error
with the warning on 195 of 200 and 297 of 300 rows, as the
consolidation reports. The intercept is at
its end and the slope undetermined (both "limit"). `se_null_drop_held()`
held only the intercept, and jc_nonest() judged the slope's direction by
its cosine in the unit-diagonal metric. The slope's curvature vanishes
there, so a move of about 1e-9 per unit read as large.

Fix: a null direction that loads only on limit parameters (an end that
is not a separation, or a parameter it leaves undetermined) is marked
(`attr(sdr$se_null, "limit")`). For such a direction, jc_nonest() asks
how far the quantity moves along it per unit (`|G n|`), against what the
parameters with a standard error move it by (`||G_k / sqrt(H_kk)||`,
about its standard error). It is NaN only above `se_pred_tol` (1e-6)
of that. Separations ("separation", a separating probability's end,
"sepother") and every other lost direction keep the cosine test.

After (`dev/boundfam-cons-fitted.R`, `dev/boundfam-out/cons/consfit-*`):

| fit | rellib-r20 | lane |
|---|---|---|
| zi ~ x, 200 rows | 0 NaN, no warning | 0 NaN, no warning, max rel. diff 1.72e-8 from 0.81.0 |
| hu ~ x, 300 rows | 0 NaN, no warning | 0 NaN, no warning, max rel. diff 2.83e-9 |
| separated bernoulli (kind b) | NaN 239 of 240, warning | NaN 240 of 240, warning |
| all-zero poisson (kind b) | 0 of 50 NaN, SEs 1e-139, silent | NaN 50 of 50, warning |

3269 (`dev/boundfam-held.R`) is unchanged: 0 NaN, within 2.0e-10 of
0.81.0. frmtmb.eam's mixture residuals (weight 6e-11, the wiener
parameters undetermined) are now finite. They move along those
directions by a negligible amount, so test-defects.R asserts finite
values and no longer requires the prediction warning.

ppdraw's `tests/testthat/test-predict-draw-edge.R`, on a scratch build
of this lane's tree with ppdraw's `R/predict-brms.R` and
`R/brms-shapes.R` (`scratchpad/bf-lib`, never a lane library or
rellib-r21; log `dev/boundfam-out/cons/pde-scratch.txt`): pass=45
fail=0 err=0 skip=0 warn=0.

### 2. gp student rows of lane gpstop's sweep

Reproduced (`dev/boundfam-gp.R <lib> <rds>`, logs
`dev/boundfam-out/cons/gp-{base,lane}.txt`):

- `exp_a_student_n40_s2`: 10.232192 on rellib-r20, 9.068636 on the
  lane.
- `pow_a_student_n120_s2`: 55.665287 against -4.946904.

Mechanism: the nonlinear-start code (`nl_fit_off_saddle()`) restarts
from moved coefficients when the default start's held warnings report
a flat direction (`nl_flat_warning()`: "Standard errors are not
available"). On rellib-r20, the default start ran student()'s nu to
exp(22.7) (pow_a), nu was "flat", and the warning sent frm() to restart;
the restart reached 55.665287. On the lane, nu at its end is a quiet
family limit, so the message replaced the warning and no restart ran.
exp_a follows the same path: the default start's 9.439723 on rellib-r20
came from that restart, and the sign multistart then reached 10.232192.
Not se_boundary_stop() and not an escape. The SE check changes no
estimate; the restart decision read its report.

Are 0.81.0's points maxima (`dev/boundfam-gpmax.R`, `gpmax.txt`, a
refit on the lane from 0.81.0's estimates)?

- pow_a: yes. From 55.665287 the refit stays there, code 0, max|grad|
  5.2e-4, no condition, every SE finite (nu's 4.08 on its link). My
  change must not move it.
- exp_a: no. From 0.81.0's point the refit climbs to 10.605195, code 1,
  max|grad| 3.86e11, false convergence; 0.81.0's own stop had code 1
  and a gradient of 1.03e11. That is a stop on a spike, not a maximum.

Fix: the quiet limit message carries the subclass `frmtmb_limit_fit`
(beside `frmtmb_boundary_fit`). `nl_fit_off_saddle()` holds it with
the default start's warnings and releases it only if that fit is
returned, and `nl_flat_warning()` counts it as a flat direction.

After (`gp-lane2.txt`): exp_a 10.232192 and pow_a 55.665287, 0.81.0's
values to the printed digits, with the restart message, and the
replaced fit's limit message is not released. New test "a family limit
on a nonlinear fit still triggers the restart" (pow_a, n = 120).

### Census of bare conditions

The release of the held limit message in `nl_fit_off_saddle()` is a
bare `message(w)`, which `test-conditions-census.R` refuses unless
exempt. The exemption is added in `dev/conditions-census-gen.R`, with
its reason beside the two `warning(w)` it joins, and the eight census
tests regenerated. Only core's copy changes in content (one line); the
seven extension copies are byte-identical to e8abed0a (`cmp`), and
git shows them as modified only by their file time.

### Runs (`dev/boundfam-fixrun4.sh`)

- test-se-limit.R on rellib-r20 (`dev/boundfam-suite5-base/`): pass=71
  fail=146 err=4. Lane: pass=237 fail=0 err=0 warn=0.
- One file per process, traced (`dev/boundfam-suite5-lane/`, every
  `lib:` line the lane's core):
  - frmtmb: 318 files, pass=23723 fail=1 err=0 skip=239 warn=0. The
    one failure was test-conditions-census.R (the bare `message(w)`
    above). After the exemption it reads pass=11 fail=0.
  - frmtmb.eam: 33 files, pass=1851 fail=0 err=0 skip=5 warn=0
    (test-defects.R with the residuals now finite).
  - frmtmb.learn: 16 files, pass=474 fail=0 err=0 skip=13 warn=0.
  - test-gp-start.R: pass=238 fail=0.
  - The core run covers every file that tests fitted() or predict()
    standard errors.
- Gated core, the 28 test-brms-*.R files (`dev/boundfam-suite5-gated/`):
  pass=4737 fail=0 err=0 skip=0 warn=0.
- `R CMD check --as-cran`, core, on the lane library: Status: 1 NOTE
  (V8), tests 15 min. That check ran before the census exemption, and
  the census test skips without the source tree. `dev/boundfam-check/`
  deleted.
- ppdraw's test-predict-draw-edge.R on the scratch build: pass=45, 0
  warnings (above).

### 3. No warning on NaN fitted values under OpenBLAS

The consolidation's OpenBLAS 0.3.26 run failed test-se-limit.R:200 (the
separated bernoulli, seed 514; its probes are
`frmtmb-wt-rel082/dev/rel082-ob-sep*.R`). Under the reference BLAS the
fit stops at code 0, the SE analysis builds a null basis, and fitted()
gives 240 of 240 NaN with the "move along a direction the fit does not
determine" warning. Under OpenBLAS it stops at code 1, the separation
warning is the convergence verdict (`se_explained`), sdreport()'s
inverse is kept with NaN rows and no null basis, and fitted() gives 240
NaN with no warning.

Reproduced here (`dev/boundfam-sepwarn.R`): under OpenBLAS 0.3.26, a
scratch build of this tree with the fix below removed gave "default:
code 1, se_explained convergence, fitted NaN 240 of 240, warnings 0".

Fix: `jc_nonest()` (R/se-check.R) is now `jc_nonest_null()`, the old
test over the null basis, OR'd with a new `jc_nan_lost()`. The new test
marks a row that reads a lost parameter (`jc$lost_pos`) whose variance
in `jc$V` is not finite. So every caller that warns through
`se_pred_warn(jc_nonest(...))` (fitted(), predict(), frm_linpred(),
hypothesis() on the joint covariance, conditional_effects()) gives the
same words on either route. A lost parameter held at zero variance (the
propagation covariance's zero row, a held family limit) is not marked.

After, the same probe on the lane library:

| route | result |
|---|---|
| OpenBLAS 0.3.26, default (code 1) | fitted NaN 240 of 240, 1 warning; frm_linpred() the same |
| OpenBLAS 0.3.26, stopped short | the same |
| reference BLAS | unchanged, NaN 240 of 240, with the warning |

New unit test "a quantity reading a NaN variance is non-estimable on any
route" (a constructed `jc`, with the absent cases: zero variance, no
lost position). The test's assertion of the warning stays.

The fix gives one more fit a warning that had been silent: in
test-mi-families.R, asym_laplace's kink parameter has no standard error
("nocurv"), and its fitted values were NaN without a word. That test now
allows the warning, with the reason in a comment.

Files changed in this step, for the re-merge: `R/se-check.R`
(`jc_nonest()`, `jc_nan_lost()`, `jc_nonest_null()`),
`tests/testthat/test-se-limit.R`, `tests/testthat/test-mi-families.R`.

Runs:

- OpenBLAS 0.3.26 and 0.3.32 with LAPACK (`dev/boundfam-fixrun5.sh`,
  `dev/boundfam-out/ob-*`, every `lib:` line the lane's core).
  test-se-limit 240, test-se-tier 112, test-predict-newdata 14,
  test-predict-re-uncertainty 71, test-linpred-layers 83 and
  test-predict-lp-basis 30 each read fail=0 err=0 warn=0 on both
  versions. `dev/ciharden-out/` deleted.
- Reference BLAS, core, 318 files one per process
  (`dev/boundfam-suite6-lane/`): pass=23727 fail=0 err=0 skip=239
  warn=1. The one escaped warning was test-mi-families.R's (above).
  After the test change it reads pass=194 warn=0
  (`dev/boundfam-suite6-ext/`).
- Every extension, one file per process (`dev/boundfam-suite6-ext/`):
  coupling 551, eam 1851, latent 687, learn 474, ode 549, sample 4050,
  spatial 478, survival 556 passes. Every one has fail=0 err=0 warn=0.
- `R CMD check --as-cran`, core: Status: 1 NOTE (V8), tests 19 min.
  `dev/boundfam-check/` deleted.
