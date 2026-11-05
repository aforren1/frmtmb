## Filed at the 0.79.0 release (2026-11-04)

Seven lanes off d5d0056d (frmtmb 0.78.0), merged onto d5d0056d
(`dev/round-20261104.md`). What each lane, its reviewer or the consolidation
filed and did not fix, with the record to read. Nothing was posted upstream.
The list in `dev/upstream-bugs.md` is renumbered: brms-45 to brms-50, RTMBdist-5
and TMB-5 are not reported.

**Closed by this release**: the survey's D1 to D8 of the 0.78.0 section
(`dev/equiv-findings.md`, section 4), and the exact `gp()` L1 gap, which was the
translator's order of the positions.

### MAJOR

- **First of the next round. `marginaleffects`' `predictions()` and
  `avg_slopes()` standard errors ignore variance components on every model with
  a group-level term.** `avg_predictions()` is off on 77 percent of 168 fits,
  `avg_slopes()` on 52 percent, levels on 67 percent; poisson with `(1 | g)` reads
  .0774 against .2019. The warning `frmtmb_interop_se_part` states it and does not
  fix it. `dev/review-predread-r2-mfxrate.R`.
- **poisson `x + 1e6`: the Laplace slope is off by 5.4 to 9.3 standard errors,
  silent on 2 of 3.** `dev/reviews/2026-11-02-impsamp.md`.
- **bernoulli `residuals_osa` near -3.7.** `dev/reviews/2026-11-02-latres.md`.
- **`b1 (1 - exp(-b2 x))` under gaussian is about 70 nats short at code 0, no
  warning, 5 of 5 seeds; 16 of 450 nl rows are short of random starts (22 before
  the release's restart fix).** -601.79
  against -530.47 at seed 101. brms's 6 starts miss it too.
  `dev/reviews/2026-11-02-nlfix.md`, "Filed, no change asked".
- **`ar()` on a step field is 25.9 nats short, student -11, under the exact
  `gp()`.** Named in NEWS, Limits, by decision. `dev/reviews/2026-11-03-ncgp.md`.

### MINOR and performance

- `fitted()`, `residuals()`, `frm_linpred(se.fit)` on a latent residual are not
  sparse: 984 s at 8000 rows.
- `fitted()` on a new RE level reads below the mean by a factor of 1.649 at sd 1.
- `a * x` with x near 1e-20 is silent: -299.33 against -272.50 on both builds
  (`dev/review-nlfix-c9.R`, seed 6201).
- `avg_slopes()` standard error on a gaussian field.
- The first-order Laplace correction's shift could fit, not warn: median log
  `sderr` error 0.106 to 0.010 on fresh designs. A feature for a later round.
- The marginal tape's Hessian for the SE check: 0.87 s against 11.3 s at 500
  positions.
- A `disc` saddle beside a `mu` ridge: -377.03 at code 1, with its warning, against
  -371.54 from a start.
- `conditional_effects()`'s band for a category under `me()`: 0.197 to 0.296 on
  category 3.
- `se__` of `conditional_effects()` is on the link scale on lognormal and hurdle
  fits, 0.068 against brms's 0.519, and 0.102 against 0.809 on poisson.
- **Move latres's `mi()` refusal after `mi_ml_frame()`**, to fire on
  `length(rows) > 0`, against a test. Its reason is false when every missing row
  is unread, and nobody has measured the case with read rows.
- Nested quadrature: an exact route for more than one RE block
  (0.73.0 refuses it).
- A me() term on another predictor than `mu` beside a latent residual on `mu`:
  one node cannot read two independent draws, and the latent's variance does not
  join that node (found at this consolidation; `R/me.R`, `me_row_loads()`).

### From the lanes' records

- refuse: the `lb = 0` of the `sd` and `sigma` class rows (filed at 0.78.0); a
  translator rule for `sigma ~ s(z)` and for `tmp_xi`; class `sd` on `car()` and
  `spde()` keeps its old spelling, `sdcar` stays refused by name.
- impsamp: crossed factors under `importance =` stay refused; a heavy-tail guard
  for hard binary designs (sd 0.67 at 512 draws); poisson `x + 1e5` false
  convergence with NaN SEs, tape condition number 9e18; MCSE by pairs.
- latres: `variables()` on draws lists `latent_state[k]`, where brms has `zerr[n]`
  and `err[n]`; a corrected estimate and logLik (Laplace + C, `theta + shift`);
  `get_joint_cov()` on a new row's design, which reads `V` of the outer parameters
  alone, and M7.
- nlfix: `simulate()` on `beta_binomial()` stops if a fitted mean rounds to 0 or 1;
  `frm_allfit(start_from_mle = FALSE)` reports bobyqa "failed" on a body's fit from
  zeros; `Beta()`, `dirichlet()` and the inflated Beta forms keep raw `lgamma()`
  differences; `a + b * x` at y about 1e6 is -3045.62 against -1658.76 linear, at
  code 0; `frm_curve()` on a linear `cov = FALSE` is an internal error, and on a
  body it omits the shift.
- predread: `posterior_predict()` on draws gives one R warning per draw on a row of
  `NA`; `frm_lp_basis()$extra_var` is finite on a row of `NA` in `x`; `predict()`'s
  `Estimate` reads `NaN` there; `predictions()` gives 9 to 13 warnings of the `me()`
  class, one per perturbation.
- ncgp: the gaussian closed form (a recommendation); the drawn route at `sigma` at
  or below 1e-3 on small n stays unconverged at parity with brms, and
  `laplace = TRUE` is the route there.
