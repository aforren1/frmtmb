## Filed at the 0.75.0 release (2026-10-24)

Four lanes off 18b536fc (frmtmb 0.74.0), merged onto 4fd134db
(`dev/round-20261024.md`). What each lane, its reviewer or the
consolidation filed and did not fix, with the record to read. Nothing
was posted upstream.

### The next round's first lane: the `mo()` coding defect

- **MAJOR. `mo()` of an integer that does not start at 0 has one
  simplex cell too many per unit of the minimum, so frmtmb fits another
  model than brms on 1..k scales.** brms 2.23.0 subtracts the minimum
  (`get_mo_values()`, brms `R/formula-sp.R` lines 603 to 627): values 1
  to 4 have `Jmo = 3`, `Xmo_1` in 0 to 3. frmtmb sizes the simplex
  from `max(x)` (`R/frame.R`, `D_mo <- max(codes)`), since 0.15.0
  (65fed080, 2026-08-31, the commit that added `mo()`), and refuses
  negative integers, which brms takes. Measured (the brmsnames review,
  `dev/review-brmsnames-mo.R` in the lane's tree, seed 5, n = 300):
  `mo(k0)` on 0 to 3 has logLik -427.181869 (df 5), Intercept 2.0834
  (SE 0.1138); `mo(k)` on 1 to 4 has the same logLik with df 6 and NaN
  standard errors; `mo(k2)` on 2 to 5 df 7, NaN. Point predictions
  agree (to 1e-6), but the ML coefficients are not identified, 3 of 4
  predictions lose their standard errors, and `frm_sample()` samples
  another model than brms with only ESS warnings (Intercept 0.912, sd
  1.330, against brms's 2.098, 0.117; a cell `simo_mok1[1]` with
  posterior mean 0.176 on a category no row is in; 2 chains, 1000
  iterations, seed 11). Likert items start at 1, so this is a wrong
  Bayesian model on the common input. The fix: subtract the minimum,
  store it for newdata as brms does (`attr(x, "min")`), and take
  negative integers. Record: `dev/brmsnames-1024-findings.md`,
  "Defects found and not fixed", item 1;
  `dev/reviews/2026-10-24-brmsnames.md`, finding 7. The user asked for
  this to be the next round's first lane.

### Default priors and pointwise likelihoods of residual terms

- **`ar()` and `cosy()` default priors are flat on the internal
  parameter.** frmtmb.sample's default for an `ar` or `cosy` coefficient
  is flat on its unconstrained parameter, where brms's is flat on the
  natural one (no Jacobian in `default_priors_for()`). `sar()` got
  brms's flat-on-its-interval default through the `bounds =` slot in
  this round, and core already has the logistic map and its Jacobian
  for `cosy` (`ac_trans_value()`, `ac_trans_logjac()`), so the fix is
  small. Record: `dev/sarfcor-findings.md`, "What needs the user";
  `dev/reviews/2026-10-24-sarfcor.md`, finding 2.
- **One `ar()` parameter has two prior slots.** `default_prior(y ~
  ar(time = t))` lists `class = "ar"` and `class = "theta", coef =
  "thetaac_1"`, the same parameter in two spellings;
  `validate_prior()` and `frm()` take a prior on each, and the fit
  applies both. brms has the `ar` row only. Pre-existing (0.74.0 does
  the same; the release review's `dev/review-rel075-log/ardouble.txt`,
  finding 3). The 0.75.0 release fixed it for `sar()`'s `lagsar` and
  `errorsar`, the parameters with `bounds` (`dev/rel075-release-fix.R`,
  item 3): their internal row is not listed and a prior on it is
  refused by the class's name. The same for `ar`, `ma`, `cosy` and
  `cortime` needs their positions within `thetaac`
  (`autocor_class_idx()`), and a test that `class = "theta"` still
  bounds an internal parameter where no class reaches it.
- **`log_lik()` for residual covariance terms.** brms gives the
  leave-one-out conditional density per row for `sar()`, `fcor()` and
  `ar(cov = TRUE)` (`log_lik_gaussian_lagsar()`,
  `log_lik_gaussian_time()`, after the loo vignette on non-factorized
  models); frmtmb refuses `log_lik()` for every residual covariance
  term (`ar(cov = TRUE)`, `cosy()`, `unstr()`, `fcor()`, `sar()`). The
  review's estimate of the cost: per draw `Q = K'K / sigma^2` with `K`
  the Jacobian of `whiten()`, `g = Q e`, `yloo = y - g / diag(Q)`,
  `sdloo = 1 / sqrt(diag(Q))`, and brms's `student_t_cov_factor()` for
  student; 80 to 120 lines in core plus frmtmb.sample's plumbing and
  tests. Seven bin-2 ported rows read it. Record:
  `dev/reviews/2026-10-24-sarfcor.md`, finding 12;
  `dev/sarfcor-findings.md`, A5.

### brms names: what the brmsnames lane found and left

Record: `dev/brmsnames-1024-findings.md`, "Defects found and not
fixed", items 2 to 6, and "What I did not do". None is a defect of the
names the lane changed.

- `s(exp(x))` is refused ("The model uses `exp(x)`, which is not a
  column of `data`"); brms fits it (`sexpx_1`, `sds_sexpx_1`). The
  smooth's frame variables are taken as names (`dpar_frame_rhs()`).
- `conditional_effects()` picks other default effects than brms on
  `exp(x):z`, `t2(x, z)`, `mo(k):z`, `cs(x)`, `sigma ~ exp(z)`, `mi()`
  and multivariate models.
- `ranef()` of a smooth-only model returns the smooth's penalized
  coefficients under the internal names `s(x).1`, ...; brms stops, and
  names them `s_sx_1[1]`, ... in `variables()`.
- `variables()` of an ML fit with `mo()` has no `simo_` entries; brms's
  has `simo_mok1[1]`, ..., and frmtmb.sample's draws have them.
- The column order of `variables()` and of the draws is frmtmb's
  template order, not brms's (every centered intercept first, then
  response by response); recorded by the earlier lane
  (`dev/brmsnames-findings.md`, "Found and NOT fixed", item 1). The
  summaries follow brms since 0.75.0: `fixef()` of the fit and of the
  draws.
- The prior table's classes sort as brms's do since 0.75.0; `confint()`
  keeps frmtmb's internal parameter rows on purpose (brms has no
  `confint()` method).
- `y ~ mo(k) + mok` with a data column `mok`: frmtmb refuses the
  duplicate name, brms lists one `b mok` row (the review's
  `dev/review-brmsnames-mok.R`). frmtmb's refusal is the safer one.

### `bhaz()` and `cox()`

Record: `dev/bhaz-findings.md`, sections 4 and 7;
`dev/reviews/2026-10-24-bhaz.md`.

- `bhaz(gr = )`, brms's stratified baseline (one simplex per level), is
  refused by name. It needs a data path for an option term, one simplex
  per level in the simplex seam, `class = "sbhaz"` with `group =`, and
  frmtmb.sample's names and default (`sbhaz[<group>, k]`).
- `bhaz()` evaluates its arguments in the formula's environment and
  does not read `data2`; brms checks that a flag's variables are in
  `data` or `data2`, then evaluates the flag in the formula's
  environment, so a variable only in `data2` fails in brms too.
- The refusal of boundary knots that do not cover the data checks `y`
  only; an interval-censoring upper bound past the knots is clamped to
  the knot, as at 0.74.0 for `cens_y2` and truncation bounds, so "must
  cover every observed time" says less than it seems to (the review's
  "Cases the lane did not test").
- Without frmtmb.survival, `bhaz()` on any family says to load
  frmtmb.survival, and with it loaded the same call on gaussian is
  refused; `dec()` has the same shape with frmtmb.eam (review finding
  10, left by instruction).
- The generic `frm_compat()` rows of an addition term now describe a
  data term and an option term; no status changed
  (`dev/rel075-compat.R`).

### `sar()`, `fcor()` and the covariance-structure seam

Record: `dev/sarfcor-findings.md` (A6, B4, "What needs the user");
`dev/reviews/2026-10-24-sarfcor.md`.

- `rr()`'s factor-space expansion, `gp()`/`hsgp()` as specials and the
  Student-t sharing of `us_t`/`diag_t` are still not on the published
  entry shape: measured (41, 65 and 38 code lines; 3, 5 and 6 commits in
  five rounds), not ported, since each would change the tape of every
  fit that uses it.
- `sar()` refuses a matrix with complex eigenvalues beyond rounding (an
  asymmetric kNN graph, common in spdep). brms fails too; a general
  complex log-determinant (`log|1 - rho lambda|` over conjugate pairs)
  is cheap and not built.
- `confint_varcorr()` gives rho a delta-method interval on the raw
  scale, so its upper end can pass `1 / max eigenvalue` near the bound.
- `sar()` with `mi()` or with `ar(cov = FALSE)`, and `newdata2`, are
  refused (brms fits the first two). A complete-graph W warns; its
  likelihood has no maximum with an intercept.
- brms 2.23.0 defects, now `dev/upstream-bugs.md` brms-37 to brms-39
  (complex eigenvalues of a row-standardized W; `sar()` with weights,
  se, cens or trunc writes code stanc rejects; `posterior_epred()`
  fails on every error SAR fit), plus brms-36 (interval-censored
  `cox()`, from the bhaz lane): not posted. The user decides.
- covr could not install frmtmb.spatial from inside an R process on
  this machine ("lazy loading failed", 13 attempts, the base commit's
  package too), while the review's single attempts succeeded; an R
  child here sometimes ends without its last line. Worth its own look.

### The saved-fit refusal and its CI step

Record: `dev/smallfix-findings.md`, "Not done, and for the user".

- frmtmb.sample's `posterior_predict(newdata =)` refusal
  (`R/methods-draws.R`) keeps the old three-case text that core
  replaced; wording it per family needs `fam_sim_ctx()`, which is not
  exported.
- The skip gate on GitHub's runners is unverified: macOS's bash 3.2,
  the Windows `find`, and the pkgdown job's saved-fit step on the
  runner's library layout. The first push will show them.
- `recover_data.frmtmb_fit()` (`R/interop.R`) turns every error into
  emmeans's character refusal; `dev/smallfix-handlers.R` lists it on
  the merged build. FIXED at the release for the missing-package case:
  smallfix's head check there was raised inside emmeans's `try()`,
  which replaced it with "Perhaps a 'data' or 'params' argument is
  needed" (the release review, finding 1). The head check is removed;
  `emm_basis.frmtmb_fit()`, outside `try()`, refuses by the package's
  name, and `ref_grid()` and `emmeans()` are in the saved-fit test.
  Still open: an error raised inside `recover_data()` for another
  reason reaches the user as text, not as its class.
- The pkgdown deploy depends on the saved-fit step: a failure, or a
  runner killed during it, stops the site build (a choice, recorded in
  the workflow).

### Found at the consolidation

- **No single model holds `cox()` and `sar()`**: `sar()` takes a
  gaussian or student response, and `cox()` is refused in a
  multivariate fit ("Families with extra parameters ('cox') are not
  supported in multivariate fits yet"; `dev/rel075-mixed-probe.R`). The
  saved-fit refusal across three seams is tested with a family, a
  structure and a residual term from three absent packages in process,
  and with a real `car()` and `sar()` fit in child processes
  (`test-extension-home-mixed.R`).
- `test-guardrail1.R`'s block "a column named like an extension's
  special" checks the first home only (`car`) and drops its assertion
  without a word when frmtmb.spatial has registered it: 36 assertions
  with every extension loaded first, 37 otherwise (the release review,
  finding 7). No defect; it could check every home that is not
  registered.
- A merge rerun restores `NAMESPACE` and `man/` to the HEAD blob, so an
  install without roxygen after it builds the old exports; the release
  library is rebuilt in full (`dev/rel075-install.sh`) after every
  merge rerun.

### Closed at 0.75.0, from the 0.74.0 section above

- `bhaz()` is accepted with frmtmb.survival loaded, and
  `cox(bhaz = list(df = 6))` is refused with the spellings to use.
- The coefficient of `exp(x)` is `expx`, as in brms.
- frmtmb.sample reads `frmtmb_covstruct_info()` instead of its own list
  of field structures; the weights check reads each structure's
  `full_rank`; autoscale says which random slopes it did not rescale.
- The `sar()` refusal is gone (`sar()` is frmtmb.spatial's, `fcor()`
  core's), and the supported-terms sentence is built from the registry.
- A saved fit that needs two extensions names both in one refusal.
- `helper-absent-package.R` has three copies and its master, checked by
  `absent-helper-identical.yaml`.
- The `mixture_mvn()` tolerances are ratios; `test-prior-dpar.R` uses
  exported functions only; the `predict(newdata =)` refusal names the
  model's reason.
- `set_prior("dirichlet(2)", ...)` is refused at one point however it
  is written; `print()` of a curve takes `digits =`.
