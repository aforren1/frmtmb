## Filed at the 0.78.0 release (2026-10-29)

Three lanes off ac4906bd (frmtmb 0.77.0 plus one test fix), merged onto
ac4906bd (`dev/round-20261029.md`). What each lane, its reviewer or the
consolidation filed and did not fix, with the record to read. Nothing
was posted upstream (brms-43 and brms-44 are in `dev/upstream-bugs.md`,
not reported).

**Closed by this release**: the MAJOR item of the 0.77.0 section,
frmtmb.sample's `posterior_predict()` of `set_rescor(TRUE)` draws, which
drew the responses independently and returned one response without
`resp`. It now draws them jointly through core's `predict_rescor_draw()`
(`dev/rescorpp-findings.md`; review `dev/reviews/2026-10-28-rescorpp.md`).
The ledger notes of `posterior_predict:69` and `:72` say so.

### The equivalence survey's defects (`dev/equiv-findings.md`, section 4)

The survey measured 172 feature pairs against brms 2.23.0 (headline
83.6 percent, Wald [78.2, 89.0], reach-weighted). Its defects, each with
a repro script, seed and brms comparison; none is fixed in 0.78.0.

- **MAJOR. D1, silent: `nl = TRUE` with an ordinal family drops the
  nonlinear body.** cumulative, sratio, cratio, acat and
  hurdle_cumulative: frm() fits the thresholds alone, without a word. On
  p049's data the logLik is -267.12, the intercept-only fit, against
  -255.52 for the same predictor written linearly and for brms's
  maximum. The compat table says nl x every ordinal family is refused,
  and frm() does not refuse. Repro `dev/equiv-repro-thresnl.R` (seed
  49); `dev/equiv-nl-sweep.R` (seed 20261029) adds hurdle_cumulative.
  Refuse it, or honor the body.
- **D2, loud: brms's `class = "sds"` prior on a smooth is refused**
  (p043). frmtmb's own spelling, `class = "sd"` with `group = "s(x)"`,
  is an L1 identity with brms's `sds` prior (constant -0.667 at three
  points, gradient 1.6e-13; `dev/equiv-variant.R` with
  `dev/equiv-variant-d2.R`). Fix: alias `"sds"` to that spelling.
- **D3, loud: `importance =` stops with "did not converge"** on the
  models it exists for (p124 binomial, a works row; p069
  `(1 + x | g)` poisson). binomial, 10 groups of 20: `importance = 200`
  and `512` stop ("the corrected negative log-likelihood ROSE"), 1000
  fits. `dev/equiv-repro-importance.R`, seeds 124 and 69.
- **MAJOR. D4, silent: `frm_curve()` and `frm_linpred(type = "link")`
  drop the `cs()` offsets** (p013). They return the shared eta; brms's
  `posterior_linpred()` at the same values returns eta plus the
  per-threshold offset, and frmtmb's own `fitted(scale = "linear")`
  already matches brms digit for digit. `dev/equiv-repro-cscurve.R`,
  seed 13.
- **D5, loud: `mm()` with an `|ID|` key, and `mm(cov = A)`, are
  refused** (p070, p087); brms 2.23.0 builds and compiles both
  (`dev/equiv-repro-refusals.R`). The reasons live only in compat notes,
  and the sampled rows also say conditional, so the rows are wrong too.
- **D6, loud: `arma(cov = TRUE)` on a poisson response is refused**
  (p167); brms builds the latent-residual model. The same refusal
  probably covers every autocor term on every family but gaussian and
  student; only this pair was measured.
- **D7, loud: `bf(..., nl = TRUE, autocor = ~ fcor(M))` is refused**
  (p158), brms's spelling, which builds and compiles.
- **D8, loud: an `mi()` response under von_mises is refused** (p169).
  NEWS 0.16.0 states the limit without a reason.
- **frm()'s default starts fail for a nonlinear coefficient under
  Gamma's and inverse.gaussian's default links.** A nonlinear
  coefficient starts at 0, where the inverse link and 1/mu^2 are
  undefined: `frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family =
  Gamma())` stops with "The nonlinear fit failed from its default
  starting values (NA/NaN gradient evaluation)", where the linear fit
  succeeds. From starts that keep mu positive the nl maximum equals the
  linear one (`dev/equiv-data/nl-sweep.tsv`; review short check 2,
  seeds 5101 and 5102). frm() has the data to start the body at the
  mean.
- **The exact `gp()` gap.** On the exact gp() pairs (p022, p023, p094,
  p109) the L1 joint differs from brms's by 94 to 1127 nats. Matching
  frmtmb's nugget to brms's size removes about 11 of p094's 94 nats; the
  remaining gap of about 83 nats (-83.81 at a nugget of 1e-9, -83.74 at
  1e-12) is unexplained. `dev/equiv-gp-gap.R` cannot attribute it in
  double precision (92 of 196 correlation eigenvalues below 1e-6). A
  candidate in the instrument's translator or in the model; the four
  pairs are scored cannot measure.
- **The 53 compat-row proposals** of the survey, NOT applied:
  `dev/equiv-data/compat-proposals.tsv` (id, pair, current status,
  proposal, reason; 14 from sampled pairs, 39 from the nl x family
  sweep). Among them the nl x family rule: the table refuses nl for every
  family but gaussian, and frm() enforces it only for the five
  multi-mu families; 35 honored rows are corrections, the 5 ordinal rows
  are D1.
- **The survey's own open items** (review finding 9 and the findings,
  section 8): the `bar_crossing` reach pattern counted nested grouping
  where the feature is the refused `x * (1 | g)` crossing; the pattern
  is corrected (83.5 with the corrected weights) but the sample is not
  redrawn. Redraw at the next frame change. An nl body cannot call
  Stan's `inv_logit()`, which brms's nonlinear vignette uses ("could not
  find function"). The translator gaps the survey did not close: mv
  mixtures, gev `xi`, the univariate `Ymi`, `mo()` outside `mu`, `cs()`
  in a mixture, `t2()` under `mvbf()`.

### Filed by the meme lane and its review

Review `dev/reviews/2026-10-28-meme.md`, "Pre-existing, outside the
lane"; repros `dev/meme-filed.R` (seed 5), `dev/meme-findings.md`,
section 7.4.

- **`update(<draws>, prior = )` replaces the first call's prior instead
  of merging it**, as brms's `update.brmsfit()` does
  (`dev/brms-suite/brms/R/update.R` lines 144 to 155). After
  `frm_sample(prior = meanme normal(0, 1) + sdme exponential(1))` then
  `update(prior = corme lkj(4))`, `meanme normal(0, 1)` and `sdme
  exponential(1)` are gone.
- **`variables()` orders the `me()` hyperparameters by term, brms by
  group**: on `me(xg, gr = g) + me(z) + me(wg, gr = g)` frmtmb gives
  `meanme_mexg meanme_mez meanme_mewg`, brms `meanme_mexg meanme_mewg
  meanme_mez`. `b_` and `bsp_` interleave differently on the dpar, mv
  and mi shapes. Same set of names on every shape.
- **`default_prior()` gives `sigma` no `lb = 0`**, on every shape;
  brms gives `lb = 0`. The same as the 0.77.0 item on the `sd` and
  `sigma` rows of the sampling table, still open.
- **`mi()` latent values are `miss_<k>` in draws**; brms names them
  `Ymi_<resp>[k]`.
- **`"only"` refuses a `me()` vector with a proper prior on some terms
  only**, naming `'meanme_1'`; brms samples it, flat on the rest. A
  refusal, consistent with class `b` (the 0.77.0 rule); the review
  accepted it. `dev/meme-findings.md`, section 5.
- **brms refuses a restated `lb = 0` on `sd` and `sdme` with a `coef`**,
  frm_sample() takes it (it bounds nothing). Harmless; review short
  check 1, S1.
- **`dev/feature-gaps.md` lines 68 and 423 are stale about `me()`**
  (not edited, historical). `dev/meme-findings.md`, section 5.

### Filed by the rescorpp lane and its review

Review `dev/reviews/2026-10-28-rescorpp.md`; `dev/rescorpp-findings.md`,
"Filed, not fixed".

- **`conditional_effects()` on multivariate draws** (pre-existing, the
  same on rellib-r15). Without `resp` it returns one effect, for the
  first response; brms returns one per response (`y1.y1_x`, `y2.y2_x`,
  `y3.y3_x`). With several responses it fails with a message that is
  not about the call: `resp = c("y1", "y2")` gives "'length = 2' in
  coercion to 'logical(1)'", and with `method = "posterior_predict"`
  "method = 'predict' needs a family with a simulator", false for a
  gaussian model. Fix: one panel set per response, as brms; until then,
  refuse several `resp` with a sentence that says so (review finding 7,
  `dev/review-rescorpp-readers.R`). The 0.78.0 release review (finding
  3 and finding 8, `dev/review-rel078-b-run/x1.txt`, `x3.txt`, `x5`) adds:
  on a `me()` model with rescor, `conditional_effects(resp = c("y",
  "y"))`, `resp = c("y2", "y2")` and `resp = c("y", "y2")` all stop with
  R's internal "'length = 2' in coercion to 'logical(1)'", where brms
  takes several responses; and **`pp_mixture()`** stops with the same
  internal message for `resp = c("y", "y")` and `c("y2", "y2")`. Neither
  reads `draws_resp_set()`, the shared reading of `resp` that refuses a
  name given twice. Fix both with it.
- **A stored `loo_R2` in saved 0.77.0 draws keeps the old weights**:
  draws saved on rellib-r15 after `add_criterion(c("loo", "loo_R2"))`
  return their stored values on 0.25.0 (`R2y2` -0.002062348), computed
  with joint-density weights for `resp`. frmtmb has no users, so no seam
  (review finding 9).
- **One-response matrices carry `NULL` dimnames where brms has
  `list(NULL, NULL)`** (pre-existing nit, review finding 12): every
  `resp = first` difference on `posterior_predict`, `predictive_error`,
  `predict(summary = FALSE)` and `log_lik` in the review's layout logs.
- **On a univariate model, `posterior_epred()` and `fitted()` with
  `resp = character(0)` stop with a raw "attempt to select less than one
  element"**, and with `resp = c("y", "y")` with "Unknown response:
  'yy'". `posterior_predict()` and `log_lik()` on the same model read
  `character(0)` as the response; brms ignores `resp` on a univariate
  model. rellib-r15 gives the same two errors (review, short check of
  fix round 1, `dev/review-rescorpp-uniresp.R`).
- **`log_lik(combine = FALSE)` without rescor is refused**; brms returns
  a draws x observations x responses array. `dev/rescorpp-findings.md`,
  section 8.
- **Two brms 2.23.0 defects**, `dev/upstream-bugs.md` brms-43 (one
  response of a student rescor model reads `nu = 1`) and brms-44 (a
  permuted `resp` stops every method that prepares several responses
  under rescor). Not reported.

### Found at the consolidation

- **Release fix R1** (meme's review, short check of fix round 2, S3):
  with both `set_prior(class = "Intercept", nlpar = "a")` and
  `set_prior(class = "b", nlpar = "a")` on `a ~ 1`, frm_sample() kept the
  later row. Fixed: the pair is resolved by specificity, as frm() does
  (`dev/round-20261029.md`, "Release fixes"). Its printed summary is
  fixed after the release review (finding 2): `prior_summary()` lists
  both rows as written.

### Filed from the 0.78.0 release review

`dev/reviews/2026-10-29-release.md`, findings 5 and 8; its evidence is in
`dev/review-rel078-b-run/` (gitignored). Each is the same on rellib-r15:
none is a regression of this release.

- **MAJOR. `me()` on new data is a plug-in: a silent narrowing of the
  predictive interval.** `posterior_epred()` and `posterior_predict()`
  with `newdata` take the observed `x` as the latent value and ignore
  `sdx`; brms 2.23.0 draws a new latent per draw from `normal(x, sdx)`
  (`dev/brms-suite/brms/R/prepare_predictions.R`, the branch "sample new
  values of latent variables"). Confirmed by
  `dev/rel078-fix2-me-newdata.R` (`dev/rel078-log/fix2-me-newdata.txt`;
  test-rescor-me.R's model, data seed 29, sampler seed 5, 800 draws): on
  row 1 the sd over draws of the expected value is 0.1658 on new data,
  the same 0.1658 with `sdx = 0`, and equal to `b0 + b x` on the draws to
  a difference of 0; brms's draw, `b0 + b (x + sdx z)` on the same draws,
  gives 0.2749, so frmtmb's sd is 0.603 of brms's and its 95 percent
  interval 0.634 wide against 1.036. The review measured 0.171 against
  about 0.27 on its own fit (`x3.txt`, part f). Nothing warns. Fix: draw
  the latent values on new data as brms does, or refuse `newdata` with a
  `me()` term; until then, say so in the help of both methods.
- **`bayes_R2()` and `loo_R2()` with `resp = character(0)` stop with
  R's internal "'data' must be of a vector type, was 'NULL'"**, where
  the other methods read an empty `resp` as every response, as brms's
  `validate_resp()` does (`posterior_epred(resp = character(0))` gives
  `y` and `y2`). Their `brms_resp_index()` matches nothing and builds an
  empty table. The same on rellib-r15, so older than this release
  (`dev/rel078-fix3-emptyresp.R`, `dev/rel078-log/fix3-emptyresp-*.txt`,
  data and sampler seed 3; the release review's short check of fix 2).
  Fix: read an empty `resp` as every response in `brms_resp_index()`.
- **`update(<draws>, prior = new)` replaces the first call's prior**;
  brms keeps the old user rows (`brms/R/update.R` lines 150 to 155). In
  `x8.txt`, `b_x` went from 9.932 to 2.094 after an `update()` that
  added only a sigma prior. `?update.frmtmb_draws` documents it, so it is
  a documented divergence; meme filed the same (section above).
- **Multivariate `cor` group rows carry `resp`**: the tables have one
  `cor group = g` row per response, where brms has one row with no
  `resp`; `set_prior("lkj(5)", class = "cor", group = "g", resp = "y2")`
  is accepted on both routes and reaches only y2's block, where brms
  refuses it (`x4-r16.txt`, `x6-*.txt`).
- **`set_prior("lkj(2)", class = "rescor", coef = "y__y2")` is
  accepted** on both routes; brms refuses it (`x4`).
- **`conditional_effects()` of a multivariate model without `resp`**
  returns the first response only (`x3`, `x5`); the rescorpp item above
  has the rest.
- **`default_prior()` table `+ set_prior()`** fails with "non-numeric
  argument to binary operator", an `Ops` method conflict; `c()` works.
  Filed at 0.77.0 too; still open.
- **Draws depend on the order of specifications in different slots**:
  the floating-point sums differ by order, pB + pb gives 4.788349 and pb
  + pB 4.797311 (`x5`). No action needed; recorded so a bitwise
  comparison of two orders is not read as a defect.
- **frmtmb.sample's Suggests lacks RTMBode.** `test-compat-preflight.R`
  skips at line 336 ("{RTMBode} is not installed") under `R CMD check
  --as-cran` only: RTMBode is installed in the user library, and the
  file passes 86 with 0 skips one file per process. frmtmb.sample
  suggests frmtmb.ode but not RTMBode, and `--as-cran` masks packages a
  package does not declare (inferred by the review, not run). The skip
  gate's allow-list for frmtmb.sample's coverage job hides it in CI. Fix:
  declare RTMBode in Suggests (finding 5).
