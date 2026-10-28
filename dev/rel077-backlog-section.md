## Filed at the 0.77.0 release (2026-10-28)

Four lanes off 435e7ce6 (frmtmb 0.76.0), merged onto 435e7ce6
(`dev/round-20261028.md`). What each lane, its reviewer or the
consolidation filed and did not fix, with the record to read. Nothing
was posted upstream. The two items of the section above (the `se()`
parity gap and saved 0.76.0 draws) stand.

### Defects found, not fixed

- **MAJOR. frmtmb.sample's `posterior_predict()` of `set_rescor(TRUE)`
  draws draws the responses independently, and without `resp` returns
  one response.** Seed 3, 300 rows, residual correlation 0.8:
  `posterior_predict(ds)` is 300 x 300, one response; the
  cross-response draw correlations at rows 1 to 5 are -0.052, -0.040,
  -0.005, 0.082 and 0.026, where core's `predict()` gives 0.795 to
  0.861. `log_lik(ds)` is exact (1.6e-15). The same on rellib-r14.
  Repro: `dev/review-mointerop-rescorpp.R` (copied from the mointerop
  tree with `dev/review-mointerop-common.R`, which it sources by the
  lane tree's path); review `dev/reviews/2026-10-27-mointerop.md`,
  finding 8 and "The rescor record". The ported rows
  `posterior_predict:69` and `:72` pass through core's `predict()` and
  say nothing about this path (their ledger note).
- **`predictions()` on a multivariate fit without `resp` returns the
  first response only.** Seed 3, 200 rows, `mvbf(y1 ~ x, y2 ~ x)`: 3
  rows equal to `fitted()`'s `y1`, no `group` column, no `y2`. It should
  return every response keyed by `group`, as for a categorical outcome,
  or refuse without `resp`. Pre-existing (rellib-r14 the same).
  Review `dev/reviews/2026-10-27-mointerop.md`, short check S3
  (`dev/review-mointerop-mv.R`).
- **`frm_sample(formula, family = gen_extreme_value())` fails** with
  "frm_linpred() needs a fitted model"; the route from a fit works.
  Pre-existing (rellib-r14 too). Review
  `dev/reviews/2026-10-27-priordraws.md`, new finding 13
  (`dev/review-priordraws-gevbase.R`).
- **`default_prior()` fails on `rr()`**: `y ~ x + rr(0 + tf | g, d = 1)`
  stops with "arguments imply differing number of rows: 1, 0", on both
  arms. Review `dev/reviews/2026-10-27-reslot.md`, finding 6
  (`dev/review-reslot-p3b.R`).
- **`fitted()` gives `Est.Error` 0, not `NaN`, at an empty category**
  of the degenerate `(mo(x) | g)` gap fit (seed 23), both arms.
  `dev/mointerop-findings.md`, fix round 1, "Pre-existing".
- **`drop1()` omits `mo()` terms**, both arms, also on a 4-category
  fit. Same record.
- **`predict()` of a gap fit draws at the estimates alone** (the
  parameter uncertainty is dropped), both arms. Same record.
- **The marginaleffects `re_formula` warning (upstream).**
  marginaleffects 1.0.0 warns that `re_formula`, `allow_new_levels`,
  `dpar` and `resp` are "not known to be supported" for a `frmtmb_fit`,
  because it keeps a fixed list of known arguments per class; they work.
  A feature request (a way for a package to declare its prediction
  arguments) would remove it. Not posted. `dev/mointerop-findings.md`,
  "Defects found, not fixed", 1.
- **`frm_sample()`'s default-prior message prints `b (flat), as brms
  leaves slopes`** when the call gives class `"b"` a prior; the prior
  applies, the line is wrong. Pre-existing (`dev/priordraws-disclosure.R`,
  `dev/priordraws-findings.md`, section 8).
- **The sampling table shows `lb` NA on `sd` and `sigma` rows**, where
  brms's `default_prior()` shows 0, on both routes. `dev/reslot-findings.md`,
  "Not done, and why".
- **`vignettes/brms-migration.Rmd`, the paragraph on prior classes**,
  still says `sdgp` and `lscale` are refused by name and that a `coef`
  on an `sd` row is refused; both have been wrong since 0.72.0 or
  earlier. Same record.
- **`summary()` and `variables()` of saved 0.76.0 draws are not
  refused** for a model with a `cov = TRUE` `ar()`, `ma()` or `arma()`
  (the old marginal sd is reported as `sigma`); only the computed
  quantities refuse. The section above has the item and the seam it
  needs.
- **`confint(parm = "x")` is refused on a saved `ar(cov = TRUE)`
  fit**, though a Wald interval for a mean coefficient does not read
  sigma. Consistent with NEWS. Review `dev/reviews/2026-10-27-sigmainn.md`,
  short check S2: check `parm`, or leave it.
- **The plain-data-frame mark of the 0.77.0 release fix R1 can refuse
  an intended `(flat)`.** `as.brmsprior()` of a data frame that holds
  frmtmb's printed `(flat)` (after the release review, on any row,
  whatever its `source` column says) is marked, and frm_sample()
  refuses it where a `(flat)` row replaces one of brms's proper
  defaults. Each of these is refused, with the message, where a brms
  user may have meant brms's `(flat)` (the release review, finding 8,
  `dev/review-rel077-b-run/r1-route.txt`, lines 57 to 63):
  `as.brmsprior(as.data.frame(<route = "sample" table>)) +
  set_prior("", class = "sigma")`; the same coerced table with
  `source = "user"` on the edited row; a coerced sample-route table
  edited to `""` on `sd`; and a data frame built by hand with `(flat)`
  on sigma. The message says to pass only the specifications you set.
  The cure is a route that survives `as.data.frame()`, which
  `test-get-prior-route.R` ("coercion drops the class and the label
  together") rules out today. `dev/round-20261028.md`, "Release fixes",
  R1 and F1.
- **`frm()` refuses `pe[, "prior", drop = FALSE]` with the wrong
  cause.** A one-column selection of an edited fit-route table keeps
  its route since release fix F1, and `frm()` then stops with
  "Duplicated prior specifications are not allowed: 'b' is given 10
  times": with only the `prior` column every row reads as class `b`.
  0.76.0 fitted it (the column selection dropped the route and the
  class, and the plain rows read the same way). A one-column table
  names no slot, so a refusal that says so would be right. Release
  review, short check of release fix 2, S3
  (`dev/review-rel077-c-run/fix2-r15.txt`, seed 11).
- **`default_prior() + set_prior()` and `set_prior() + default_prior()`
  fail** in R's `Ops` dispatch ("Incompatible methods" between
  `Ops.data.frame` and `+.frmtmb_priorlist`, then "non-numeric argument
  to binary operator"), on 0.76.0 too; brms 2.23.0's `+` combines them.
  `c()` takes the table in either order since the release review
  (finding 3). Release review, `dev/review-rel077-b-run/r1c-rellib-r15.txt`.
- **PARITY GAP: a cluster sandwich over the series of a residual
  correlation term.** Under `ar(time = t, gr = h, cov = TRUE)`,
  marginaleffects' `vcov = ~ h` is refused by `vcov_cluster()` with
  every other sandwich spelling ("does not support the residual
  correlation term ... Use frm_bootstrap()"). The scores do factor by
  series there, so `~ <the ar group>` (or a grouping that nests it)
  has a valid cluster sandwich. Release review, finding 2
  (`dev/review-rel077-b-run/x3.R`, `x3.txt`).
- **A residual-correlation term of order 2 or 3 under `prior = "flat"`**:
  `prior_summary()` prints `ar` as `uniform(-2, 2)` or `uniform(-3, 3)`
  with no bounds, not as `(flat)`. The round trip gives identical draws,
  and brms has no such model. Release review, finding 8
  (`dev/review-rel077-b-run/x1.txt`, line 52).
- **Release fix R2's message**: it reads "a list of 1 set_prior()
  specifications" for `list(sp1)` and still suggests `sp1 + sp2` there;
  `list("flat")` and `list()` get the generic message. Release review,
  finding 8 (`s1r2-rellib-r15.txt`).
- **`frm(prior = "flat")` is refused with the generic message**, where
  `frm_sample()` takes the string; a pointer to `frm()`'s
  maximum-likelihood default would help. Pre-existing. Release review,
  finding 8 (`s1r2-rellib-r15.txt`, line 24).
- **`predictions(newdata =)` on an `ar(cov = FALSE)` fit stops behind
  marginaleffects' generic message** ("newdata repeats a time within a
  group"), on rellib-r14 too; and the `"HAC"` refusal says a frmtmb fit
  "does not carry" a time order, but an `ar()` fit does. Release review,
  finding 8 (`x4.txt`).
- **The bin-2 adapters find brms's suite through `../../brms-suite`**
  when `FRMTMB_PORT_ROOT` is unset: a stage folder at another depth
  then fails `distributions:446`, `:479`, `:527` and `:549` with "object
  'ndraws_vec' not found". With the variable set, as
  `dev/brmsport-bin2-all.sh` sets it, the file passes 92. Release
  review, finding 8.
- **`check-frmtmb-sample.yaml` calls the skip gate with no allow-list**
  and its install steps name no r-universe repository, yet RTMBode
  arrives in CI (run 40182637732). If it stops arriving,
  `test-compat-preflight.R:336` skips and the gate fails the job.
  Pre-existing; the same cause as the local skip-gate failure. Release
  review, finding 8.
- **The record `dev/round-20261028.md` was the one changed file stored
  LF** (every other is CRLF; harmless with `core.autocrlf=true`).
  Release review, finding 8.
- **Optimizer precision on `mo()` fits**, pre-existing:
  `frm(y ~ mo(x) + z)` on two values stops at "relative convergence (4)"
  about 2.5e-6 relative from `lm()`'s coefficients, where
  `frm(y ~ I(x - 1) + z)` reaches them; `newton_polish()` does not move
  it. Statistically nil (1e-5 standard errors). `dev/mointerop-findings.md`,
  "Defects found, not fixed", 2 (`dev/mointerop-twofit.R`).

### Not done, with the reason

- **brms's `meanme` and `sdme` prior classes** (a user question): brms
  needs proper priors there for `sample_prior = "only"` on an `me()`
  model; frmtmb has no prior slot for them, so `"only"` refuses `me()`
  models naming the two classes. The change touches `set_prior()`'s
  class lists, the prior table, the resolver, the brmsprior translation
  and frm_sample()'s defaults, all in `R/priors.R`. Decide whether to
  add them. `dev/priordraws-findings.md`, section 14, item 5.
- **Prior draws for `gp()` terms.** brms 2.23.0 names them
  `prior_sdgp` and `prior_lscale__1_gpx` by an accident of its renaming,
  and a second gp term would give duplicate names. Not emulated; no
  column. `dev/priordraws-findings.md`, section 8.
- **`car()`'s dependence parameter is flat on its logit** under the
  defaults and `prior = "flat"`, not brms's `(flat)` on (0, 1): it
  needs a seam slot for a registered structure's (flat) chart. Under
  `prior = "flat"` such a posterior is improper. Documented in
  `?frm_sample` and NEWS. `dev/reslot-findings.md`; review
  `dev/reviews/2026-10-27-reslot.md`, finding 5.
- **Two random-effect terms on one grouping factor** with different
  ids (`(1 | g) + (0 + x | g)`) get one `prior_sd_g` statement here;
  brms makes one per id, both named `prior_sd_g`. `dev/priordraws-findings.md`,
  section 8.
- **`prior_samples()` is not defined in frmtmb.sample**; brms's
  deprecated `prior_samples()` reaches `prior_draws()` while brms is
  loaded. Same record.
- **`variables()` of an ML fit lacks `simo_` entries**: the hypothesis
  layout would need the simplex and its delta-method rows.
  `brmsfit-methods:978` stays a divergence anyway. `dev/mointerop-findings.md`,
  section 7.
- **The ported helper `brms_shim_update()` drops `sample_prior`** with
  the sampling arguments, so a ported `update(fit, sample_prior = )` row
  would not reach it. No such row is ported. `dev/priordraws-findings.md`,
  section 8.
- **`dev/brmsport-bin2-survey.tsv` still carries "not ported (cost)"**
  for the four rescor rows; the ledger's run overrides it. Same as the
  0.76.0 survey note on `sar()`. `dev/mointerop-findings.md`, "Defects
  found, not fixed", 3.
