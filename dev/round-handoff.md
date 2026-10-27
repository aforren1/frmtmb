# Handing a round to a new session

Written 2026-09-24 at the 0.63.0 release, rewritten 2026-09-28 at the
0.64.0 brms-parity release, 2026-09-29 at the 0.65.0 and 0.66.0
releases, 2026-09-30 at the 0.67.0 release, 2026-10-06 at the
0.68.0 release, 2026-10-07 at the 0.69.0 release, 2026-10-13 at the
0.70.0 release, 2026-10-14 at the 0.71.0 release, 2026-10-20 at the
0.72.0 release, 2026-10-21 at the 0.73.0 release, 2026-10-23 at the
0.74.0 release, 2026-10-24 at the 0.75.0 release and 2026-10-27 at the
0.76.0 release. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: the library has been lost EIGHT times, and the last three
losses each followed processes being killed, not low disk.

## The 0.76.0 round, 2026-10-27

`dev/round-20261027.md` is this round's record. Each lane has its
findings (`dev/mocode-findings.md`, `dev/acprior-findings.md`,
`dev/loglik-findings.md`) and its review with the fix-round checks
(`dev/reviews/2026-10-27-<lane>.md`). Three lanes off `70476dc3`
(0.75.0, pushed), merged onto `70476dc3` in this order:

- `mocode`: `mo()` codes an integer from its minimum, as brms does
  (BREAKING since 0.15.0); fractional, missing, factor and out-of-range
  new data refused; `conditional_effects()` and emmeans hold `mo()` as
  brms does; prediction and marginaleffects standard errors include the
  simplex, `NaN` at an empty category; a saved 0.75.0 fit coded the old
  way is refused.
- `acprior`: frmtmb.sample's autocorrelation defaults flat on the
  natural parameter; a prior on class `sigma` under a first-order
  `cov = TRUE` `ar/ma/arma` read at brms's innovation sd; the cosy chart
  on brms's [0, 1]; declared bounds kept under a user prior; one prior
  slot per residual correlation parameter (BREAKING). New in the
  sampling API: `autocor_flat_priors()`, `autocor_bound_chart()`,
  `autocor_innov_sd()`.
- `loglik`: `log_lik()` of draws for residual covariance terms (brms's
  leave-one-out conditional) through `autocor_row_loglik()` and the
  seam slot `frmtmb_register_autocor(jacobian =)`, which `sar()` fills;
  the 13 bin-2 `sar()` and `fcor()` rows ported, 8 of them with the new
  ledger class "adapted-input".

Versions: frmtmb **0.76.0** (minor, breaking); frmtmb.sample
**0.23.0** and frmtmb.spatial **0.3.0** (minor), each needing frmtmb
0.76.0, measured by `dev/rel076-floor.sh`; survival 0.2.0, latent
0.7.1, eam 0.13.0, coupling 0.6.1, ode 0.7.2 and learn 0.7.3 unchanged.
The release library is `C:/Users/adf44/source/r/rellib-r14`;
`rellib-r13` keeps 0.75.0.

### Decisions made at the consolidation, for the coordinator

- **The sampling-API page keeps both lanes' names** (acprior's three,
  then loglik's one); the only textual conflicts outside NEWS.
- **`brmsfit-methods:924` is a divergence by policy**, as mocode's
  review asked: the block's newdata hold `Exp = 4` on every row, and
  brms builds that model only under testmode. Three neighbouring rows
  of the same block (923, 927, 937) keep their divergence outcome but
  now record the refused refit as their observation.
- **The bins 1 and 2 figure is 619 of 823**, not the 620 the brief
  expected: 607 + 13 (the bin-2 port) - 1 (row 924).
- **Two consolidation tests** in frmtmb.sample's `test-loglik-autocor.R`:
  `cosy()` with `mo()` and `ar(cov = TRUE)` through `frm_sample()` and
  `log_lik()`, against the dense conditional of each draw. The `sar()`
  prior tables with both changes are checked by a probe
  (`dev/rel076-overlap.R`), since no test file loads frmtmb.sample and
  frmtmb.spatial together.
- **The floors are measured without an install** (`dev/rel076-floor.sh`
  loads copies of rellib-r14's builds over rellib-r13), because this
  consolidation installs only into rellib-r14.
- **The lanes' dev/ scripts are tracked, the reviewers' are not** in
  the release tree, as at 0.75.0; the loglik reviewer's three output
  folders of ledger tables are left out of its lane commit too.

### The release review and its fixes

`dev/reviews/2026-10-27-release.md`: FIX, no blocker or major. The
coordinator asked for findings 1 to 5 and 7 as release fixes
(`dev/rel076-release-fix2.R`, sourced by the resolve step); finding 6
(three ledger reasons the run no longer reaches) is filed:

1. NEWS lists the refusal of a two-valued integer `mo()`.
2. The record's comparison of the two tmbstan objectives was invalid
   (tmbstan's process-global pointer, `dev/rtmb-pitfalls.md` item 23);
   the reviewer's measurement replaces it, and the backlog item is
   resolved: one ulp in the summation order, amplified by the sampler.
3. `ac_bound_theta()` takes a bound within 64 ulps of a computed end as
   that end, so `ub = 1` on a row-standardized W is accepted, and its
   refusal prints the ends in 15 digits.
4. `?frm_sample` says the sigma intercept's prior under `sigma ~ x` is
   read at the log marginal sd.
5. `log_lik(newdata =)` names `posterior_predict(newdata =)` only where
   it answers (`draws_pp_structured()`).
6. The release commit carries the review and its scripts.

### What the 0.76.0 round is evidence for

**Two tmbstan fits in one process share one objective.** tmbstan's
`sampling()` sets process-global pointers to the fit it runs, so
`log_prob()` of an earlier stanfit evaluates the latest one. The
consolidation compared two runs that way, found them bitwise equal at
800 points and called the cause of their different draws unknown; the
release review evaluated each right after its own run and found the
ulp. A comparison that cannot fail is not evidence: build the case
where the two must differ first.

**A lane's log can be older than its last fix round.** Three of the
lanes' logs that the equality rerun compares against (mocode's
`edge-lane.txt`, acprior's `agree.txt` and `agree-up.txt`) were written
before those lanes' fix rounds and differ from the merged build. Rerun
on each lane's own final library, the merged build's output is the
same line for line (`dev/rel076-equality-lanelib.sh`). Compare with the
lane's final build, not with the newest log of that name.

**An expected figure in a brief is arithmetic, not a measurement.**
The brief added the 13 bin-2 rows to 607 and did not take
off the bin-1 row the same round turns into a divergence. Regenerate,
then report what the ledger says.

**`cd dir && ( ... ) &` runs the cd in the background.** One
consolidation command started the bin-2 pipeline that way and then ran
`dev/ciharden-scan.sh` in the session's own directory, the main
checkout, which made an empty run there (two gitignored files, removed;
record, "The round's process"). Put the cd inside the background group
or use absolute paths.

## The 0.75.0 round, 2026-10-24

`dev/round-20261024.md` is this round's record. Each lane has its
findings (`dev/smallfix-findings.md`, `dev/bhaz-findings.md`,
`dev/brmsnames-1024-findings.md`, `dev/sarfcor-findings.md`) and its
review with the short-check sections (`dev/reviews/2026-10-24-<lane>.md`).
Four lanes off `18b536fc` (0.74.0 plus two CI commits), merged onto
`4fd134db` (18b536fc plus one CI commit to `pkgdown.yaml` and
`dev/release/build-docs.R`) in this order:

- `smallfix`: one refusal names every package a saved fit lacks
  (`stop_extension_parts()`), `keep_missing()` in every handler that
  could swallow it, `require_fit_packages()` at the head of methods; the
  `mixture_mvn()` test debt; the Dirichlet refusal at the model
  (BREAKING) and `print(curve, digits =)`; the CI skip gate for core.
- `bhaz`: `frmtmb_register_aterm(option =)`; brms's `bhaz()` in
  frmtmb.survival, and `cox(knots =, Boundary.knots =)`.
- `brmsnames`: brms's coefficient names where the parser's alias
  leaked (`.frm_expx` is `expx`), the prior table in brms's names and
  order, `fixef()` order of nonlinear and multivariate rows (BREAKING).
- `sarfcor`: `frmtmb_register_autocor()`, core's `fcor()`,
  frmtmb.spatial's `sar()`, `frmtmb_autocor_pars()` and the `bounds =`
  prior slot (part A); `frmtmb_covstruct_info()`, `full_rank` and the
  autoscale message (part B). The user kept both parts.

Versions: frmtmb **0.75.0** (minor, breaking); frmtmb.survival
**0.2.0**; frmtmb.spatial **0.2.0**; frmtmb.sample **0.22.0** (minor:
its NEWS has breaking changes and new features); frmtmb.latent
**0.7.1** (tests only); eam 0.13.0, coupling 0.6.1, ode 0.7.2 and learn
0.7.3 unchanged. All four changed extensions need frmtmb 0.75.0,
measured by `dev/rel075-floor.sh`. The release library is
`C:/Users/adf44/source/r/rellib-r13`; `rellib-r12` keeps 0.74.0.

### Decisions made at the consolidation, confirmed by the coordinator

The coordinator confirmed the decisions below, except one: the
consolidation kept smallfix's head check in
`recover_data.frmtmb_fit()`, and the release review showed that it
broke emmeans (next section).

- **One refusal across three seams.** sarfcor's residual-term check
  joins smallfix's collector (`seam_missing_parts()`), the "session"
  sentence names the seam that registers the part, and
  `ac_block_entry()` refuses through the same builder.
- **Two integration fixes the merged tree needed**, found by the
  consolidation's own test: `autocor_matrix()` and
  `frmtmb_autocor_pars()` check the whole fit first, and
  `frmtmb_covstruct_info()` uses `keep_missing()`.
- **The consolidation's test** of a saved fit with parts from three
  absent packages, and of a real saved `car()` and `sar()` fit. The
  brief's "cox + car + sar" cannot be one model: `sar()` takes a
  gaussian or student response, and `cox()` is refused in a
  multivariate fit.
- **frmtmb.sample is 0.22.0**, a minor bump, and **frmtmb.latent
  0.7.1** with the floor 0.75.0, because a test pins 0.75.0's message.
- **The upstream entries are renumbered**: brms-36 is bhaz's
  interval-censored `cox()`, brms-37 to brms-39 are sarfcor's; the
  references in sarfcor's findings, review and scripts follow.
- **The bin-2 `sar()` and `fcor()` rows are not ported** in this pass
  (record, "The ported brms suite").

### The release review and its fixes

`dev/reviews/2026-10-24-release.md`: FIX, no blocker. The coordinator
asked for every finding but 7 (a test-count nit, filed) as a release
fix (`dev/rel075-release-fix.R`, sourced by the resolve step):

1. (major) `emmeans()` and `ref_grid()` of a fit with a missing part
   stopped with emmeans's "Perhaps a 'data' or 'params' argument is
   needed": smallfix's head check in `recover_data.frmtmb_fit()` is
   raised inside emmeans's `try()`. The head check is removed;
   `emm_basis.frmtmb_fit()` refuses by name; both calls are in the
   saved-fit test.
2. `frmtmb_covstruct_info("sar")` names a residual correlation term.
3. One `sar()` parameter, one prior slot: no internal row in
   `default_prior()`, and a prior on it is refused by its class.
   `ar()`'s same double slot is filed.
4. frmtmb.survival's README and vignette show `bhaz()`.
5. The saved-fit step's comments name both fits.
6. `class = "lagsar"` on an error SAR names the model's `errorsar`.
7. `dev/rel075-floor.sh` installs only outside the tree and keeps
   Git's bash and PATH; rerun.

### What the 0.75.0 round is evidence for

**A head check inside someone else's `try()` is lost.** smallfix put
the fit-level check at the head of every method that could meet a
missing package, `recover_data()` among them; emmeans calls that one
inside `try()` and discards what it raises. The consolidation read the
code's own comment, kept the check, and recorded that it "refuses
first": it did, and nobody saw it. A check is placed where its caller
lets it through, and a test calls the method the way users do
(`emmeans::ref_grid()`), not the method directly.

**A safety check's refusal is reported, not worked around.** When the
tool refused a deletion written with shell variables, the
consolidation reran it with literal paths. The rule is to report the
denial and leave the deletion to the user; the record files it as an
incident.

**A merge can compose two lanes into a gap neither lane had.** Each
lane's saved-fit test passed in its own tree; only a test written for
the merged tree showed `autocor_matrix()` naming the `sar()` term alone
in a fit that also lacked a structure's package. Every integration of
two seams gets a test at the consolidation, written against the merged
tree and seen to fail without the resolution (`dev/rel075-seenfail.R`).

**A merge rerun restores generated files.** `dev/rel075-merge.sh`
restores `NAMESPACE` and `man/` to the HEAD blob, so an install without
roxygen after a rerun builds the old exports; the extensions then fail
to load ("'frmtmb_register_autocor' is not exported"). Rebuild the
release library in full after every merge rerun.

**A brief's example can be impossible.** "A saved fit with cox + car +
sar" cannot be fitted; measure that first, then test the nearest
construction and say why.

## The 0.74.0 round, 2026-10-23

`dev/round-20261023.md` is this round's record. Each lane has
`dev/<lane>-findings.md` and `dev/reviews/2026-10-23-<lane>.md`. It is
the first round under the core boundary rule (`dev/core-boundary.md`,
adopted the same day): every lane moved a feature between packages,
measured its guardrail-3 numbers first, and the user decided at the
consolidation to keep all five moves. Three lanes off `19b162b9`
(0.73.1), merged in this order:

- `splitsurv`: curve inference (`frm_curve()`, `frm_curve_deriv()`,
  `frm_curve_feature()`) from frmtmb.spline into core, with the family
  hook `post$curve_check`; the new extension frmtmb.survival with
  `cox()`, `cox_baseline()`, `royston_parmar()` and `rp_floored()`; the
  seams `frmtmb_family(simplex_pars =, aterm_data =)`;
  `R/extension-home.R`, so a saved fit names the extension its family
  needs; frmtmb.spline retired.
- `mvnlatent`: `mixture_mvn()` from core into frmtmb.latent; the
  exports `us_sigma()` and `us_theta()`; the mixture component
  interface documented in `?frmtmb-extension-api`.
- `spatial`: `frmtmb_register_covstruct()` and the covariance structure
  registry; `R/brms-vocabulary.R` and `test-guardrail1.R` (every brms
  name resolves in core, in a named extension, or is refused by name);
  the new extension frmtmb.spatial with `car()` and `spde()`;
  `frmtmb_local_registries()` exported.

Versions: frmtmb **0.74.0** (minor, breaking); frmtmb.survival
**0.1.0** and frmtmb.spatial **0.1.0** (new); frmtmb.latent **0.7.0**;
frmtmb.sample **0.21.1**; frmtmb.learn **0.7.3**; frmtmb.spline
retired; eam 0.13.0, coupling 0.6.1 and ode 0.7.2 unchanged. Every
floor was measured (`dev/rel074-floor.sh`): survival, spatial, latent
and sample need frmtmb 0.74.0; learn keeps 0.68.0. The release is in
lockstep: frmtmb.latent 0.6.2 does not load with frmtmb 0.74.0. The
release library is `C:/Users/adf44/source/r/rellib-r12`; `rellib-r11`
keeps 0.73.0.

**A release fix.** frmtmb.survival's first R CMD check --as-cran was
`Status: 1 ERROR, 1 NOTE`: `test-saved-fits.R` could not build its
"not installed" library stack under R CMD check. The coordinator
approved a release fix (`dev/rel074-release-fix.R`): one helper both
packages copy (`helper-absent-package.R`, frmtmb.spatial's two stacks
and probe), and an "installed but unloadable" arm for frmtmb.survival
that runs under R CMD check. The check is now `Status: 1 NOTE` (V8);
`dev/round-20261023.md`, "R CMD check --as-cran".

**The release review** (`dev/reviews/2026-10-23-release.md`, FIX, no
blocker) is applied as release fixes (`dev/rel074-release-review-fix.R`;
the record's "After the release review"): the README names
frmtmb.spatial and says `car()` and `spde()` are there; the README and
`vignette("brms-migration")` give one way to install every extension
(the user decided on 2026-10-23: no meta-package); the README's first
example runs; a `cox()` fit saved by frmtmb 0.73.0 is refused by name
with "Refit the model with library(frmtmb.survival) loaded"
(`old_core_family_home()`); NEWS tells users to remove frmtmb.spline;
the one-session runner watches the covariance-structure registry; a CI
job keeps the two helper copies identical; and the nits. After them:
441 of 441 files, 29752 assertions, 0 failures; frmtmb `Status: 1
NOTE`, frmtmb.sample `Status: OK`.

### Decisions made at the consolidation, accepted by the user

- **frmtmb.sample's floor rises to frmtmb 0.74.0** although its code
  loads on 0.73.0: four of its test files fail there (they read
  `frmtmb:::cs_entry()`, `frmtmb_family(simplex_pars =)` and
  frmtmb.latent 0.7.0's `mixture_mvn()`), and the rule is that a file
  that fails only on the old core names the floor.
- **The ported brms suite ran with every extension loaded**
  (`FRMTMB_PORT_EXTENSIONS=all`, a new switch of `dev/brmsport-run.R`
  and `dev/brmsport-bin2-run.R`), as guardrail 4 asks.
- **`dev/ciharden-scan.sh` has a plain tier** (`CIHARDEN_UNGATED=true`)
  and runs the nine suites.
- **The NEWS sentences over 25 words are split** (19), as at 0.73.0,
  with four edits of fact (`dev/rel074-news-split.R`).
- **`mixture_mvn()`'s row of `dev/core-boundary.md`** keeps the user's
  wording ("none; no string spelling reaches it").

### What the 0.74.0 round is evidence for

**A lane measures what it can load.** splitsurv compared the compat
table with core and the one extension it changed; the merged table with
every extension loaded also differs on 7 notes of frmtmb.learn, the
note splitsurv reworded there (`dev/rel074-compat.R`). A move's compat
check loads every extension, as `test-compat-extensions.R` does.

**Two absent-package checks compose, but one at a time.** A fit with a
family from frmtmb.survival and a structure from frmtmb.spatial is
refused by a package's name in every method that needs one, and never
fails on an internal (`dev/rel074-saved-mixed.R`); the message names
one missing package per call.

**R CMD check is a different library stack.** A test that builds its
own `.libPaths()` for a child process passed one file per process and
in one session, and failed only under R CMD check --as-cran, where the
tests see one temporary library holding the package and every
dependency. A lane that writes such a test runs R CMD check after it,
and probes the stack before it asserts, as frmtmb.spatial's does.

**A gate can hide a test from the skip count.** A file that defines its
gated tests only when the gate is set reports no skip in the plain
tier, so "rerun the files with a skip" missed 634 assertions; the
release reran every file with every gate on.

**A new package's first roxygen pass links to another package.** In a
fresh `man/`, `[cox()]` resolved to `brms::cox`, because the package's
own topic did not exist yet; the second pass links locally. Run
roxygen twice on a new package and compare.

## The 0.73.0 round, 2026-10-20 to 2026-10-21

`dev/round-20261021.md` is this round's record, with the verification
of the release tree under the reference BLAS and under OpenBLAS 0.3.26
and 0.3.32 (BLAS and LAPACK, 4 threads). Each lane has
`dev/<lane>-findings.md` and `dev/reviews/2026-10-20-<lane>.md`. Four
lanes off `4854e425` (0.72.0), merged in this order:

- `loodraws`: brms's `loo_R2()`, `loo_epred()`, `loo_linpred()`,
  `loo_predict()`, `loo_predictive_interval()` and
  `loo_model_weights()` on draws; brms's `r_eff` rule on a subset of
  draws; `posterior_average()` weights from `loo_model_weights()`; ten
  generics in core (`R/loo-generics.R`) with their refusals on a fit.
- `cdfs`: the CDFs of 23 more families (`R/families-cdf.R`), the
  log-scale window code of `row_lpdf()` and the `cdf_log` slot, brms's
  reading of a discrete censoring bound, brms's two `trunc()` checks,
  the exgaussian closed form, `simulate()`'s inverse-CDF fallback.
- `nb2dm`: `brmsfamily("negbinomial2")` and `dirichlet_multinomial()`,
  and the log-scale Dirichlet draw.
- `fuzzfix`: quadrature refuses a row in more than one block (the
  exact iterated route is `dev/fuzzfix-nested-quad.patch`), and the
  refusal of a weighted likelihood with no maximum
  (`check_weights_unbounded()`).

Versions: frmtmb **0.73.0** (minor, breaking); frmtmb.sample
**0.21.0** (minor), floor frmtmb 0.73.0 (measured: 0.21.0 does not
install over 0.72.0, `dev/rel073-floor.sh`). The six other extensions
keep their versions and floors: no lane changed a file of theirs. The
release library is `C:/Users/adf44/source/r/rellib-r11`; `rellib-r10`
keeps 0.72.0. In brms's ported suite bins 1 and 2 pass 607 of 823
(73.8%; 574 at 0.72.0), 607 of 706 (86.0%) without the decisions; bin
1 alone 392 of 494, unchanged.

### Decisions made at the consolidation, for the user to confirm

- **The emmeans leak is fixed in the generator**, with
  `withr::local_package("emmeans")` at the file's top level, which
  detaches at the end of the file (measured; `teardown_env()` would
  not).
- **rstan's Rcpp module classes are tagged as a load** in
  `dev/release/run-session.R` ("[with load of a Stan model (rstan)]"),
  measured to come from rstan alone; every other global stays
  untagged.
- **`loo_predict` joins core's `_pkgdown.yml` index**, which lane
  loodraws had left out (`pkgdown::check_pkgdown()` failed on it).
- **The NEWS sentences over 25 words are split** (28), as at 0.72.0.
- **`?posterior_average` links its weights to `sample-loo_model_weights`**
  by name: roxygen had resolved the bare `[loo_model_weights()]` to
  `loo::` or to the local re-export, depending on the library it ran
  against.

The release review (`dev/reviews/2026-10-21-release.md`, SHIP AFTER
FIXES) is applied ("After the release review" in the record): the
weights and asym_laplace check now runs on the uncensored rows of a
`cens()` model and carries `zero_inflated_asym_laplace()`; the four
core redirects of `kfold`, `reloo`, `loo_moment_match` and
`loo_subsample` point at `reference/loo_predict.html`; lane adhess's
filed poisson `trunc(lb = 1)` seed-19 item is closed by lane cdfs's
log-survivor normalizer (measured by the review).

### What the 0.73.0 round is evidence for

**A lane that adds an Rd file to another package must add its site
index entry.** Lane loodraws moved ten generics into core with a new
`man/loo_predict.Rd`; its core suite and R CMD check passed, and only
`pkgdown::check_pkgdown()` (the site workflow's own check) stopped on
the missing topic. A lane that adds or renames a topic runs
`check_pkgdown()` for the package it touched.

**A `teardown_env()` cleanup is a run's, not a file's.** The obvious
fix for an attached package, `withr::local_package(..., .local_envir =
teardown_env())`, keeps it attached for every later file of a session;
only the file's own environment detaches it at the file's end
(`dev/rel073-local-package.R`). The one-session gate is the instrument
that shows the difference.

**A three-way numbering clash is cheap only if one pass renumbers.**
Three lanes each appended a "brms-30". Renumbering one map at a time
(30 to 31, then 31 to 32) would have shifted cdfs's entries twice;
`dev/rel073-resolve.R` maps each text once, and two passages of a
review that discuss the numbering itself keep their numbers with a
note.

## The 0.72.0 round, 2026-10-14 to 2026-10-20

`dev/round-20261020.md` is this round's record, with the verification
of the release tree under the reference BLAS and under OpenBLAS 0.3.26
and 0.3.32 (BLAS and LAPACK, 4 threads). Each lane has
`dev/<lane>-findings.md` and `dev/reviews/*-<lane>.md`. Nine lanes off
`47cf659d` (0.71.0 plus the test-order fix), merged in this order:

- `testiso`: test isolation (`local_frmtmb_registries()`, the
  compatibility contract with all seven extensions loaded), the
  one-session release gate (`dev/release/run-sessions.sh`), and the
  coverage workflow (core alone, plus one job per extension). No
  package code.
- `catmulti`: brms's `categorical()` and `multinomial()` (`refcat`,
  `refcat = NA` with the reference model's degrees of freedom, the bare
  `multinomial()`), `fitted()` and `posterior_linpred()` of every
  location predictor, `logistic_normal()` in `mvbf()`, the cumulative
  tied-threshold `NaN`.
- `prilabel`: per-coefficient `sd` priors, brms's labels and order of
  nested and multi-membership groups (and of the theta parameters), one
  Student-t latent scale per factor.
- `sampleprior`: brms's default priors in frmtmb.sample, `seed =`,
  `sd_`/`cor_` draw columns; core `uniform()`, `lognormal()`, support
  guards on the tape, the logit Jacobian.
- `bin2port`: bin 2 of brms's suite on 0.71.0's families, with value
  checks (dev only).
- `ndaterm`: every addition term on `newdata`, coded as the fit coded
  it; mixture and `hmm()` term rows; the double softmax of
  `predict(newdata =)`.
- `illse2`: Cook's distance, D1 and D2 through the covariance factor;
  REML and nonlinear bodies in the design's coordinates; bounded
  Newton steps; the 400-parameter cap.
- `adhess`: the exact AD Hessian for every fit without random effects;
  RMS-scaled steps with them; "nocurv"; frmtmb.coupling's refusal.
- `drop1raw`: `drop1()`, `step()`, `update()` and `influence()` refit on
  the data the fit holds, with the fit's bases.

Versions: frmtmb **0.72.0** (minor, breaking); frmtmb.sample **0.20.0**
(minor), floor frmtmb 0.72.0; frmtmb.eam **0.13.0** (minor), floor
frmtmb 0.72.0 and Suggests frmtmb.sample (>= 0.20.0); frmtmb.latent
**0.6.2** (patch), floor frmtmb 0.72.0; frmtmb.coupling **0.6.1**
(patch), floor kept at frmtmb 0.61.0 (measured: its suite passes on
0.71.0; update the two together). Each floor was measured by
installing the extension alone over rellib-r9 (`dev/rel072-floor.sh`).
frmtmb.spline, frmtmb.learn and frmtmb.ode keep their versions. The
release library is `C:/Users/adf44/source/r/rellib-r10`; `rellib-r9`
keeps 0.71.0. In brms's ported suite bins 1 and 2 pass 574 of 823
(69.7%; 549 at 0.71.0), 574 of 706 (81.3%) without the decisions; bin
1 alone 392 of 494 (391).

### Decisions made at the consolidation, for the user to confirm

- **frmtmb.coupling keeps its floor** (frmtmb >= 0.61.0) where the
  brief and lane adhess named 0.72.0: installed alone over 0.71.0, its
  11 test files pass. A floor cannot stop coupling 0.6.0 meeting core
  0.72.0, so its NEWS lead says to update the two together.
- **frmtmb.eam suggests frmtmb.sample (>= 0.20.0)**: lane sampleprior's
  `test-sample-ndt-default.R` in frmtmb.eam fails with frmtmb.sample
  0.19.0's defaults.
- **The two consolidation fixes** (`dev/round-20261020.md`, "Found at
  the consolidation"): `confint_varcorr()` no longer spreads a
  "nocurv" NaN to rows that do not move along it, and lane
  sampleprior's draw-name test takes prilabel's `variables()` order.

### What the 0.72.0 round is evidence for

**A lane's reruns after a fix round are a subset, and the subset can
miss its own regression.** Lane adhess's full suites passed before its
fix round 1; the fix (`sdr_nan_nonpos()`, a NaN row for a non-positive
variance) was rerun on the files the lane chose, and
`test-review-v29.R` was not among them. On the merged build it failed
10 expectations, and on lane adhess's build alone too
(`dev/round-20261020.md`, "Found at the consolidation"). A fix round
that changes package code reruns the whole suite, as the round's
first build did.

**A test that writes another lane's order down is a merge conflict
git cannot see.** Lane sampleprior's new test pinned `variables()` in
0.71.0's order as "brms's names"; lane prilabel changed that order to
brms's. Each lane's suite passed; the merged one did not. When a test
asserts an order or a label, run it against the brms call it claims
to mirror rather than against the build of the day.

**Measure a floor; do not take it from the lane.** Lane adhess said
frmtmb.coupling needs frmtmb 0.72.0; installed alone over 0.71.0, its
whole suite passes, so its floor stays (`dev/rel072-floor.sh`). The
other three extensions whose floors rose fail on 0.71.0, as their
lanes said.

## The 0.71.0 round, 2026-10-13 to 2026-10-14

`dev/round-20261014.md` is this round's record, with the verification
of the release tree under the reference BLAS and under OpenBLAS 0.3.26
and 0.3.32 (BLAS and LAPACK, 4 threads). Each lane has
`dev/<lane>-findings.md` and `dev/reviews/2026-10-13-<lane>.md`. Seven
lanes off `d2ef6312` (0.70.0), merged on top of `253b3ef7` (0.70.1,
lane ci701), in this order:

- `dblbar`: every `||` term as brms reads it, `(x | gr(g, cor =
  FALSE))`, one block coded with the formula's contrasts;
  `gr(cor = )`, `mm(cor = )`; refusals of multi-column structures over
  `gr(cov/prec)`.
- `famsimplex`: `dirichlet()`, `logistic_normal()`,
  `brmsfamily("dirichlet2")`; a formula for the first location; a
  constant on a location refused (core and the extension families).
- `famuni`: `frechet()`, `gen_extreme_value()`, `discrete_weibull()`;
  frmtmb.eam's `wiener()` links in brms's order; the GEV `tmp_xi`
  chart for frmtmb.sample; `nan_unless_pos()` (`dev/rtmb-pitfalls.md`
  item 22).
- `postfix`: `frm_multiple()` methods by Rubin's rules,
  `pp_check(resp = )`, draws' `loo()` and `waic()` in brms's order,
  `conditional_effects(method = "predict")` on draws, the nonlinear
  default-start error (the nonlinear scan was removed).
- `optstart`: least-squares starts for location-scale fits,
  `restart_result()`, `lp_design_kappa()`, `fd_hessian()`'s retake.
- `illse`: accurate SEs on ill-conditioned designs through a factor of
  the covariance; the SE check on fits with a convergence warning; the
  flat simplex direction.
- `remeasure`: the vignette port and bin 2 of brms's suite (dev only).

Versions: frmtmb **0.71.0** (minor, breaking); frmtmb.sample **0.19.0**
(minor), floor frmtmb 0.71.0 (the new sampling-API exports and the GEV
chart); frmtmb.eam **0.12.0** (minor, breaking `wiener()` argument
order), floor frmtmb 0.70.0 (measured: its suite passes on 0.70.0);
frmtmb.latent **0.6.1** (patch, NEWS-only Breaking entry); frmtmb.spline
**0.10.2** (patch, no new floor: the new version works on 0.70.0 and
0.10.1 refuses loudly on the new core; update them together). The
other three extensions keep their versions. The release library is
`C:/Users/adf44/source/r/rellib-r9`; `rellib-r8` keeps 0.70.0. In
brms's ported suite bins 1 and 2 pass 549 of 823 (66.7%, 534 at
0.70.0), 549 of 709 (77.4%) without the decisions; bin 1 alone 391 of
494 (389).

### Decisions the user made in this round

- **`com_poisson` keeps frmtmb's mean parameterization.** brms 2.23.0
  fits `mu` as the COM-Poisson rate and post-processes another pmf
  (brms-24 in `dev/upstream-bugs.md`); `?brmsfamily` and
  `?frmtmb-families` say that frmtmb's `mu` is the mean.
- **Upstream bugs are written into `dev/upstream-bugs.md` only**
  (brms-24 to brms-27 this round); the user files them.

### What the 0.71.0 round is evidence for

**A trial merge with one lane is not a trial merge with the release.**
Lane illse's review ran a trial merge with lane optstart and found
`test-se-conv.R` failing; illse fixed that test. Its next test, added
in its fix round 2, had the same exposure (a premise that depends on
`frm()`'s default start) and failed 5 expectations on the seven-lane
tree. And illse's tree was on 0.70.0, so no trial merge of it had met
lane ci701's commit, with which `R/se-check.R` conflicted in four
blocks. Both were resolved at the consolidation
(`dev/round-20261014.md`, "Found at the consolidation"). A test that
pins a stop short of the maximum should give every parameter its
start, and the start must leave only the direction under test short.
Check that by removing the fix: the consolidation's first rewrite of
illse's test started every parameter at 0, so the intercept was short
too, and the test passed with illse's fix stubbed out; the release
review found it (`dev/reviews/2026-10-14-release.md`, 6a), and the
test now starts at `lm()`'s maximum with `iter.max = 0` and fails
without the fix (`dev/rel071-seconv3.R`). A test that fails on the
old build is not yet shown to pin the change you mean.
Attribute every difference a battery shows on the merged build
by running the battery on each lane's build too: here every verdict
difference of illse's false-alarm battery from 0.70.0 was lane
optstart's (the merged build equals optstart alone on 820 of 820 fits
under both BLAS), a battery optstart's own review never ran.

## The 0.70.0 round, 2026-10-07 to 2026-10-13

`dev/round-20261013.md` is this round's record, with the verification
of the release tree under the reference BLAS and under OpenBLAS 0.3.26
(BLAS and LAPACK, 4 threads, the ubuntu-24.04 runners' configuration)
and 0.3.32. Each lane has `dev/<lane>-findings.md` and
`dev/reviews/2026-10-08-<lane>.md`. Three lanes off `48c1e6c9`
(0.69.0), merged on top of `2884ccbd` (0.69.1), in this order:

- `moparity`: brms features that were refused or diverged: `mo(x) * f`
  and `(mo(x) | g)`, Dirichlet priors on `simo` and `sbhaz`,
  `update(prior =)` merging by slot, the `cox()` baseline in the
  `mo()` chart, floors for the dense blocks with several scales,
  `parnames()`, and frmtmb.sample's `cox()` baseline and non-centered
  exact `gp()`.
- `optrap`: fits that stopped short of the maximum: a unit for the
  identity-link `mu` of ten location-scale families, the restart that
  raises, `obj_settle()` on every discarded run, `newton_polish()`
  for fits without random effects, the `mo()` search.
- `seflags`: which SEs a fit reports: an exact design alias, the
  prediction test in the unit-diagonal metric, the sd boundary verdict,
  curved ridges, the degree-7 polynomial, an ill-conditioned design
  named as one (scoped down in punch round 2c: such a direction keeps
  0.69.0's NaN).

Versions: frmtmb **0.70.0** (minor, breaking); frmtmb.sample
**0.18.0** (minor, breaking), floor frmtmb 0.70.0 (it reads
`simplex_frame_terms()`); frmtmb.eam **0.11.3** (patch), floor frmtmb
0.70.0 (`test-defects.R` requires the prediction warning only 0.70.0
gives; rstan added to Suggests). The other five extensions keep their
versions and floors: no lane changed a file of theirs. The release
library is `C:/Users/adf44/source/r/rellib-r8`; `rellib-r7` keeps
0.69.1. In brms's own ported suite, bin 1 passes 389 of 494 (390 at
0.69.0): `brmsfit-methods:326` held at 0.69.0 by luck.

### Decisions the user made in this round

- **`mo(x) * f` and `(mo(x) | g)` follow brms 2.22.0's
  parameterization**, because brms 2.23.0 writes wrong Stan code for
  these terms (every column multiplied by the first covariate) and
  stops on a factor of three or more levels crossed with the main term
  (brms#1828). frmtmb fits the model brms 2.22.0 fits.
- **A Dirichlet prior on a MAP fit is a density on the simplex's
  unconstrained coordinates with the chart's Jacobian**, frmtmb's
  convention for every constrained parameter: the mode under
  `dirichlet(a)` is the mode of `dirichlet(a + 1/2)` on the weights,
  so `dirichlet(1)` is not the absence of a prior, and a concentration
  below 1/2 is refused (`?frm`, "Monotonic effects").
- **`update(prior =)` merges the new prior into the stored one by
  slot, as brms does**, and a row of a prior table with an empty prior
  returns its slot to the default, as in brms; `prior = NULL` or a
  zero-row prior keeps the stored one. A new specification is still
  checked as `frm()` checks it, where brms drops it silently.
- **Response units and `newton_polish()`**: the identity-link `mu` of
  ten location-scale families gets a unit (the MAD of the
  least-squares residuals, rounded to a power of two), and a fit
  without random effects takes Newton steps on the exact Hessian while
  it is positive definite and a step predicts a gain above 1e-6; a
  step is kept only if it gains at least 0.1 of its prediction, the
  steps are undone if they end where the Hessian does not factor, and
  a fit with random effects takes them too when its fixed design is
  near collinear (`predictor_design_kappa() > 1e5`). The last three
  are lane optrap's consolidation fixes of two defects the merged
  build showed (below).
  BREAKING: starting values and the optimizer's units change for those
  families.
- **An ill-conditioned design keeps 0.69.0's NaN**, with a message
  that names the design as ill-conditioned (and, where it applies,
  that the optimizer stopped short, with the Newton gain). Lane
  seflags' attempt to report accurate SEs there (punch round 2b) gave
  wrong finite numbers through emmeans, marginaleffects, `vcov()` read
  whole and the ordinal category route; "Accurate SEs on
  ill-conditioned designs" is its own lane (`dev/test-backlog.md`,
  Open, high). The coordinator chose the scope-down on the review's
  recommendation; the user did not pick it explicitly.
- **0.69.1 was committed straight to main** (`2884ccbd`, lane ci691,
  without a separate review round), and this round was consolidated on
  top of it as 0.70.0.

### What the 0.70.0 round is evidence for

**Each lane was right, and two of their unions were not, until lane
optrap fixed both before the release.** Every
reviewer's own reruns reproduce on the merged build, row for row (the
moparity attacks, the `cox()` study bit for bit, optrap's scale,
`mo()`, high-R^2 and wide studies, seflags' routes and false-alarm
battery). But the reruns on the merged build found defects that no
lane's own runs could show: lane optrap's unit makes raw polynomials
with `(1 | g)` stop short in silence at degree 6, and lane seflags'
correct degree-7 SEs removed the false "flat" warning that had covered
the same shortfall at degree 7; and lane optrap's `newton_polish()`
accepts a step whose actual gain is a tiny fraction of the predicted
one, which on a near-Poisson negative binomial runs the shape toward
its infinite MLE and leaves every SE NaN behind a large-gradient
warning. Lane seflags' false-alarm battery, run against the merged
build, found the second; the optrap review's polynomial probe, rerun
with an exact profile maximum, found the first. Lane optrap fixed both
in its consolidation fixes (a 0.1 gain gate with an undo, and Newton
steps on near-collinear random-effect fits), and the release tree took
them by a rerun of `dev/rel070-merge.sh`, then was verified again
(`dev/round-20261013.md`, "Found at the consolidation, fixed"). Run
each review's batteries on the merged build, not only the suites, and
keep the merge a script, so a lane's fix after the merge is one rerun.

## The 0.69.0 round, 2026-10-06 to 2026-10-07

`dev/round-20261007.md` is this round's record, with the verification
of the release tree under the reference BLAS and under OpenBLAS 0.3.32
(BLAS and LAPACK, 4 threads), the Ubuntu 26.04 runners' version; the
consolidation review (`dev/reviews/2026-10-07-release.md`) ran it under
0.3.26 as well, the 24.04 runners'. Each lane has
`dev/<lane>-findings.md` and `dev/reviews/2026-10-07-<lane>.md`.
Four lanes merged, in this order:

- `ciharden`: the test suites no longer depend on platform rounding (a
  scan of all eight suites under five BLAS and thread configurations,
  `dev/ciharden-scan.sh`, with an OpenBLAS build of R on Windows,
  `dev/ciharden-openblas.sh`); `test-perf.R` counts tape nodes and
  bytes, not time; the ported suite is seeded per block; two
  frmtmb.sample defects the runner showed; CI pinned.
- `surface`: the brms post-processing calls a ported script meets on
  an ML fit answer or refuse by name (`stancode()`, `standata()`,
  `pp_mixture()`, `add_criterion()`, `plot()`, `update()` with a stale
  prior, `frm_multiple()` pooled by Rubin's rules).
- `setier`: which standard errors a fit can report (tier 3 asked
  first, a curvature probe), boundary fits, separation, the nonlinear
  ridge with random effects.
- `optima`: fits that stopped short of the maximum (`mo()`, the probit
  tails, nlminb's rejected trial point, `cs()` ordinal mixtures).

Versions: frmtmb **0.69.0** (minor, breaking); frmtmb.sample
**0.17.0** (minor, breaking), floor frmtmb 0.69.0 (it imports the four
generics and the `mo_` helpers). The other six extensions keep their
versions and floors: no `R/` change, and the two test edits (frmtmb.eam
`test-sampling.R`, frmtmb.learn `test-families.R`) allow a warning
without requiring it and pass on 0.68.1. The release library is
`C:/Users/adf44/source/r/rellib-r7`; `rellib-r6` keeps the 0.68.1
build. In brms's own ported suite, bin 1 passes 390 of 494 (389 at
0.68.1).

### Decisions the user made in this round

- **Unseen `gp()` positions are grouped at 15 significant digits, as
  brms 2.23.0's new-data rule groups them** (`pos_rowkey()`, brms's
  paste key). Lane ciharden's exact key was reverted in its punch
  round (review B2); the 0.68.1 backlog item asking for an exact key
  is closed by this decision.
- **frmtmb.sample puts brms's default `dirichlet(1)` on every `mo()`
  simplex and reports brms's `simo_` weights. BREAKING**: the sampler
  put no density on a simplex, so the fit's coordinates were its
  prior; `prior = "flat"` keeps `dirichlet(1)`, the flat density on the
  simplex.
- **A boundary fit is reported as a message, as lme4 reports a
  singular fit**: "Boundary (singular) fit: ...", class
  `frmtmb_boundary_fit`, not the standard-error warning; `check_se =
  "stop"` makes it an error and `"ignore"` silences it. When nlminb
  stops with a non-zero code there, the gradient is within `grad_tol`
  and every lost parameter is at its edge, the message replaces the
  convergence warning, as lme4's checkConv() does.
- **The quadrature carve-out**: under `quadrature = TRUE` no boundary
  verdict is given, because the upward check that tells a boundary
  from a stop short of the maximum cannot run on a quadrature
  objective; such an sd keeps the "flat" standard-error warning, as
  at 0.68.1 (setier, punch 2b).
- **The Stan job runs on Ubuntu 26.04 in `ubuntu-next.yaml` only**
  (weekly and on demand, beside 24.04, its own cache key per image),
  never on push; `brms-likelihood.yaml` is unchanged.
- **CI is pinned to `ubuntu-24.04`**, every workflow with the dated
  comment "Pinned 2026-10-06". When `ubuntu-next.yaml` is green on
  26.04, move the pins (search the workflows for that comment).
  `ubuntu-latest` moves to 26.04 between 2026-10-19 and 2026-11-19.

### What the 0.69.0 round is evidence for

**A clean three-way merge still needs the site checked.** No lane
touched `_pkgdown.yml`, every lane's own check passed, and the merged
tree failed `pkgdown::check_pkgdown()`: lane surface's four new core
topics were in no index, and two `redirects:` entries for retired core
pages (`reference/pp_mixture.html`, `reference/stancode.html`) would
have overwritten the live pages surface brought back. Run
`check_pkgdown()` in every lane that adds or moves a topic.

**Two lanes can each be right and their union platform-dependent.**
On the merged tree, two results fail when OpenBLAS (0.3.26 or 0.3.32)
also provides LAPACK, the ubuntu-24.04 and 26.04 runners'
configuration, and pass with either version as BLAS only: a fuzz spec
errors where base fits it (lane optima's build errors too), and the
ported row `data-helpers:7` gets a false alarm on five estimable
predictions, because optima's change makes fixture 1 converge there
and the 0.68.0 estimability test then reads a lost basis with the
intercept's share zeroed (the consolidation review, B1). Neither file
runs on CI; both are filed. The lesson: a lane's OpenBLAS pass must
route LAPACK too, `dev/ciharden-openblas.sh <ver> lapack`, because
`dev/optima-openblas.sh` and `dev/cifix-openblas.sh` route BLAS only,
and lane optima's pass with the first could not show either result.
Run each lane, and the merged tree, under the runners' configuration
(`dev/ciharden-scan.sh`), not only the reference.

## 0.68.1, the CI fix, 2026-10-06

0.68.1 (with frmtmb.spline 0.10.1, floor `frmtmb (>= 0.68.1)`) fixes
the Ubuntu failures of 0.68.0: check-frmtmb.spline
(`test-difference.R:354` and `:490`), R-CMD-check on both Ubuntu jobs
and test-coverage (`test-se-check.R:383`). Lane cifix,
`dev/cifix-findings.md`; review `dev/reviews/2026-10-06-cifix.md`
(MERGEABLE). What it changed:

- core: a fit with random effects and a lost standard error gets a
  positive semi-definite joint covariance (`joint_cov_repair()`), so
  its predictions no longer read a negative coefficient variance, and
  a prediction along a lost direction gets NaN and one warning, as
  without random effects (BREAKING for those fits only);
- core: each distinct unseen `gp()` position is kriged once, so rows
  at one position are bit-identical on every BLAS;
- spline: one standard-error arithmetic in core's clamp order, one
  source for the rows without a standard error (`nonest` as well as
  `se_nonest`), and a row both grids of a difference share evaluated
  once;
- tests: the separation fit gets a budget that converges on both BLAS,
  and two fits whose sd runs to zero allow the platform-dependent SE
  warning.

`dev/cifix-openblas.sh` builds a copy of R with OpenBLAS 0.3.26, the
runners' version, which reproduced every CI symptom on this machine.
The one unverified step is the real Ubuntu runner: only a push shows
it green. Watch `test-cumulative-cs.R:132` and
`test-ordinal-mixture.R:751` when `ubuntu-latest` moves to Ubuntu 26
from 2026-10-19. Found-not-fixed items are in `dev/test-backlog.md`,
"Filed at 0.68.1"; the largest is the SE check's tier 1 accepting a
pure-noise Hessian row (117 fits in 58 test files), for a planned
round. On this hybrid CPU, compare results within one process or pin
the core type (`dev/rtmb-pitfalls.md` item 21).

## The 0.68.0 round, 2026-10-05 to 2026-10-06

**The user decided on 2026-10-06 to hold the release for lane nanse**,
which was merged onto the four-lane tree after its review
(`dev/round-20261005.md`, "Lane nanse merged"). Its check runs by
default: every `frm()` call builds the Hessian at the optimum and
warns when a parameter has no standard error, naming it and why
(`frmtmb_control(check_se = "ignore")` skips it, `"stop"` makes it an
error). The cost is 8 to 24 percent of the fit on models without
random effects (the exact Hessian), 19 to 23 percent on mixed models
with up to about 23 outer parameters (a finite-difference Hessian),
and 3 percent or less on larger mixed models, where the check waits
for the first use of a standard error (lane nanse's findings,
"Cost").

`dev/round-20261005.md` is this round's record, with the verification
of the release tree. Each lane has `dev/<lane>-findings.md` and
`dev/reviews/2026-10-05-<lane>.md` (`2026-10-06-nanse.md` for nanse).
Five lanes merged:

- `fixes`: `emmeans()` on transformed predictors, `trunc()` and `se()`
  in `conditional_effects()`, `posterior_linpred(incl_thres = TRUE)`,
  per-threshold `Intercept` rows, `confint()` names for the ordinal
  internal parameters, compat rows for `disc` and the threshold
  structures, brms's `s()` null-space basis (BREAKING), and a warning
  at the optimum when a nonlinear model's likelihood is flat.
- `gpby`: `gp(x, by = )` with brms's `cmc`, `gr`, `iso` (default
  `TRUE`, BREAKING) and prior classes `sdgp` and `lscale`; the exact
  `gp()`'s kriging covariance at new positions
  (`frm_lp_basis(extra_cov = TRUE)`, `frm_extra_cov_deriv()`), its draw
  in frmtmb.sample, and frmtmb.spline's curves over it.
- `ordmix`: ordinal mixtures (`order = "none"` and `"mu"`),
  `hurdle_cumulative()` with `thres(gr = )` and `cs()`, the
  degenerate-component and flat-direction warnings.
- `vigport`: the brms vignette port measured again (records only); its
  defects are filed.

- `nanse`: a fit says which parameters have no standard error and
  why, keeps the others' (sdreport()'s inverse lost all of them), gives
  NaN predictions along a lost direction, and `frm_allfit()` starts
  where `lme4::allFit()` does (vigport's defects 8 and 9, lane fixes'
  consolidation items 1 and 3).

Versions: frmtmb **0.68.0**; frmtmb.sample **0.16.0**, floor frmtmb
0.68.0 (`ord_thres_linpred()`, `gp_brms_natural()`, the `sdgp` and
`lscale` classes, the kriging draw, ordinal mixtures);
frmtmb.spline **0.10.0**, floor 0.68.0 (`frm_lp_basis(extra_cov =
TRUE)`, `frm_extra_cov_deriv()`); frmtmb.ode **0.7.2**, floor 0.68.0,
no code change (its test requires the flat warning); frmtmb.learn
**0.7.2**, floor 0.68.0, no code change (its test requires nanse's SE
warning). The release
library is `C:/Users/adf44/source/r/rellib-r6`; `rellib-r5` keeps the
0.67.0 build. In brms's own ported suite, bin 1 passes 389 of 494
(387 at 0.67.0).

### Decisions the user made on 2026-10-06, for this round

- **The exact `gp()` nugget stays `1e-6`.** Over 180 fits (the gpby
  re-check, r1) `1e-6` converges on every one, `1e-8` warns on 6 of 45
  and `1e-12` converges on none. Per-draw parity with brms's `1e-12`
  jitter, if wanted later, belongs to a sampled fit, not to maximum
  likelihood.
- **`posterior_linpred(incl_thres = TRUE)` on `hurdle_cumulative()`
  stays refused**, a divergence: brms returns the hurdle probability
  beside the threshold predictors times `1 - hu`, the predictor of
  nothing (`dev/upstream-bugs.md`, brms-21).
- **`cs()` on `cumulative()` is lifted**, as brms 2.23.0 fits it, with
  brms's warning; `hurdle_cumulative()` and cumulative mixture
  components follow the same rule. A row whose offsets cross two
  thresholds keeps frmtmb's `NaN` `fitted()` and `NA` draw, where brms
  returns a negative probability (brms-22).
- **`conditional_effects(method = "predict")`'s `estimate__` is the
  median** of the plug-in predictive draws, the same draws its band is
  a quantile band of, as brms reports its median. BREAKING.
- **`predict()` on an ordinal `cs()` row whose thresholds cross**:
  the row is `NA` only where the plug-in estimate crosses, as
  `fitted()` is `NaN` there; elsewhere the proportions are over the
  replicates that do not cross, and the call warns once with the
  replicates dropped per row (the release review's m3).
- **vigport's defects 8 and 9** (silent NaN standard errors on
  `mo(income) * age`; `frm_allfit()` refits a nonlinear model with
  `start = NULL`) are addressed in lane nanse, with fixes' consolidation
  items 1 and 3.

### What the 0.68.0 round is evidence for

**A test runner can answer for a missing variable.** Lane gpby's new
`test-gp-by.R` asks for a prediction without the `by = f` column and
expects frmtmb's refusal. Under `dev/release/run-tests.R`, whose global
`f` held the test file's path, R's lookup rule found that `f` and the
test saw "New factor levels" instead. R CMD check has no such global.
The runner now runs inside `local()`; a runner must leave the global
environment empty.

**Four clean lanes still need their interplay written down.** The
three `R/priors.R` conflicts were the visible part; the compat cells,
the `incl_thres` refusal and the mixture labels merged clean and were
each wrong on the merged tree until resolved. Each resolution has a
test seen to fail on a build with that one resolution undone
(`dev/rel068-mutants.R`).

## The 0.67.0 round, 2026-09-30

`dev/round-20260930.md` is this round's record, with the verification
of the release tree. Each lane has `dev/<lane>-findings.md` and
`dev/reviews/2026-09-30-<lane>.md`. Three lanes of brms parity merged:

- `ordinal`: brms's `disc` on every ordinal family,
  `threshold = "equidistant"` and `"sum_to_zero"`, and every brms link
  on `acat()`, with brms's draw names in frmtmb.sample.
- `ceplot`: `plot()` of conditional effects and of a hypothesis
  returns plot objects and takes brms's arguments; brms's effect check;
  new levels on crossed and `mm(by = )` terms in
  `conditional_effects()`; `fitted(sample_new_levels = "old_levels")`
  with one choice per grouping factor and per call; `parnames()` on a
  fit; `nsamples(incl_warmup = TRUE)`; brms's order in
  `posterior_samples(pars = )`.
- `formrobust`: addition-term expressions, `update()` that keeps the
  parameter formulas, bernoulli on any two values, `ar()` and the other
  autocorrelation functions, `autocor()`, `drop_unused_levels`, offsets
  in the grids of `conditional_effects()` and `emmeans()`, and the
  ARMA fill of a missing response on newdata.

Versions: frmtmb **0.67.0**; frmtmb.sample **0.15.0**, floor frmtmb
0.67.0 (it imports `parnames()` and calls `ord_delta_info()`,
`arma_cond_fill_dpars()`, `arma_cond_fill_epred()` and
`response_codes_newdata()`). The other six extensions are unchanged.
The release library is `C:/Users/adf44/source/r/rellib-r5`; `rellib-r4`
keeps the 0.66.0 build. In brms's own ported suite, bin 1 passes
387 of 494 (368 at 0.66.0).

### Decisions the user made on 2026-09-30, for this round

- **The ordinal default of `conditional_effects()` stays the
  per-category display**, a deliberate divergence from brms. brms's
  default is the expected category number, and its own warning calls
  that display likely invalid for ordinal families and asks for
  `categorical = TRUE`, which is frmtmb's default. brms's warning comes
  with `categorical = FALSE`. Ledger row `brmsfit-methods:217` is a
  divergence.
- **A family-default fixed `disc` is hidden, as brms hides it.** A
  `disc` held at 1 shows nowhere: not on the `Links:` line, not as a
  `Fixed dpar:`, not in `fixef(flatten = TRUE)` or `coef()`. A `disc`
  the user fixes at another value still shows.
- **Addition terms read data variables only**, as brms requires. A
  variable found only in the formula environment is refused by name,
  because predictions on newdata and refits read the term again, where
  an outside value can be missing or changed. Predictor formulas and
  offsets keep R's rule of looking in the formula environment.
- **Filed as an idea, not built: snapshot outside constants at fit
  time**, so that a predictor or offset that reads a value outside the
  data reads the same value in a later prediction or refit.
- **bernoulli codes any two values by level order, as brms does**, and
  warns when both values lie strictly between 0 and 1, since that looks
  like a proportion (lme4#682). Effect coding such as -0.5 and 0.5 does
  not warn.
- **The ARMA fill of a missing response**: frmtmb.sample's draws sample
  the fill, as brms does, and the maximum likelihood `fitted()` fills
  with the expected value, a deliberate divergence, since a fit has no
  draws to carry the fill through.

### What the 0.67.0 round is evidence for

**A worker's report can fail to arrive; read its findings file.** The
findings file is the record the lane is judged on. When a report does
not reach the coordinator, the lane's state is still in
`dev/<lane>-findings.md`, and the consolidation read each lane's
intent there: for every conflict, and for the ledger rows each lane
said it flipped.

**A small punch round can reverse the direction of a change.** Lane
formrobust's first `emmeans()` offset change went one way, its review
found brms does the other, and the punch round turned it around. Run
brms on the case, rather than only reading its code: the reversal came
from a run (`dev/formrobust-rev-log/brms-emm.txt`), not from the
source.

**A merge of two clean lanes can make a wrong answer that neither
lane has.** formrobust put offset variables in the model frame for the
grids; ceplot's effect check read every variable of the terms. Together
they accepted `effects = "time"` for a variable only an offset reads,
which brms refuses, while each lane alone stopped or answered for
another reason. Run the interplay of two lanes that touch one function
on the merged tree, with brms beside it (`dev/rel067-effoffset.R`).

**A lane's record of a ledger row can be wrong in the direction it
expects.** formrobust recorded `brmsfit-methods:955` as a class-name
divergence because `update()` "refits". On the lane's own library the
refit stops at its starting values. Read the message the merged build
records for every row whose message moved
(`dev/rel067-msgdiff.R`), not only the rows a builder stops on.

## The 0.66.0 round, 2026-09-29

`dev/round-20260929b.md` is this round's record, with the verification
of the release tree. Each lane has `dev/<lane>-findings.md` and
`dev/reviews/2026-09-29-<lane>.md`. Six lanes of brms parity merged:

- `sampfix`: frmtmb.sample reads `frm_sample(laplace = TRUE)` draws in
  their own layout (they returned NaN or finite wrong numbers) and
  refuses by name what reads the integrated values; `pp_check()`'s
  `loo_*` types; `stanfit = NULL` draws; brms's names for ordinal
  draws; one-parameter models.
- `fams2`: `xbeta()`, `zero_inflated_beta_binomial()`,
  `hurdle_cumulative()` and `cse()`.
- `aterms2`: `subset()`, `index()`, `mi(x, idx = )`, `rate()` and
  `cat()`; `nobs()` is brms's.
- `formula2`: dpar equations (`sigma1 = "sigma2"`), `cmc`, `y ~ .`, a
  list of families, and `update()` with a `bf()` delta.
- `defects`: 38 rows of the ported brms suite, among them two silent
  wrong answers (`fitted(scale = "linear")` with `cs()`, residuals of a
  missing `mi()` response).
- `postfit2`: `conditional_smooths()`, `make_conditions()`,
  `update_adterms()`, `posterior_average()`, `conditional_effects()`'s
  `surface`, `too_far`, `select_points` and `spaghetti`, a Wald band
  for a nonlinear predictor, and brms's per-row rule for which group a
  display row reads under `re_formula = NULL` (a silent wrong answer
  found by lane sampfix). `frm_bootstrap()` no longer redraws smooths.

Versions: frmtmb **0.66.0**; frmtmb.sample **0.14.0**, floor frmtmb
0.66.0 (new sampling-API exports, the ordinal inverse maps, equated
parameters); frmtmb.eam **0.11.2**, floor 0.66.0 (its `test-family.R`
reads `frmtmb:::row_aterms`; no code change); frmtmb.learn **0.7.1**
and frmtmb.ode **0.7.1** (tests only: `test-stan-identity.R` and the
scale tier's `test-scale.R` let no warning escape). frmtmb.coupling,
frmtmb.latent and frmtmb.spline are unchanged. In brms's own ported
suite, bin 1
passes 368 of 494 (306 on lane defects alone, 275 at 0.65.0).

### Decisions the user made on 2026-09-29, for this round

- **`default_prior()` and `validate_prior()` skip the response check**,
  as brms's do: a `Beta()` model on a response outside (0, 1) gets its
  table, and `frm()` still refuses the fit. An ordinal response the
  family cannot read is still refused, as brms's `extract_nthres()`
  refuses it.
- **`update(fit, data = )` is kept** and refits on the new data, as
  `stats::update()`, lme4 and glmmTMB read it. brms refuses that
  spelling. Recorded as a divergence (`brmsfit-methods:927`).
- **Laplace conditional draws stay refused.** `frm_sample(laplace =
  TRUE)` samples an approximate posterior whose use is checking the
  approximation; a call that reads the integrated values (group-level
  coefficients, a smooth's coefficients, `mi()` values) is refused by
  name, and the values are not filled in from `N(b_hat, H^-1)`.
- **`frm_bootstrap()` never redraws a smooth, `gp()` or `hsgp()`
  term**, one indexed by a grouping factor included, and
  `conditional_effects(band = "boot")` inherits the rule. This
  withdraws the whole-model default of 2026-09-24 (below).
- **Upstream bugs are listed in `dev/upstream-bugs.md`** for the user
  to file. No session files them.
- **Nobody uses the package, so break compatibility freely.** This
  round did: `nobs()`, `frm_bootstrap()`, `fitted(scale = "linear")`
  on `cs()` fits, the prior-table spellings and the ordinal draw names
  all changed, and NEWS says so under "Breaking changes".

### What the 0.66.0 round is evidence for

**Run every extension suite that can reach a changed path.** Lane
aterms2 added `subset()` and `index()` to core's compatibility table
and ran core and frmtmb.sample twice. The change broke 5 assertions in
frmtmb.eam's `test-family.R`, which the reviewer found only by running
all eight suites. `dev/lane-rules.md` now says so.

**A `test_file()` runner must attach the package.** Lane runners that
called `test_file(package = p)` without `library(p)` reported failures
in `test-conditions.R`, `test-data2.R` and `test-id-kron.R` that are
the runner's. Measured on the release build
(`dev/round-20260929b.md`): unattached, `test-conditions.R` has 3
failures and 3 errors and the other two 1 error each; attached, all
three pass. `dev/lane-rules.md` now says so.

**A small punch round can add large code, so give it a final check.**
fams2's second punch round added a "some rows" rule to the `kappa`
warning, and the final check found it warned on well-identified fits
(K1). postfit2's punch rounds replaced the whole new-level
construction of `conditional_effects()` with a per-row, per-term plan
(`R/ce-levels.R`, new), and each re-check found a new blocker in it.
Neither would have been seen without a review after the last punch
round.

**A clean textual merge can call a function another lane removed.**
sampfix's laplace probe in frmtmb.sample's `conditional_effects()`
called `ce_draw_new_levels()`, which postfit2 removed, and the two
lanes' hunks merged without a conflict. The probe is now written on
postfit2's `ce_plan_eval()`. Grep for every removed export after a
merge; the suite finds the rest.

## The 0.65.0 round, 2026-09-29

`dev/round-20260929.md` is this round's record, with the verification
of the release tree. Each lane has `dev/<lane>-findings.md` and
`dev/reviews/2026-09-29-<lane>.md`. Seven lanes merged:

- `gradcheck`: the convergence warning no longer fires on correct fits
  (289 of 720 before, 0 after); bound-aware Newton decrement.
- `csfactor`: `cs()` on a factor builds brms's treatment dummies;
  `x + cs(x)` and `mo(m) + cs(m)` refused by a rank test.
- `resmooth`: `re_formula = NA` keeps every smooth; ordinal `Est.Error`
  now carries the smooth's coefficients.
- `thresrefit`: `influence()` refits keep the fitted threshold layout
  and report failed refits.
- `arcovsample`: frmtmb.sample `log_lik()` and `loo()` for
  `cov = FALSE` ARMA; the Student-t rescor row density at large `nu`.
- `records`: the pkgcheck failures and this round's records.
- `docsci`: `.github/workflows/pkgdown.yaml` and
  `dev/release/build-docs.R`. Its deploy job skips until the owner
  switches Pages to the workflow builder; the cutover runbook is in
  `dev/docsci-findings.md`.

One more lane was opened DURING consolidation, because the release
check found frmtmb.spline's `curve-inference` vignette failing against
the new `re_formula = NA` rule:

- `splinecurve`: `allow_new_levels` through `frm_curve()`,
  `frm_curve_deriv()` and `frm_curve_feature()`, so the population
  curve of a factor-smooth model is read at an unseen level; a refusal
  of a difference between unseen `(1 | g)` levels, whose standard error
  would otherwise be exactly 0.

Versions: frmtmb **0.65.0**; frmtmb.sample **0.13.0**, floor frmtmb
0.65.0 (it calls `arma_cond_resp()` and `arma_cond_dpars()`);
frmtmb.spline **0.9.0**, floor 0.65.0 (the new `re_formula` rule and
`allow_new_levels`); frmtmb.eam **0.11.1** (a Suggests line). The
other four extensions are unchanged.

## The 0.64.0 round, 2026-09-25

`dev/parity-round-20260925.md` is the 0.64.0 round's own record. Each
lane's findings file has its validation, its numbers and its scripts.

Nine lanes merged:

- `me`: `me()` noise-free predictors and `set_mecor()`.
  `dev/me-findings.md`
- `thres`: `thres()` for the ordinal families.
  `dev/thres-findings.md`
- `grby`: `gr(g, by = )` and `mm(g1, g2, by = )`.
  `dev/grby-findings.md`
- `arcov`: `ar()`, `ma()` and `arma()` without `cov = TRUE`.
  `dev/arcov-findings.md`
- `mv`: three or more responses summed with `+`; `student()` with
  `rescor`; ordinal families in a multivariate model.
  `dev/mv-findings.md`
- `icpt0`: `0 + Intercept` and `bf(center = FALSE)`.
  `dev/icpt0-findings.md`
- `fams`: `hurdle_negbinomial()` and `zero_one_inflated_beta()`.
  `dev/fams-findings.md`
- `emm`: emmeans on nonlinear and multivariate fits.
  `dev/emm-findings.md`
- `sratio`: `sratio()` thresholds unordered, found by lane thres.
  `dev/sratio-findings.md`

Versions: frmtmb **0.64.0**; frmtmb.sample **0.12.0**, which floors on
frmtmb 0.64.0 because it needs `rescor_row_loglik()` and the new
`fit_extras(resp = )`. Every other extension is unchanged from the
0.63.0 release: frmtmb.eam **0.11.0**, frmtmb.latent **0.6.0** and
frmtmb.learn **0.7.0** keep their floor at frmtmb 0.63.0;
frmtmb.coupling **0.6.0**, frmtmb.spline **0.8.0** and frmtmb.ode
**0.7.0** keep theirs at 0.61.0.

### Verified at the release commit

Verification ran on a LINUX CONTAINER, not on the Windows box, which is
why `R CMD check` was `--no-manual` there: the container has no LaTeX.
The Windows run with the manual was done on 2026-09-28 on the fresh
install (`dev/release/check.log`; core and eam rerun in
`dev/release/check-rerun.log` after drmTMB and makeindex were
installed): core, frmtmb.eam, frmtmb.latent, frmtmb.ode and
frmtmb.spline at 1 NOTE, the environmental V8 one; frmtmb.coupling,
frmtmb.sample and frmtmb.learn OK. The ungated suite on the same tree
was 305 files, 17,893 assertions, 0 fail, 0 error
(`dev/suite-baseline.md`). The reference log-likelihoods of
`dev/machine-library.md` reproduce exactly on the container:
-295.602189818 and -332.137876329.

- Suite, core and all seven extensions, one file per process: 305
  files, 17,815 assertions. Two files fail, and both fail identically
  on the unmodified base build there: `test-pp-check-types.R` (4), from
  conda's bayesplot 1.15.0, and frmtmb.sample `test-reparam.R` (1).
- Gated: 54 files, 3,767 assertions, 0 fail, 0 error, with Stan
  compiling. `test-drmtmb-agreement.R` (131) and `test-fuzz.R` ran
  separately behind their own gates, 0 skip.
- Ported brms suite: the regenerated verdicts file and every generated
  test file are byte-identical to the committed ones. Bin 1 passes 275
  of 494, up from 252.
- `R CMD build`: OK, vignettes included.
- `R CMD check --as-cran --no-manual`: in-check tests 9,098 pass, 4
  fail, all four `test-pp-check-types.R`. The other warnings and notes
  are the container's: the locale, a missing `qpdf`, no CRAN index and
  the clock check. One warning was real, an undeclared `ordinal::` in
  `test-thres.R`; `ordinal` is now suggested.

The Windows machine got a FRESH R INSTALL on 2026-09-28, R 4.6.1, with
the user library restored from `dev/restore-cran-2026-09-28.txt` and
the manifest `dev/win-library-4.6-manifest-2026-09-28.csv` (412
packages, 402 of them from CRAN). Several Windows-side records were
stale after that and are corrected in this round: the StanHeaders pin
is gone (`dev/lane-rules.md`, `dev/machine-library.md`), and
`codemeta.json` is regenerated at 0.64.0.

**Defects the merge created, found and fixed on the branch.** A lane
can only test its own feature, so these appeared where two lanes met or
in test files a lane did not run: `ordinal_ncat()` read the first
response's family; `cs()` in a multivariate ordinal response had zero
coefficients; `lp$center` was read with `$`; two refusals shared one
message; a nonlinear body with `m[, 1]` stopped with 'argument "ei" is
missing'; nine manual port verdicts were stale. Both
`test-parity-integration.R` pins were seen failing on a build without
the fixes. The table with each cause is in
`dev/parity-round-20260925.md`.

## What is next, in order

**The docs cutover is DONE, 2026-09-29.** GitHub Pages builds from
`.github/workflows/pkgdown.yaml`: every push to main that touches the
documented sources rebuilds and deploys the eight sites, and `docs/` is
no longer tracked (it is in `.gitignore`; a local build still writes
there by default). The first deploy was of `c9d9126b`, and all 22 checks
of runbook step 5 passed on the live site. To roll back, re-track a
built `docs/` and switch Pages to `legacy` from `main:/docs` (runbook
step 7). A released site beside a development site is recorded in
`dev/docsci-findings.md` as later work, for when a tagged release
exists that people install.

**Filed at 0.66.0, in `dev/test-backlog.md`** under "Filed at the
0.66.0 release": what the six lanes left, led by a NaN the laplace
probe can pass, `predict(newdata = )` with NA responses under
`cov = FALSE` ARMA, and the bernoulli recoding of a two-valued
response.

**Filed at 0.65.0, in `dev/test-backlog.md`:** closed at 0.66.0 by
lane sampfix: `pp_check()`'s four `loo_*` types,
`posterior_predict()` and `posterior_epred()` on `laplace = TRUE`
draws, `frm_sample(laplace = TRUE)` with no random effect, and
`variables()` on the draws of an ordinal fit. Still open:
`fitted(newdata = )` on a multivariate ordinal fit errors;
`mo(m) + m` is unidentified and nothing reports it;
the finite-difference `Est.Error` at `newdata` off a `gp()`'s fitted
positions costs one evaluation pair per coefficient; standard errors
are not bound-aware at a constrained optimum.

**Then, in order:**

- **Priors on `me()` hyperparameters.** Classes `meanme`, `sdme` and
  `corme` are refused, not implemented (lane me).
- **`gr(g, by = f, cov = A)`.** Refused by name today, because brms
  correlates the by-levels through `A`, which is not a by-split
  (lane grby).
- **Latent-residual AR for non-gaussian families**, which is brms's
  other autocorrelation form (lane arcov).
- **Phases 4 and 5 of `dev/extension-gaps-plan.md` are untouched.**

**Also open, and filed in `dev/test-backlog.md`:** REML with `me()` in
`mu` is approximate and registered as conditional, and the same
argument applies to `mi()`, which is registered as working; autoscale
on covariance structures other than `us`, `diag` and Student-t;
`frm_sample()` on a one-parameter model; mixture sampling defaults for
`sigma` and `theta`; six prior spellings where frmtmb and brms disagree
on accept or refuse; missing `default_prior()` rows; a negative-hazard
penalty for `royston_parmar()`; `frm_curve()` refusing on `gamma1 ~ x`;
two core seams for frmtmb.eam (a log-difference `lcdf` slot, a refusal
of NA `dec()` on censored rows); `lba()` and `gddm()` censoring and
`lba()`/`rdm()` contaminant, not built; base under-coverage of
`log sd(log bs | s)`.

**Upstream reports, drafted and NOT filed**: RTMB `log_pnorm_both`'s
derivative (`dev/phase3b-rtmb-report-pnorm.md`); TMB `TanhOp::reverse`
(`dev/phase3b-rtmb-report-tanh.md`); drmTMB's REML and `beta_sigma`
(`dev/drmtmb-findings.md`); TMB's macOS binary and OpenMP. One more is
ready to draft: brms's `posterior_predict_hurdle_negbinomial()` does
not draw the zero-truncated negative binomial (lane fams).

**From the drmTMB comparison** (`dev/drmtmb-findings.md`), what is
left after this round: `fcor()`, boundary-corrected variance-component
tests, heritability and ICC accessors. `hurdle_negbinomial`,
`zero_one_inflated_beta` and `gr(g, by = f)` shipped at 0.64.0; a
negbinomial CDF for `trunc()` has not.

## Decisions the user made on 2026-09-29

- **Redirects** for the 104 retired core reference pages live as a
  `redirects:` block in `_pkgdown.yml`, not as static files.
- **The docs workflow has no `release:` trigger.** A push to main
  deploys; a tagged release does not redeploy an older tree. The Pages
  cutover's one-build gap is acceptable.
- **`y ~ x + cs(x)` stays refused**, in `frm()` and in `frm_sample()`
  alike. brms samples it and the prior sets the split between the global
  and category-specific parts, which the data cannot inform; `cs(x)`
  alone gives the same predictive distribution. A deliberate departure
  from brms.
- **A partial `re_formula` beside a factor smooth is ACCEPTED** and
  keeps the smooth, as brms does. This lifts the 0.63.0 refusal, whose
  reason (`NA` drops the smooth) stopped being true at 0.65.0. The
  refusal is gone entirely: its other documented case, a `car()` or
  `spde()` field, was never reachable, because a field carries the bar
  `(1 | loc)` and a formula names it like any group-level term.
- **A family added after several multivariate responses fills only the
  responses that have none.** brms gives it to every response, which
  silently overwrites an ordinal response's family; frmtmb keeps its
  own rule, a documented departure.

## Decisions the user made on 2026-09-24

- **Prior spellings brms refuses:** refuse `class = "b"` where the
  predictor has no slope, `rescor` with `resp`, and `resp = "y"` on a
  univariate model. Keep `cor` with `resp` and `Intercept` with `nlpar`
  as frmtmb extensions.
- **Several-location extension families** (`lca`, `hmm`, `lba`, `rdm`,
  `mixture_mvn`) need `dpar` on `b` and `Intercept` priors, frmtmb's own
  rule where brms has no such family (organizer's call, not objected to).
- **`frm_bootstrap()` keeps its whole-model default**, redrawing group
  effects and smooths; `re_formula = NULL` conditions on them. No new
  argument. WITHDRAWN 2026-09-29 for smooths: `frm_bootstrap()` never
  redraws a smooth, `gp()` or `hsgp()` term (see the 0.66.0 round).
- **Residual correlation at newdata counts fitted time levels**, a
  documented departure from brms, whose position-based reading is not
  consistent under marginalization.
- **Random-effect smooths under `re_formula = NA` should match brms**
  (keep them). Filed, not built.
- **Keep reporting nonzero optimizer convergence codes.**
- **`autoscale` defaults to `NULL`** (BREAKING), engaging below sd 1e-3,
  or 0.05 for a random-slope column, and falling back to the plain fit
  when its pre-fit fails.
- **eam contaminant window:** `contaminant_range =`, or with a `trunc()`
  upper bound, the fastest response to the deadline; otherwise refused.
- **eam left and interval censoring** use the boundary `dec()` names.
- **`rp_floored()`** warns when every flagged row is censored and refuses
  on a flagged event row.

**Settled earlier, do not reopen:**
- The tiebreaker: match brms, UNLESS it is clearly obvious that brms
  should be doing it the other way.
- `re_formula = ~(1 | nosuch)` stays refused. (The companion rule, that
  a partial `re_formula` stays refused on a fit with a factor-smooth
  term, was REVERSED by the user on 2026-09-29; see above.)
- `propagate_error` is the name, spelled out.
- The non-generic name collisions with brms stay, because `::` is
  sufficient in both load orders.
- A gratia older than 0.9.0, loaded before frmtmb.sample, stops it
  loading; the user accepted that.
- `inverse.gaussian` defaults to brms's `1/mu^2`.
- Duplicate priors on one slot, and duplicate group-level effects, are
  refused as brms refuses them.
- `frm_simulate(newparams =)` takes brms names only.
- `variables()` keeps frmtmb's order; frmtmb objects carry no brms class.
- `REML = TRUE` integrates the `mu` coefficients only.

## How a round runs here

Lanes in manual git worktrees off main, one item or one coherent group
per lane. A worker writes. A reviewer whose job is to FALSIFY rather
than confirm reads the same worktree against a shared reference build
of the base commit. Punch rounds go back to the worker, and a reviewer
re-checks.

    git worktree add ../frmtmb-wt-<name> -b wt-<name> <sha>

Two punch rounds is the cap, a third only for a blocker. The previous
release's library served as the reference build again this round.

Consolidation: commit each lane on its branch, merge, set versions,
roxygenise, THEN install (that order), verify the tree is unchanged by
roxygenise, run the tiers, run each suite in one session three ways
(below), commit the release and the docs separately, remove the
worktrees, prune branches, regenerate `dev/suite-baseline.tsv`. Do not
install into the release library while any lane still uses it as its
base build.

**Run each suite in one session: forward, reverse and preload.** The
tiers run one test file per R process. covr and `R CMD check` run a
whole suite in one session, in alphabetical order. 0.71.0's
test-coverage job failed only there: one file loaded frmtmb.eam, and a
later file assumed that `dec()` was not registered. With every
extension in the release library, run:

    NOT_CRAN=true bash dev/release/run-sessions.sh rel<ver> 8

It runs all eight suites through `dev/release/run-session.R` three
ways: alphabetical order, reverse order, and preload (alphabetical,
with all seven extensions loaded before the first file): 24 sessions.
Preload is the session in which no file loads an extension, so every
registry difference it reports is a test's own, and a test that
assumes an extension is absent fails there. The CI coverage job
measures core with no extension installed (`test-coverage.yaml`), so
this gate is the only run of core's suite in one session with the
extensions present.

The run passes when the last line of `dev/release/rel<ver>.log` reads
`ran 24 of 24`, every SESSION line reads `fail=0 err=0`, and every
LEAK line ends in `[with load of ...]`. A LEAK line without it is a
test that left a registry, the notice state, an option, an environment
variable, the search path or the global environment changed for the
files after it; a global or search-path LEAK is never tagged. When the
release library changes, change the `LIB` line of `run-session.R`
with the one in `run-tests.R`.

**Run `R CMD check` on a QUIET machine.** Its examples-timing NOTE
measures load here; it did not appear at the 0.63.0 release, run alone.

**Memory.** On the machine of 2026-09-28, which has 63 GB of RAM and a
fast disk, test runs need no process cap. Keep ONE TEST FILE PER R
PROCESS and run as many at once as the work needs. The history stays
because the failure mode is real on a smaller box: in the 0.63.0 round
one lane running 11 R fitting processes at once crashed the machine of
that time, and 34 of its fits had already failed with `std::bad_alloc`
(`dev/lane-rules.md`).

## What the user has settled

- **The user pushes.** No session pushes, and no git operation without
  their say so.
- **Nobody uses this package yet.** Break backward compatibility
  freely: ship the refusal, bump, and say plainly in NEWS what stops
  working.
- **Nothing pins StanHeaders.** The `pinlib` at 2.32.10 was retired on
  2026-09-17. rstan 2.32.7 compiles against the user library's
  StanHeaders 2.39.1 because
  `C:/Users/adf44/Documents/.R/Makevars.win` carries
  `CXX17FLAGS += -std=gnu++17`, and tmbstan 1.2.1 no longer samples a
  standard normal. That file is outside `%LOCALAPPDATA%`, so a library
  loss does not touch it (`dev/tmbstan121-findings.md`).
- **The R user library stays under `%LOCALAPPDATA%`** by the user's
  decision, knowing it will recur. RTMB 2.0 comes from r-universe after
  a restore, because CRAN's Windows binary is 1.9.
- Core is the user's lane except where they ask otherwise.

## What the 0.65.0 round is evidence for

Seven lanes, each with an adversarial review and one or two punch
rounds.

**The reviews found false sentences more often than false numbers.** In
four of the five code lanes the first review blocked on documentation:
a manual page that promised every short fit would warn, a help page
that said `cs()` is not available in a multivariate fit, a help page
that described which rows of an ARMA group are shifted, and a findings
file whose reasoning was wrong in both directions while its outcome was
right. Each was measured false by construction. Review the prose
against the code with the same care as the code.

**A guard failed closed again.** The withholding of the gradient
verdict on importance-corrected fits also withheld the clean line on
every such fit, well converged or not, where 0.64.0 printed it. The
worker had argued in a code comment that withholding was the safe side.
The rule stands: construct the case where the guarded thing is ABSENT,
and a default must never be less diagnostic than what it replaces.

**A count settles what a reason cannot.** A reviewer attributed a
28.9x slowdown to a batching routine returning NULL. The worker counted
batches and model evaluations and showed the routine never returned
NULL: it split per column position, which on a smooth is one batch per
coefficient. The reviewer withdrew its cause on the page.

**The `tail -f` trap bit the consolidation itself, twice.** A watcher
tailing the four tier logs made the PowerShell drivers' `Add-Content`
throw, exactly as `dev/lane-rules.md` records. Stopping the watcher did
NOT release the logs: a monitor that is stopped or expires leaves its
`tail` process running, and ten such orphans from the session's
earlier watchers were still holding files open, so the relaunched
gated tier ran all 39 files, logged none of them, and printed "GATED
ran 39 of 39". Kill `tail.exe` before trusting a log, and wait on the
driver PROCESS, never on the log.

**A lane that changes core behavior must run every extension's check,
not only their suites.** The resmooth lane updated frmtmb.spline's help
pages and `test-curve.R` for the new `re_formula = NA` rule and ran the
extension suites green, but never built frmtmb.spline, whose
`curve-inference` vignette drew a population curve at `re_formula = NA`
on a grid with no grouping column. The release check found it; the fix
needed a new `allow_new_levels` argument through the `frm_curve()`
family, a lane of its own (`dev/splinecurve-findings.md`).

**A refusal can outlive its reachable cases.** Lifting the factor-smooth
case of the partial-`re_formula` refusal left its other documented
case, a `car()` or `spde()` field, which a new guard test showed was
never reachable: both carry a `(1 | loc)` bar. The help page had said
so since 0.63.0. A guard test for a surviving branch is what found it.

**A fresh machine has no git identity.** Commits could not be made
until the user configured one, so the consolidation built and verified
the integrated tree in `../frmtmb-wt-release` first and committed
afterwards.

## What the 0.64.0 round is evidence for

Nine lanes, consolidated and verified on a Linux container.

**A lane cannot test the seam it shares with another lane.** Six
defects came out of the merge itself, not out of any lane, and two of
them were wrong answers rather than refusals. The consolidating session
found them by running the whole suite on the merged tree, which is the
only place they exist.

**Verification on one platform leaves a hole.** `--no-manual` on the
container skipped the manual sections, which is where an unescaped `%`
in Rd has bitten this project before, and the pkgcheck failures this
lane fixed were on a Windows-side check that the container run never
reached. A record round after a platform move is not optional.

## What the 0.63.0 round was evidence for

Five lanes, each with at least one adversarial review and one or two
punch rounds; one machine crash and one power loss.

**Review still finds silent wrong answers at about one per lane.** Most
of the defects that round listed were found by a reviewer, not by a
tier, and several predated the round. Until that rate falls, a result
should not be trusted without a check against brms or an exact
reference.

**A shared seed is a hidden replicate count.** Phase 3b reported a
drift-intercept interval covering 84.8 percent over 231 fits. The arms
drew their random effects right after the same `set.seed()`, so the 231
fits were 80 seeds, and an interval that knew every true drift did no
better. Fresh seeds on the released build cover 65 of 70. Coverage
claims now use seeds independent across arms, or say they are shared.

**A fix can be worse than the defect it fixes.** The autoscale default
closed a 62-unit silent stall and opened a path where a separated fit
reported success with no warning. The rule that closed it is general:
a default must never be LESS diagnostic than the setting it replaces.

**Interrupted runs are a standing condition, not an accident.** A crash
and a battery shutdown each cut every lane mid-tier that round. What
held: treat every log after the cut as void, check each saved result
reads back, verify each lane library against its source before trusting
it.

## Worktrees

None should be live once 0.67.0 is committed. The 0.67.0 round used
three lane worktrees off `57c25589` (`wt-ordinal`, `wt-ceplot`,
`wt-formrobust`) and one integration worktree, `wt-release` off
`57c25589`, where the three were combined and verified. Their private
libraries are `C:/Users/adf44/source/r/wt-<lane>-lib`; the release
library is `C:/Users/adf44/source/r/rellib-r5`, and `rellib-r4` holds
the 0.66.0 reference build. Remove a worktree and prune its branch only after its
work is merged and its evidence is committed on main under `dev/`.
