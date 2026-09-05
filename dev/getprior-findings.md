# get_prior route= lane, findings

## Environment
- worktree C:/Users/adf44/source/r/frmtmb-wt-getprior (wt-getprior @ a233c3c)
- The user library C:/Users/adf44/AppData/Local/R/win-library/4.6 is HALF BROKEN:
  350 directories, only 138 valid installs (RTMB, TMB, reformulas, testthat,
  lme4, brms all absent as far as R is concerned). Not mine to repair.
- Private lib gp-lib holds brms, tmbstan, rstan, posterior, loo, coda,
  bayesplot and the worktree core. Remaining deps come from the sibling
  lane's is-lib, read-only. Env: scratchpad/gp-env.sh.

## Defect reproduced (dev/getprior-load-order-probe.R)
Both get_prior(formula, data) and get_prior(fit) differ before/after
library(frmtmb.sample): (flat) everywhere becomes
student_t(3, 288.7, 59.3) Intercept / student_t(3, 0, 59.3) sigma+sd /
lkj(1) cor. Confirmed FALSE for identical() on both calls.

## API as shipped
get_prior(formula, data = NULL, family = NULL, data2 = list(),
          route = c("fit", "sample"))

Help text for `route`:
  Which route's defaults the `prior` column reports. `"fit"` (the
  default) reports the defaults `frm()` applies, which are flat in
  every slot; it consults no registry, so its answer does not depend
  on which packages are attached. `"sample"` reports the defaults
  `frmtmb.sample::frm_sample()` applies, and refuses when no package
  has registered any. The returned object records the route and
  `print()` names it on its first line.

Returns class c("frmtmb_prior_rows", "data.frame") with attr "route".
print.frmtmb_prior_rows prefixes one line:
  route = "fit": the prior defaults frm() applies
  route = "sample": the prior defaults frm_sample() applies

Refusal (R/sampling-api.R require_prior_defaults()):
  route = "sample" reports the prior defaults frm_sample() applies,
  and no loaded package states them. They live in frmtmb.sample: run
  library(frmtmb.sample) and call again (install it from
  extensions/frmtmb.sample in the frmtmb repository). For the defaults
  frm() applies, use route = "fit".

Test seam: swap_prior_defaults(providers = list()) in R/sampling-api.R,
internal (@noRd), returns the previous provider list.

## Results so far
- extension suite (one process, load_all, core from gp-lib2):
  881 pass, 0 fail, 0 error, 2 skip (opt-in FRMTMB_BRMS_FIT_TESTS).
- new core file test-get-prior-route.R: 29 pass, 0 fail.

## Verification after the fix
- dev/getprior-load-order-probe.R against the patched core:
  get_prior(formula) identical before/after attach: TRUE
  get_prior(fit)     identical before/after attach: TRUE
  (with the INSTALLED frmtmb.sample, i.e. an extension build that
  predates this change still cannot perturb the fit route.)
- route = "sample" with core alone: refused, message names
  frmtmb.sample and library(frmtmb.sample) and route = "fit".
- route = "sample" with frmtmb.sample: the brms defaults, unchanged
  from what the old unconditional lookup produced.

## Call sites, and what each was decided to be
core R/:  no internal caller of get_prior(); mentions in roxygen only.
core vignettes/: none.
core tests (all FIT route, unchanged):
  test-data2.R:259, test-nlf.R:175, test-prior-compat.R:171/204/373,
  test-priors-autocor-classes.R:182/192/196, test-v19.R:78/89/96.
  Every one asserts WHICH SLOTS a design offers (class/coef/group/
  dpar/nlpar/resp columns); none reads the `prior` column. The fit
  route is what they meant and what they still get.
  test-brms-agreement.R:463 and test-prior-compat.R:439 call
  brms::get_prior(), not ours. Untouched.
  test-bracket-access.R:18 is a comment mentioning get_prior. Untouched.
frmtmb.sample:
  tests/testthat/helper-brms.R:4 is a COMMENT listing brms functions;
  there is no get_prior() call in that file. Nothing to change.
  vignettes/sampling.Rmd:67 MEANT the sampling defaults: rewritten to
  show both routes and to pass route = "sample" for the sampling one.
  R/zzz.R:13 comment rewritten.
frmtmb.ddm tests/testthat/test-surface.R:121 (outside my surface):
  fit route, asserts dpar/class columns only, needs no change.

## Check note
First R CMD check pass errored only on "Packages suggested but not
available: clubSandwich, fmesher, marginaleffects, mclust, metafor",
an artifact of the private library, not of the change. Installed them
into gp-lib and re-ran the check against the same tarball.

## Runner hazard, recorded
TaskStop on a background bash job kills the shell but NOT the loop it
spawned. A stale gp-suite.sh kept running and interleaved into the same
summary and log directory, which produced one bogus
"test-core-boundary.R CRASHED rc=0" line. Killed by PID (8284) after
confirming sibling lanes (islap, rv9) had their own R processes that
must not be touched. The live loop's own line for that file is clean.

## Test results
Core, one file per process, load_all, NOT_CRAN=true:
  102/102 files, 5272 pass, 0 fail, 0 error, 6 skip.
  Skips: test-brms-agreement.R 3, test-case-studies.R 1,
  test-famgaps.R 1, test-fuzz.R 1 (all pre-existing opt-in gates).
  New file test-get-prior-route.R: 29 pass, 0 fail.
  test-message-uniqueness.R 6 pass, test-bracket-access.R 7 pass.
frmtmb.sample, whole suite in one process, load_all, core from gp-lib2:
  881 pass, 0 fail, 0 error, 2 skip (FRMTMB_BRMS_FIT_TESTS opt-in).
  New file test-prior-route.R: 9 assertions, all pass.
frmtmb.sample R CMD build: vignettes render OK after the sampling.Rmd
  rewrite.

## R CMD check
Pass 2 (full Suggests): 1 WARNING, mine:
  "checking for unstated dependencies in examples ... WARNING
   'library' or 'require' call not declared from: 'frmtmb.sample'"
  Cause: the two-call example guarded the sample route with
  requireNamespace("frmtmb.sample"). Core cannot Suggest an extension
  it does not depend on, so the sample-route half moved into
  \dontrun{}. Pass 3 re-run to confirm.

## R CMD check, final
Pass 3, after moving the sample-route example into \dontrun{}:
  R CMD build OK, R CMD check --as-cran --no-manual with
  _R_CHECK_CRAN_INCOMING_=false: Status: OK.
  examples 23s OK, examples --run-donttest 42s OK,
  tests testthat.R 242s OK, re-building vignettes 115s OK.

## Deliberate omissions
- DESCRIPTION not touched: core cannot Suggest frmtmb.sample (an
  extension that depends on core and is not on CRAN), so the
  sample-route example is \dontrun{} rather than requireNamespace().
- frmtmb.ddm's get_prior() call left alone: fit route, and outside the
  assigned surface.
- Core's own get_prior() tests left alone: all fit route, none reads
  the `prior` column.
- swap_prior_defaults() left unexported: a test seam, not API.
- The `route` value is validated by match.arg(), adding no new message
  template; test-message-uniqueness stays green.

## Review fixes (GO-WITH-FIXES, dev/review-getprior.md)
Reviewer's own comment edit at R/priors.R print method (ROW vs COLUMN
subsetting) kept verbatim.

1. Column subset: KEPT the current behavior, added the test.
   A `[.frmtmb_prior_rows` that carried `route` through column subsets
   would let the label survive onto a selection that has dropped the
   `prior` column, i.e. a table with no route left to describe. Extra
   S3 surface for a worse claim. tests/testthat/test-get-prior-route.R
   now asserts gp[, c("class","coef")] keeps the class, loses the
   attribute, and prints through the unlabeled branch.
2. as.data.frame.frmtmb_prior_rows() added (R/priors.R:853), drops the
   attribute with the class. print() now delegates to it instead of
   stripping by hand. The extension test drops its bare() helper and
   uses as.data.frame().
3. dev/getprior-load-order-probe.R walks up to find
   extensions/frmtmb.sample, so the load_all fallback resolves from the
   repo root and from dev/. Verified in both, with the extension absent
   from the library path so the fallback is actually taken.
4. \dontrun{} example replaced by prose in the route section naming
   route = "sample" and what it prints. Fit-route example stays
   runnable. man/get_prior.Rd has no \dontrun.

Reviewer report corrections accepted: the extension test uses no `:::`;
"survives subsetting" is rows only; the 5272 pass count is
environment-dependent and is quoted with its environment.
