# Lane smallfix: four small backlog items

Round opened 2026-10-24, base 18b536fc (frmtmb 0.74.0 plus two CI-only
commits). Worktree `frmtmb-wt-smallfix`, library `wt-smallfix-lib`
(this lane's frmtmb, frmtmb.latent, frmtmb.survival and frmtmb.spatial,
installed from the worktree), then `rellib-r12` (the base build, which
also serves the five extensions this lane does not change), then the
user library. "Before" is `rellib-r12` alone.

Each item is a separate diff in `dev/smallfix-out/` (gitignored, local;
`bash dev/smallfix-diffs.sh` writes them again and checks that every
changed file is in exactly one diff):

| Item | Diff |
|---|---|
| S1, a saved fit names every missing package | `S1-saved-fit-packages.diff` |
| S2, mixture_mvn test debt | `S2-mixture-mvn-tests.diff` |
| S3, two refusal messages | `S3-refusals.diff` |
| S4, the CI skip gate for core | `S4-ci-skip-gate.diff` |
| NEWS of all four | `news.diff` |
| this record and the lane's runners | `harness.diff` |

The items touch disjoint files except `NEWS.md` and frmtmb.latent's
`NEWS.md` (which bullet is whose is said per item) and, since fix
round 1, three R files: S1's handler rewrite touches `R/priors.R` and
`R/curve.R` (S3's) and `R/predict-brms.R` (S2's). Their hunks are split
between the two diffs by `dev/smallfix-split-hunks.R`. The script
checks that each item diff applies alone to a clean archive of HEAD,
and that all of them together give this worktree (0 files differ).
Two dependencies stay inside an item: S2's latent test of the
`predict(newdata =)` message needs S2's `R/predict-brms.R` hunk, and
S1's new core test needs S1's three-copy helper. No item needs another
item.

## S1. A saved fit that needs two extensions names every missing package

### What changed

`R/extension-home.R`. The two checks keep their trigger: the family
check fires only when a family's package (or an old core family, or a
mixture component's package) cannot be used, and the structure check
only for a block's structure. A method that needs nothing of a missing
package still answers (`logLik()`). Once one fires, the refusal is
built from every part of the object that this session cannot use:
families (`family_missing_parts()`) and structures
(`structure_missing_parts()`), through one builder,
`stop_extension_parts()`. Class `frmtmb_extension_missing` is kept.

- One part: the message is word for word the old one for a missing
  family, a missing, session or stale structure, and an old cox()
  family whose home is installed (the review compared the texts:
  `dev/reviews/2026-10-24-smallfix.md`, finding 4). It is NOT the old
  one for an old family whose home is not installed ("Install
  frmtmb.survival and refit ...", where 0.74.0 said "Refit ...") or,
  since fix round 1, installed but not loadable ("frmtmb.survival is
  installed but cannot be loaded in this session. Reinstall ...").
  Every existing test that matches a one-part text still passes
  (frmtmb.survival's `test-saved-fits.R`, frmtmb.spatial's
  `test-guardrail.R`, core's `test-covstruct-seam.R`).
- Several packages: "This model needs 2 extension packages that are not
  installed or cannot be loaded in this session, so their functions are
  not available: frmtmb.survival (the cox family of this model) and
  frmtmb.spatial (the car() structure of this model). Install them and
  run library(frmtmb.survival); library(frmtmb.spatial)".
- An old core family (`old_core_family_home()`): its two sentences
  first, then the other missing parts, then one instruction: "Install
  <every package not installed> and refit the model with
  library(<home>); library(<others>) loaded". The home is asked whether
  it is installed (`system.file()`), not loaded, so reading an old fit
  does not register an extension's families as a side effect.
- A structure registered in a session, or one whose package is loaded
  but does not register it, keeps its own sentence and instruction.

"One install line": the refusal names every package in one
instruction. It does not write an `install.packages()` or `remotes`
call: no extension is on CRAN, and a package recorded by
`frmtmb_register_covstruct(package =)` can be a third party's. The
README holds the install commands, as the guardrail-1 refusals assume
("Install it and run library(...)").

Found while testing: `predict()` of such a fit first warned "the
covariance of the estimates could not be recovered from the Hessian",
and then refused. `fit_draw_space()` (`R/brms-shapes.R`) read the
refusal from `vcov()` as a covariance that failed. It now calls both
checks first, so `predict()` refuses at once
(`dev/smallfix-predict-probe.R`: before, "WARNING in predict" and then
the refusal).

`tests/testthat/helper-absent-package.R` is new in core: a third copy
of `dev/rel074-helper-absent-package.R`, which now takes several
packages (`unloadable_package_libs(c(a, b))`,
`absent_package_libs(c(a, b))`); one package behaves as before. The two
extension copies are updated to the same bytes, and
`absent-helper-identical.yaml` checks three copies (its body run
locally: "3 copies, 1 distinct hash(es)",
`dev/smallfix-out/helper-identical.txt`).

### The record, rerun

`dev/rel074-saved-mixed.R` (seed 20261023, the 0.74.0 cox() and car()
fit, logLik -278.375375098, conv 0), read in four libraries holding
this lane's (after) or rellib-r12's (before) frmtmb, with rellib-r12's
extensions (`bash dev/smallfix-saved-mixed.sh before|after`, logs
`dev/smallfix-out/saved-mixed-<arm>.txt`, tabulated with awk). 13
methods; what each refusal names:

| library | before | after |
|---|---|---|
| neither extension | 8 name frmtmb.spatial, 4 frmtmb.survival; logLik ok | 12 name both; logLik ok |
| frmtmb.survival only | 12 name frmtmb.spatial; logLik ok | the same |
| frmtmb.spatial only | 9 name frmtmb.survival; logLik, ranef ok; 2 other errors | the same |
| both | 10 ok, 3 other errors | the same |

The other errors are the record's: `VarCorr()` "does not contain
covariance matrices" (a car() field is not a grouping-factor block) and
cox()'s own refusals of `simulate()` and `fitted()`. No method gives
"could not find function".

### Tests

`tests/testthat/test-extension-home-mixed.R` (new, 4 blocks, 166
expectations):

- Two absent test packages own the family and an `mxiid()` structure
  (in process): each of 9 methods refuses with the exact two-package
  message; with the family core's, with the one-package message.
- The family's package alone absent and the structure registered:
  every refusal names the family's package only. The guard absent
  (nothing missing): all 9 methods answer.
- An old core cox() family with an absent structure package: the exact
  combined message from both checks. The guard absent: a core block
  passes the structure check, and the family refusal names no
  structure.
- The real case in child processes: a cox(df = 4) and
  car(type = "escar") fit (seed 20261024, 3 x 3 lattice, 108 rows)
  saved in a child with both extensions, then read with both
  unloadable (stubs, which R CMD check can build) and both not
  installed: 12 calls, 11 refused, each naming both packages and the
  one instruction. Then each package alone unloadable: every refusal
  names that package and not the other. ("One absent" uses the stub,
  because a stack without one extension but with the other cannot be
  built where both share a library, as in rellib-r12 and on CI.) It
  skips with "the saved-fit arm needs these packages installed: ..."
  where either extension is not installed: core's CI jobs (the
  extensions are removed from the coverage job, and R CMD check
  --as-cran hides undeclared packages; see S4).

Seen to fail on rellib-r12 (`dev/smallfix-out/before.log`, the new
file on the base build): pass 111, fail 55, for example "Expected
`conditionMessage(e)` to be identical to `msg2` ... actual: The mxiid()
structure of this model is from the extension package
frmtmbabsenttestpkgb", and in the child arm "CALL summary :: TRUE ::
The cox family of this model is from the extension package
frmtmb.survival, which is not installed ...".

`test-extension-home.R`: one expectation now matches "efit the model
with library(frmtmb.survival) loaded", because where frmtmb.survival is
not installed (R CMD check) the instruction starts "Install
frmtmb.survival and refit". Count unchanged, 31.

### NEWS

Core "Bug fixes": the first two bullets.

### Fix round 1 (review `dev/reviews/2026-10-24-smallfix.md`)

**Finding 1, every method.** The first pass made the two checks build
one refusal, but a method that reached one part on its own still named
one package: `hypothesis()` through `cs_entry()` with no fit,
`residuals()` through `response_mean()` on the family alone, and
`drop1()` through the parser's per-special refusal. I walked every S3
method frmtmb registers for `frmtmb_fit` (NAMESPACE, 92 generics):
`dev/smallfix-methods-probe.R` calls each on the saved cox() and car()
fit of block 4 (seed 20261024), read where both extensions are
unloadable stubs and where both are loaded (logs
`dev/smallfix-out/methods-probe-{none,all}-r0.txt` before, `-none.txt`
and `-all.txt` after). Before: 12 methods refused naming both, 5 named
one (`hypothesis`, `ngrps`, `parnames`, `residuals`, `variables`),
`drop1` gave the parser's car() refusal, `cooks.distance` said "all 108
deletion refits failed", `conditional_smooths` and `pp_mixture` gave
their own unrelated refusals first, and `fixef` answered with 6 NA.
The fit-level
checks, as one helper `require_fit_packages()`, now head
`hypothesis()`, `residuals()`, `drop1()`, `fixef()`, `ngrps()`,
`parnames()`, `variables()`, `influence()` (so `cooks.distance()`),
`profile()`, `conditional_smooths()`, `pp_mixture()`, `plot()`,
`emm_basis()`, `recover_data()` and every fit of `anova()`
(`dev/smallfix-head-checks.R` inserted the first 14). After: 22 refuse
naming both, none names one, none answers with NA. The 17 that answer
read only stored values (`logLik`, `nobs`, `formula`, `family`,
`model.frame`, `coef`, `sigma`, ...); 27 others stop for their own
reason (posterior-only methods, `anova()` with one fit). 13 generics
of insight, marginaleffects and lme4 were not called by the probe
(their generics are not frmtmb's); their methods go through `vcov()`,
`predict()` or `frm_linpred()`, which refuse naming both. 13 more need
a setup the probe does not build (`kfold`, `loo_subsample`, ...): the
posterior ones refuse a fit for their own reason, and `influence`,
`profile`, `plot`, `emm_basis` and `recover_data` have the head check.

**Finding 2, the handlers.** `dev/smallfix-handlers.R` walks the
installed namespace's call graph (codetools) back from the functions
that raise the refusal (331 functions can reach one) and lists every
`tryCatch()` or `try()` in `R/` whose expression calls one of them and
whose error handler does not signal the refusal again: 72 handlers
(`dev/smallfix-out/handlers-before.txt`). A new helper,
`keep_missing(value)` (`R/extension-home.R`), returns the handler's
old fallback for any other error and signals this refusal again; 70
handlers were rewritten to it by `dev/smallfix-keep-handlers.R`
(simple fallbacks `NULL`, `integer(0)`, `NA_real_`, the condition
itself), and two by hand (`curve_linkinv()`'s `try()`, and
`fit_end_checks()`'s warning handler). The scan after: 0. Per file and
function (count of handlers):

- `R/allfit.R`: `frm_allfit` (1)
- `R/autoscale.R`: `autoscale_prefit` (1), `unit_theta_start` (1)
- `R/bootstrap.R`: `boot_refits` (1)
- `R/brms-shapes.R`: `fit_draw_space` (3)
- `R/conditional-effects.R`: `ce_boot_draws` (1)
- `R/confint.R`: `diagnose` (1), `diagnose_singular` (2),
  `flat_par_display` (1), `smooth_edf` (1)
- `R/curve.R`: `curve_linkinv` (1)
- `R/families-extreme.R`: `gev_fit_check` (1)
- `R/families-simplex.R`: `lncor_summary_rows` (1)
- `R/families.R`: `xbeta_kappa_se` (1)
- `R/fit-end.R`: `fit_end_checks` (1)
- `R/fit.R`: `bound_held_by_step` (1), `check_convergence` (1),
  `headroom_drop_rows` (2), `headroom_reparam_hessian` (1),
  `newton_polish` (1), `nl_flat_message` (1), `quad_breakdown_message`
  (1), `sd_up_escape` (3)
- `R/influence.R`: `influence.frmtmb_fit` (1)
- `R/insight.R`: `get_predicted.frmtmb_fit` (1)
- `R/interop.R`: `coef_vcov_root` (1), `cov_held_as_factor` (1),
  `emm_basis.frmtmb_fit` (1), `interop_vcov` (1)
- `R/me.R`: `me_hyper_table` (1)
- `R/methods-fit.R`: `fixef.frmtmb_fit` (1), `ranef.frmtmb_fit` (1),
  `summary.frmtmb_fit` (1), `summary_fixed_frame` (1),
  `summary_gp_frame` (1), `summary_mo_frame` (1),
  `summary_random_list` (1), `summary_spec_frame` (2)
- `R/predict-brms.R`: `fd_eta_chain` (1), `fitted_old_levels` (1),
  `predict_cross_rows` (1), `predict_new_level_spec` (1),
  `re_row_support` (1), `re_used_b` (1)
- `R/predict.R`: `fitted_point_se` (1), `residual_point_se` (1)
- `R/priors.R`: `prior_drop_unmatched` (1), `resp_missing_refusal` (1)
- `R/refit-data.R`: `fit_prior_density` (1)
- `R/sandwich.R`: `cluster_root` (1)
- `R/se-check.R`: `se_boundary_blocks` (1), `se_boundary_labels` (1),
  `se_boundary_names` (1), `se_boundary_stop` (2), `se_check` (1),
  `se_edge_sd` (1), `se_flush_deferred` (1), `se_sd_scale` (1),
  `separation_check` (1)
- `R/se-reparam.R`: `cov_products_lossy` (1), `joint_cov_whole` (1),
  `reparam_basis` (1), `se_reparam_cov` (1), `se_reparam_defer` (1)

Two counts, not to be confused: 72 handlers can reach the refusal (the
scan), and `R/` has 90 `keep_missing()` call sites outside its
definition (`grep -rn "keep_missing(" R/`). The 18 more are simple
handlers inside a listed function whose expression cannot reach the
refusal, which the rewrite changed with the rest (it works per
function); they still fall back as before, since `keep_missing()`
differs only for this one class. Many sit on
the fitting path (`R/fit.R`), where a fresh fit has every package by
construction; they are rewritten so that no handler that can reach the
refusal is left, not because a fit can meet it there.

Test: `test-extension-home-mixed.R` gains `fixef`, `hypothesis`,
`residuals`, `drop1`, `ngrps` and `variables` in the in-process method
list and `fixef`, `hypothesis`, `residuals`, `drop1` in the child block
(16 calls, 15 refused, all naming both), and a block where an old core
cox() family (the 0.73.0 fit) is refused by `fixef()`, `vcov()`,
`summary()` and `hypothesis()` whether or not frmtmb.survival is
installed, with the guard absent (finite `Est.Error` on the fit as it
was). This defect is in 0.74.0 too (the review's
`s1b/new-coxcar--base_none.txt`: "CALL fixef :: OK Est.Error: NA,NA");
NEWS has a bug-fix bullet for it.

**Finding 3, the instruction.** "Install it" is used only when the
sentence naming the package is the one before the instruction;
otherwise the packages are named ("... before using the fit. Install
frmtmbabsenttestpkga and run library(frmtmbabsenttestpkga)"). The
install list and the `library()` list now come in one order (old
families' homes first, then the missing packages). Test: "the
instruction names the package after another sentence".

**Finding 4, the record and the stub case.** The "word for word" claim
above is corrected. An old family's home that is installed but cannot
be loaded is told apart without loading it: an installed package has
`Meta/package.rds`, and a copy without it gets "frmtmb.survival is
installed but cannot be loaded in this session. Reinstall
frmtmb.survival and refit the model with library(frmtmb.survival)
loaded". Test: "an old family's home that cannot be loaded is to be
reinstalled", in child processes (stub, whole installed copy, not
installed), each asserting the home is still not loaded afterwards.
A missing (not old) package keeps its one-part text, "is not installed
or cannot be loaded ... Install it and run library(P)", which the
extensions' saved-fit tests match; "install" there covers reinstall.

**Finding 5, a CI job for the real case.** The pkgdown job, which
installs every extension, gets a step after the compatibility-rows
step: `Rscript .github/scripts/saved-fit-check.R`, which runs
`test-extension-home-mixed.R` and fails on any failure, error, skip or
an empty result, as that step does. On the runner the extensions sit
in the one library that also holds frmtmb's dependencies, so the
helper's "not installed" stack (every library that holds neither
extension) would drop frmtmb's dependencies, and that arm would skip.
The script moves frmtmb.survival and frmtmb.spatial into a library of
their own beside that one (a rename) for the test's duration, puts it
first in `.libPaths()`, and moves them back in a `finally`. The child
processes are started by `child_with_libs()` with `include.site =
FALSE` and the parent's whole `.libPaths()`, so R_LIBS_SITE cannot put
a package back and loses nothing. Locally, in a library holding this
lane's frmtmb and both extensions, with the user library behind
(`dev/smallfix-out/saved-fit-check.txt`): 7 blocks, 256 expectations,
0 failed, 0 skipped, exit 0, both packages back in place. After the
short check the `finally` also removes the moved-to library when it is
empty (left in place, with its contents, when a package could not be
moved back), and the step's comment says that a failure blocks the
Pages deploy, as the compatibility-rows step's does, and why that is
accepted. Rerun on fresh scratch copies
(`dev/smallfix-out/saved-fit-check2.txt`): 7 blocks, 256
expectations, 0 failed, 0 skipped, exit 0, both packages back, and no
`frmtmb-saved-fit-ext-lib` left. The trigger paths gain the test, its
helper and the script. The per-file
merge with main's uncommitted `pkgdown.yaml` edit is clean (`git
merge-file` on base 18b536fc, this lane's file and main's working
file: 0 conflicts; `dev/smallfix-out/merge/`). Not verified: the step
on the runner.

**Finding 6.** `absent-helper-identical.yaml` hashes the master,
`dev/rel074-helper-absent-package.R`, with the three copies (4 files)
and lists it in `paths:`. Its body run locally: "4 files, 1 distinct
hash(es)" (`dev/smallfix-out/helper-identical2.txt`).

## S2. mixture_mvn test debt in frmtmb.latent

### Tolerances

`extensions/frmtmb.latent/tests/testthat/test-mixture-mvn.R`. Every
absolute tolerance is now a ratio to something the run measures. For
each, the right answer and a constructed wrong answer were measured
against the same bound (`dev/smallfix-mvn-ratios.R`, the tests' seeds
and data; log `dev/smallfix-out/mvn-ratios.txt`):

| was | now | right | wrong (constructed) | bound |
|---|---|---|---|---|
| \|ll + direct\| < 1e-5 | \|ll - direct\| / \|ll\| | 1.08e-13 | 9.38e-5 (direct ML, mixing logit 0.05 off) | 1e-7 |
| rowSums(P) = 1, tolerance 1e-10 | max\|rowSums(P) - 1\| / eps | 1 | 4.5e6 (P x (1 + 1e-9)) | 16 |
| colMeans, tol 1e-4 | max\|diff\| / se of the column mean | 1.05e-5 | 0.37 (mixing weights 0.01 off) | 1e-3 |
| faithful means, 0.3 and 2 | max\|mean - ref\| / SE | 0.133 | 0.743 (the spherical model, EII) | 0.5 |
| slopes, 0.3 (seed 43) | \|slope - truth\| / SE | 1.22 | 12.2 (x shuffled) | 4 |
| slopes, 0.3, eight models (seed 47) | the same | 2.66 | 11.3 (x shuffled, VVV) | 4 |
| gated >= flat - 1e-6 | (ll_flat - ll_gate) / \|ll_flat\| | 0 | 2.23e-6 (a maximum 1e-3 low) | 1e-8 |
| eight models equal at a point, 1e-10 | max\|diff\| / (\|ll\| eps) | 0 | 1.68e9 (VEI shape 1e-3 off) | 1024 |
| nesting, -1e-6 | (ll_sub - ll_super) / \|ll_sub\| | 0 | 1.46e-6 (1e-3 low) | 1e-8 |
| mclust logLik, 1e-4 | \|diff\| / \|mclust ll\| | 5.64e-16 | 0.113 (mclust EII against VVV) | 1e-9 |
| mclust means, 1e-4 | max\|diff\| / max\|mclust\| | 0 | 2.56e-3 | 1e-5 |
| mclust covariances, 1e-4 | the same | 4.98e-16 | 0.479 | 1e-5 |
| mclust probabilities, 1e-4 | the same | 1.22e-15 | 0.681 | 1e-5 |

Each bound sits between the two, so each test fails on its wrong
answer. Two wrong answers are synthetic, a maximum 1e-3 below the one
it must not be below, because a nested fit that does worse is what
those tests exist to catch. The weakest margin is the faithful means
(0.133 against 0.743 at 0.5); the right answer's ratio is mostly the
rounding of the three-figure references (2.04 against 2.036). The
slope bound of 4 standard errors is a z-bound at a fixed seed; the
largest of 36 values is 2.66. Standard errors are `fixef()`'s
`Est.Error`. Kept as it was: `logLik(f_x) > logLik(f_0) + 20`, a
likelihood-ratio statistic (40 on 4 df), not a tolerance.

The mclust "right" values are near eps because each frmtmb fit starts
at mclust's solution and does not move; the log-likelihood agreement
(two independent evaluations of the same density) is the real check.

### The three core internals

`test-prior-dpar.R`'s `pd_reached()` reads which parameters a prior
reaches through exported functions only: `frm(dry_run = "objective")`
for the spec and frame, `resolve_prior_input()` (exported for
frmtmb.sample) for the resolver, and `par_template()` for the names. No
new export was needed. The two routes give identical answers on 7
priors over lca(), hmm() and mixture_mvn() (`dev/smallfix-pd-reached.R`,
log `dev/smallfix-out/pd-reached.txt`: "agree 7 of 7"). The file's
count is unchanged (13). No `frmtmb:::` is left in frmtmb.latent's
tests (`grep -rn ":::"` finds only the census test's pattern for it).

### The predict(newdata =) refusal

`R/predict-brms.R`: the refusal names the reason for the model at hand.
A family with its own structured simulator: "predict(newdata = ) is not
supported for family '<family>': it draws the whole response at once,
through its own simulator, and a prediction at new rows draws one row
at a time. Drop newdata to predict the rows the model was fitted on"
(raised under the family's package). A residual correlation: "... its
residuals are correlated (ar(time = t, gr = g, cov = TRUE)), and a draw
of them follows the rows the model was fitted on ...". Before and after
on mixture_mvn, hmm, ar(cov = TRUE), unstr() and a plain gaussian
(which answers in both): `dev/smallfix-predict-newdata.R`, logs
`dev/smallfix-out/predict-newdata-{before,after}.txt`.

Tests: latent's `test-mixture-mvn.R` (+2 in the predict block), core's
`test-predict-newdata.R` (a new block, +4, with the guard absent). Both
fail on rellib-r12 (`before.log`: "Expected `conditionMessage(e)` to
match string "is not supported for family 'mixture_mvn""; actual "its
draws are structured (a hidden state sequence, ...").

Not fixed, the same defect: frmtmb.sample's
`posterior_predict(newdata =)` refusal (`R/methods-draws.R:1167`) has
the same three-case text. Wording it per family there needs
`fam_sim_ctx()`, which is not exported, and the brief named core only.

### NEWS

Core "Bug fixes", third bullet; frmtmb.latent "Other changes", both
bullets.

## S3. Two refusal messages

### The Dirichlet prior on a parameter that is not a simplex

brms 2.23.0, run (`dev/smallfix-dirichlet-probe.R`,
`dev/smallfix-dirichlet-store.R`; logs
`dev/smallfix-out/dirichlet-probe-before2.txt`,
`dirichlet-store-before.txt`): `set_prior("dirichlet(2)", ...)` is
accepted for `class = "Intercept", dpar = "sigma"`, `class = "sigma"`,
`"b"`, `"Intercept"`, `"b", dpar = "sigma"` and `"sd"`.
`validate_prior()` and `stancode()` accept each on the model where
brms has that class (`class = "sigma"` on `bf(y ~ x)`, where sigma is
one number; the first probe ran it on `bf(y ~ x, sigma ~ 1)`, where
brms has no class `sigma` and says so, "The following priors do not
correspond to any model parameter"; corrected in fix round 1, review
finding 7). `brm()` stops at the Stan compile for each ("Ill-typed
arguments supplied to function "dirichlet_lpdf": (real, int)").
frmtmb 0.74.0 refused the frmtmb classes in `set_prior()`, and
`class = "sigma"` in `validate_prior()` and `frm()`.

The point chosen: where the prior meets the model, `validate_prior()`
and `frm()` (and `update(prior =)`, `frm(dry_run = "objective")`,
`frm_sample()` and the formula method of `par_template(prior =)`,
which share the resolver), before any fitting. Review finding 8 says
`par_template()` takes no prior; that is its fit method only. The
formula method takes one and refuses this row
(`dev/smallfix-partemplate.R`, log `dev/smallfix-out/partemplate.txt`:
"par_template(formula, prior =) : REFUSED: dirichlet() is the density
of a simplex ..."; "par_template(fit, prior =) : REFUSED:
par_template() has no argument `prior`"). The reasons:

1. Only the model can say whether a name is a simplex. A family can
   declare one (`frmtmb_family(simplex_pars =)`), and `set_prior()` does
   not see the family, so `class = "sigma"` cannot be refused there.
   One point for every spelling must then be the model.
2. `set_prior()` stores `class = "Intercept", dpar = "sigma"` and
   `class = "sigma"` in the same class and dpar fields; only the
   `natural` flag tells them apart. A refusal keyed to the spelling,
   not to the parameter, is what made two points.
3. brms has one point for every spelling too, the compile. frmtmb's is
   earlier and says why.

The cost: `set_prior("dirichlet(2)", class = "b")` is no longer refused
when written, only at `validate_prior()` or `frm()`. Marked BREAKING in
NEWS.

`R/priors.R`: the `set_prior()` block is removed, and
`prior_simplex_check()` repeats the spelling as written ("not class =
"Intercept", dpar = "sigma"", where 0.74.0 said "not class = "sigma""
for both). The link-scale spelling of a declared simplex,
`class = "Intercept", dpar = "tw"`, which `set_prior()` refused before
and the resolver would otherwise read as the simplex, is refused by
name ("names the intercept of a linear predictor, and "tw" is a simplex
this model's family declares, which has none. Write class = "tw"").

Tests: `test-simplex-pars.R`, new block "a Dirichlet on no simplex has
one refusal point, however written": five spellings, each accepted by
`set_prior()` and refused by `validate_prior()` and `frm()` with its
own spelling; the guard absent (a declared simplex fits under
`class = "tw"`); the link-scale refusal. On rellib-r12 the block errors
at its first `set_prior()` ("dirichlet() is the density of a simplex;
...", `before.log`). Two expectations that pinned the write-time
refusal (`test-dirichlet-prior.R`, `test-simplex-pars.R`) now assert
the model-time one. Counts: `test-simplex-pars.R` 37 to 74,
`test-dirichlet-prior.R` 24 to 24.

### print() of a curve given digits =

Accepted. `print.frmtmb_curve()` and `print.frmtmb_feature()` take
`digits = NULL` and pass it to `print.data.frame()` for the rows; the
header keeps its fixed precision. The reasons: a curve is a data frame,
whose print takes `digits`; frmtmb's other prints of a table of numbers
take it (`print.frmtmb_hypothesis()`, `print.frmtmb_multiple()`,
`print.frmtmb_pooled_anova()`); only the prints of a fit and its
summary refuse it, because there it is brms's argument for a summary
that frmtmb renders with fixed rounding. Any other argument is still
refused by name. `frm_curve_feature()`'s result is the same kind of
object, so its print changed with it.

Test: `test-curve.R` "print methods say what was checked" (+6): the
printed rows equal `print(as.data.frame(x)[rows, ], digits =)` at 3 and
7, for a curve and a feature; the two outputs differ; `foo = 1` is
refused. On rellib-r12: "print() has no argument `digits`. It takes: x"
(`before.log`). Count 85 to 91.

### Fix round 1 (review findings 7, 8 and 9)

- Finding 7: the brms probe now fits `class = "sigma"` on `bf(y ~ x)`
  (see above; log `dirichlet-probe-before2.txt`).
- Finding 8: corrected above; the formula method of `par_template()`
  takes `prior =`, the fit method does not.
- A brms prior object with a Dirichlet row on class `b`
  (`dev/smallfix-brmsrow.R`, logs `dev/smallfix-out/brmsrow-{before,
  after}.txt`): before, `validate_prior()` and `frm()` refused it in
  the translation step, "A brms prior table has 1 row with no faithful
  frmtmb spelling: row 1 (dirichlet(2), class = "b"): dirichlet() is
  the density of a simplex; ..."; after, both refuse it in the
  resolver with the common text, "... not class = "b"". A side effect
  of removing the `set_prior()` block; NEWS says so.
- `class = "theta"`: brms 2.23.0's `default_prior()` for
  `mixture(gaussian, gaussian)` lists `dirichlet(1)` on class `theta`,
  and its Stan code declares `simplex[2] theta; // mixing proportions`
  with `lprior += dirichlet_lpdf(theta | con_theta)`
  (`dev/smallfix-theta-probe.R`, log `theta-probe-lane.txt`). frmtmb's
  `theta` is the internal covariance vector. The refusal on
  `class = "theta"` (or `theta1`, ...) now adds "In brms, class =
  "theta" of a mixture() takes a dirichlet() on the mixing
  proportions;" and the advice `brms_prior_class_refusal()` gives for
  the brms route (a proportion with its own predictor is class
  "Intercept" with its dpar). Test: `test-simplex-pars.R` (+1).
- `digits =`: the help of `frm_curve()` and `frm_curve_feature()`
  (their Value sections) says that `print()` passes it to
  `print.data.frame()` for the rows, as significant digits with
  `getOption("digits")` when not given, and that the header keeps its
  own precision. (`print.frmtmb_hypothesis()` rounds to decimal
  places and `print.frmtmb_multiple()` uses `signif()`; the curve
  prints follow the data frame they print.)

### NEWS

Core "Breaking changes" (the Dirichlet, with the brms-object and theta
notes) and "Bug fixes", last bullet (digits).

## S4. The CI skip gate for core

### What the latest runs skipped

Logs downloaded with `gh run view <id> --log` (read only) into
`dev/smallfix-ci-log/` (gitignored). Every "{pkg} cannot be loaded" and
"{pkg} is not installed" skip, per job:

| run | job | skips |
|---|---|---|
| R-CMD-check 40085548721 (8413abf9) and 40088018335 (18b536fc) | windows-latest, macos-latest, ubuntu-24.04 release, ubuntu-24.04 devel | `{evd} is not installed (1): 'test-families-extreme.R:191:3'`, on each |
| test-coverage 40085548802 (8413abf9) | test-coverage (core) | `{evd} is not installed` |
| | coverage frmtmb.sample | `{RTMBode} is not installed` (allowed there) |
| | coverage frmtmb.spatial | `{fmesher} cannot be loaded (1): 'test-spde.R:266:3'` |
| test-coverage 40088018358 (18b536fc) | test-coverage (core) | `{evd} is not installed` |
| | coverage frmtmb.sample | `{RTMBode} is not installed` (allowed) |

No "cannot be loaded" skip on any core job, so no install step to fix
in core. The spatial one was the broken install that 18b536fc fixed
(fmesher's system libraries); in its run the spatial gate passed.

### No allow-list for core (fix round 1, review finding 10)

The first pass allowed `{evd} is not installed` on the four R CMD check
jobs and installed evd in the coverage job. The review showed the skip
was avoidable: extraDistr, already in Suggests, has the same `dgev()`.
The GEV test now reads `extraDistr::dgev(mu =, sigma =, xi =)`
(`dev/smallfix-gev.R`, log `dev/smallfix-out/gev.txt`, the test's grid
of 48 rows, 40 inside the support):

- evd against extraDistr, inside: max |difference| / max(1, |evd|)
  1.169e-15; outside both are -Inf.
- frmtmb's taped density against extraDistr under the test's own
  tolerance (`64 * eps * max(1, |ref|) / min(1, z)`): passes, the worst
  row at 0.027 of its tolerance.
- Two wrong densities fail it: sigma's log-Jacobian dropped (FALSE),
  and xi with its sign flipped (FALSE).

`test-families-extreme.R` passes 108 (fail 0, skip 0), and the GEV
test runs under R CMD check now (core check below: no `{evd}` skip).
`absent: 'evd'` and `any::evd` are gone: no core job allows any
package to be absent.

What the first pass found about evd stays true and explains the old
skip: under `--as-cran`, R 4.6.1's `tools:::.check_packages` sets
`R_check_suggests_only <- TRUE`, and the tests see only declared
packages (`dev/smallfix-ascran-env.R`, log `ascran-env.txt`, deparsed
lines 5922 and 6292-6294), so evd, installed in the user library here,
still skipped in the first core check.

The extensions are removed from the core coverage job on purpose; the
tests that need one skip with their own messages ("not installed:
frmtmb.coupling, ...", "the saved-fit arm needs these packages
installed: ..."), which the gate does not read.

### The changes

- `R-CMD-check.yaml`: a gate step after check-r-package,
  `bash .github/scripts/skip-gate.sh check`, nothing allowed. The check
  directory is check-r-package's default, `check/` under its working
  directory, which is the repository root for core (the logs show
  `check_dir = ("check")` and "using log directory
  '.../frmtmb/frmtmb/check/frmtmb.Rcheck'" on all four).
- `test-coverage.yaml` (core job): a gate step after the Codecov
  upload, as in the extension jobs, nothing allowed.
- `skip-gate.sh`: `mapfile` replaced by a read loop, because the gate
  now runs on macOS, whose `/bin/bash` is 3.2 (the review confirmed
  the runner's `shell: bash` is `/bin/bash`). Fix round 1 (review
  finding 12): an output with no testthat summary line ("[ FAIL n |
  WARN n | SKIP n | PASS n ]"), empty or cut off, fails the gate; it
  read as clean before.

### What was verified locally, and what was not

- Fixtures (`dev/smallfix-gate-fixtures.sh`, log
  `dev/smallfix-out/gate-fixtures.txt`), each with the status it must
  give: no output file (1), an empty file (1), a skip block with no
  summary even with the package allowed (1), a cut-off output (1), a
  clean output (0), a "not installed" skip (1, and 0 with it allowed), a
  "cannot be loaded" skip (1, and 1 with it allowed). "fixtures wrong:
  0".
- The script on the downloaded logs (`dev/smallfix-gate-replay.sh`
  writes each job's "Skipped tests" block and its summary line as a
  testthat.Rout; log `dev/smallfix-out/gate-replay2.txt`), Git Bash
  5.2.37, runs 40088018335 and 40088018358: each R-CMD-check job and the
  core coverage job fail on `{evd} is not installed`, as they must: those
  runs used the evd reference, which this change replaces. The eight
  extension coverage jobs pass (frmtmb.sample with RTMBode allowed, as
  its workflow has it). The first pass's replay of 40085548721 and
  40085548802 is in `gate-replay.txt`; there the 8413abf9 spatial job
  fails on "{fmesher} cannot be loaded" even with fmesher allowed.
- On this lane's four local check directories, after the final core
  check: see "Runs".
- The YAML of the four changed workflows parses with yaml 2.3.12, with
  the values meant: four matrix entries with no `absent`; the gate steps'
  `run` and `shell: bash`; coverage `extra-packages: any::covr,
  any::xml2`; the identity job's four files and paths; pkgdown's new
  step and its three new trigger paths (`dev/smallfix-yaml-parse.R`,
  log `dev/smallfix-out/yaml-parse2.txt`).
- Not verified: anything on GitHub's runners. In particular the gate
  under macOS's bash 3.2 (this machine has none), and the pkgdown step
  on the runner's library layout (see S1, fix round 1).

### NEWS

Core "Other changes".

## Runs

All with `NOT_CRAN=true`, the lane library first. Every `lib:` line
names `frmtmb 0.74.0 C:/Users/adf44/source/r/wt-smallfix-lib` (308
files of this lane's four packages, 134 of the five unchanged
extensions from rellib-r12 over this lane's core).

### One file per process, every suite

`bash dev/smallfix-par.sh full dev/smallfix-out/full.jobs 8` (a copy of
`dev/release/run-par.sh` that writes to `dev/smallfix-out/`; log
`full.log`): 442 of 442 files ran.

| package | files | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 275 | 19137 | 2 | 0 | 201 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 32 | 1838 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 12 | 681 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 430 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 65 | 3024 | 0 | 0 | 9 | 0 |
| frmtmb.spatial | 7 | 328 | 0 | 0 | 2 | 0 |
| frmtmb.survival | 14 | 405 | 0 | 0 | 0 | 0 |

The 2 core failures were this lane's, in `test-message-uniqueness.R`:
`frm_stop(paste(s, collapse = " "))` made the literal " " a message
template that another call has too. Fixed (the message is built into a
variable first) and rerun with the files it reaches (`fix1.log`,
`fix2.log`): `test-message-uniqueness.R` pass 6 fail 0,
`test-extension-home-mixed.R` 166, `test-extension-home.R` 31,
`test-covstruct-seam.R` 131, `test-conditions-census.R` 11, all fail 0.

### Counts of the changed test files, before and after

Before: the base commit's test files on rellib-r12 (`orig.log`). After:
this tree on the lane build (`full.log`).

| file | before | after | why |
|---|---|---|---|
| frmtmb `test-extension-home.R` | 31 | 31 | one expectation's text |
| frmtmb `test-extension-home-mixed.R` | new | 166 | S1 |
| frmtmb `test-predict-newdata.R` | 12 | 16 | S2, +4 |
| frmtmb `test-simplex-pars.R` | 37 | 74 | S3, +37 (five spellings, two routes) |
| frmtmb `test-dirichlet-prior.R` | 24 | 24 | one expectation moved to the model |
| frmtmb `test-curve.R` | 85 | 91 | S3, +6 |
| frmtmb.latent `test-mixture-mvn.R` | 308 | 310 | S2, +2; tolerances rewritten, none dropped |
| frmtmb.latent `test-prior-dpar.R` | 13 | 13 | route only |
| frmtmb.survival `test-saved-fits.R` | 56 | 56 | helper only |
| frmtmb.spatial `test-guardrail.R` | 104 | 104 | helper only |

No expectation was dropped.

### One session per suite, three orders

`bash dev/smallfix-sessions.sh sess 8` (a copy of
`dev/release/run-sessions.sh`; log `sess.log`): 27 of 27 sessions,
fail 0, err 0, warn 0 in every one; 42 LEAK lines, all tagged, none
untagged. Core: files 265/275, pass 19138, 19238, 19332 (forward,
reverse, preload), skip 201, 201, 203. The spread is the 0.74.0
record's (18917, 19017, 19111; `dev/round-20261023.md`). Each
extension's three sessions give the same counts as its
one-file-per-process run.

### The gated tier, on the files the changes reach

`FRMTMB_BRMS_FIT_TESTS=true` (`gated.log`): `test-brms-priors.R` 103,
`test-simplex-brms.R` 15, `test-brms-methods.R` 1006,
`test-brms-agreement.R` 191, `test-mo-factor-re.R` 44; fail 0, skip 0.

### The library-stack tests, three ways

| test | one file per process | R CMD check --as-cran, NOT_CRAN | covr, install_path set |
|---|---|---|---|
| frmtmb `test-extension-home-mixed.R` | 166, skip 0 | in-process blocks run; the child block skips ("the saved-fit arm needs these packages installed: frmtmb.survival, ..."), because --as-cran hides undeclared packages | 166, skip 0 |
| frmtmb.survival `test-saved-fits.R` | 56 | package check fail 0 (the "not installed" arm skips, as at 0.74.0) | 56 |
| frmtmb.spatial `test-guardrail.R` | 104 | package check fail 0 (the same) | 104 |

covr: `dev/smallfix-covr.R` calls `covr::package_coverage()` with
`install_path` and runs the one file inside covr's batch process, where
R_LIBS_SITE is covr's whole stack (`covr-*.txt`).

### R CMD check --as-cran

Built and checked in `dev/smallfix-check/<package>/` (`bash
dev/smallfix-check.sh`). The .Rcheck trees were deleted after copying
`check.txt`, `00check.txt` and `testthat.Rout.txt`.

| package | Status | tests |
|---|---|---|
| frmtmb | 1 NOTE (HTML manual: V8 unavailable) | FAIL 0, WARN 0, SKIP 229, PASS 18953 |
| frmtmb.latent | 1 NOTE (the same) | FAIL 0, WARN 0, SKIP 4, PASS 673 |
| frmtmb.survival | 1 NOTE (the same) | FAIL 0, WARN 0, SKIP 4, PASS 366 |
| frmtmb.spatial | OK | FAIL 0, WARN 0, SKIP 3, PASS 282 |

The core check also confirms S4's argument on this machine: evd is in
the user library, and the GEV test still skips as "{evd} is not
installed" under --as-cran. The gate on the four local check
directories (`gate-local-checks.txt`): core passes with evd allowed and
fails without it; the three extensions pass with nothing allowed.

### OpenBLAS

The latent tests are exactness tests whose tolerances changed, so they
ran under OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh
0.3.26 lapack`: 616 LAPACK and 74 BLAS routines routed; 4 threads):
`test-mixture-mvn.R` pass 310, fail 0 (`openblas-mvn.txt`). The ratio
measurements under the same R (`mvn-ratios-openblas.txt`) separate
right from wrong as on reference BLAS; only T1 (1.07e-13), T9
(4.51e-16) and the probabilities (1.33e-15) moved, far below their
bounds. The OpenBLAS copy of R was deleted.

### Fix round 1 runs

The tables above are the first pass. After fix round 1, all with the
lane library first; every `lib:` line names this lane's frmtmb.

- Every suite, one file per process (`full2.log`): 442 of 442 files,
  fail 0, err 0, warn 0 in every package. Core 275 files, pass 19229,
  skip 201; coupling 551, eam 1838, latent 681, learn 430, ode 549,
  sample 3024, spatial 328, survival 405 (skips as in the first pass).
  Core is +92 on the first pass: `test-extension-home-mixed.R` 166 to
  256 (+90: the six methods in the in-process list in three blocks,
  four in the child block, and three new blocks), `test-simplex-pars.R`
  74 to 75 (theta), the two failures of `test-message-uniqueness.R`
  gone.
- One session, three orders, core, survival and spatial (`sess2.log`):
  9 of 9; core pass 19229, 19329, 19423, skip 201, 201, 203;
  survival 405 and spatial 328 in each order; fail 0, err 0, warn 0;
  21 LEAK lines, all tagged.
- The library-stack test, `test-extension-home-mixed.R`: one file per
  process 256 (skip 0); covr with `install_path` 256 (skip 0,
  `covr-frmtmb-r1.txt`); R CMD check --as-cran with NOT_CRAN: runs, and
  only the real-fit block skips ("the saved-fit arm needs these
  packages installed: frmtmb.survival, ..."), as designed, since core
  declares no extension; the pkgdown step's script, locally: 256, skip
  0.
- R CMD check --as-cran, core (`dev/smallfix-check/frmtmb/`):
  `Status: 1 NOTE` (the V8 HTML-manual note); tests FAIL 0, WARN 0,
  SKIP 228, PASS 19011. No `{evd}` skip any more (229 to 228), and the
  gate passes on this check directory with nothing allowed.
- OpenBLAS 0.3.26 with LAPACK, `test-families-extreme.R` with the
  extraDistr reference: pass 108, fail 0 (`openblas-extreme.txt`); the
  OpenBLAS copy of R was deleted again.
- Workflows parsed (`yaml-parse2.txt`): `R-CMD-check.yaml`,
  `test-coverage.yaml`, `absent-helper-identical.yaml`, `pkgdown.yaml`.
- Probes rerun on this build: the methods probe (`methods-probe-*.txt`),
  the handler scan (0 left), the theta probe (`theta-probe-lane.txt`),
  the brms Dirichlet probe with `class = "sigma"` on `bf(y ~ x)`
  (`dirichlet-probe-before2.txt`: brms accepts all six spellings in
  `set_prior()`, `validate_prior()` and `stancode()`, and every one
  fails the Stan compile; `dirichlet-probe-after2.txt`: frmtmb accepts
  all six in `set_prior()` and refuses all six in `validate_prior()`
  and `frm()`).

## Bumps and floors

- frmtmb: minor. Behavior changes: the Dirichlet refusal point
  (BREAKING), the refusal texts, `print(curve, digits =)`.
- frmtmb.latent: patch (tests and NEWS only). `test-mixture-mvn.R` now
  pins the `predict(newdata =)` text, so its tests need the next
  frmtmb; its code does not. A Depends floor raise is for the
  consolidation.
- frmtmb.survival, frmtmb.spatial: no bump needed (a test helper only,
  which reads no new core function).
- The ported brms ledger: no row of `dev/brmsport-verdicts.tsv` names
  dirichlet or digits. I expect no row to flip, but did not regenerate
  the ledger.

## Overlaps with other lanes

Read from their worktrees (read only) on 2026-10-24:

- sarfcor changes `R/extension-home.R`: it adds, at the end of
  `require_structure_packages()`, a loop that refuses a registered
  residual correlation term (`sar()`) whose package is absent
  (`ac_block_entry()`). This lane rewrites that function, which now
  returns early when no structure part is missing, so a per-file merge
  conflicts. Resolve by keeping this lane's body and running sarfcor's
  loop at its top, before the early return (it must run when no
  block is missing); better, make the term a part in
  `structure_missing_parts()`,
  so its refusal joins the combined message.
- brmsnames changes `R/priors.R` (`set_prior_one()` near line 712,
  `family_simplex_targets()` and `prior_simplex_check()` near 4012 to
  4064) and `R/brms-shapes.R` (near line 441). This lane changes
  `set_prior_one()` at 709 to 723, `prior_simplex_check()` at 4010 to
  4060, and `fit_draw_space()` at 424. The `R/priors.R` hunks likely
  touch the same lines; merge `prior_simplex_check()` by hand.
- `NEWS.md`: all three other lanes add to the development section.
- `.github/workflows`, frmtmb.latent: no other lane changes them.
  `pkgdown.yaml` has an uncommitted edit in main's working tree; this
  lane's hunks merge with it cleanly (S1, fix round 1).
- Corrected after the short check: sarfcor also changes
  `R/predict-brms.R` and `R/priors.R`, and both merge cleanly with this
  lane's. sarfcor's `require_structure_packages()` now ends with
  `stop_autocor_missing(autocor_missing_parts(frame))`, so at
  consolidation its autocor parts join `structure_missing_parts()` and
  the combined message, in place of the loop the first entry above
  describes. `dev/smallfix-handlers.R` must be rerun on the merged
  tree: the other lanes' handlers are not in this lane's scan.
- Fix round 1 widens S1: the handler rewrite and the head checks touch
  25 files of `R/`. Of those, sarfcor also changes `R/autoscale.R`,
  `R/families.R`, `R/fit.R`, `R/predict.R`, `R/predict-brms.R` and
  `R/priors.R`, and brmsnames `R/brms-shapes.R` and `R/priors.R`.
  Each S1 hunk there is one
  handler line (`error = keep_missing(...)`) or one inserted
  `require_fit_packages()` line, so a conflict is resolved by keeping
  the other lane's line and making the same one-word change to it. A
  new handler another lane adds that can reach the refusal should use
  `keep_missing()`; `dev/smallfix-handlers.R` finds any that do not.

## Not done, and for the user

- frmtmb.sample's `posterior_predict(newdata =)` has the same stale
  three-case text (S2); not changed.
- The gate on GitHub's runners is unverified (S4). The first push will
  show it; watch macOS bash and the Windows `find`.
- `dev/smallfix-ci-log/` (96 MB of downloaded logs, gitignored) is kept
  as the source of S4's table; delete it when the round closes.
