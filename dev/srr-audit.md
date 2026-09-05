# srr tag audit before the rOpenSci submission

Every `@srrstats` and `@srrstatsNA` tag in `R/`, `tests/testthat/` and
`vignettes/` was read against the code or prose it sits on, and given a
verdict. The package has been through 39 versions since most of the
tags were written, so the audit looks for drift, not for absence.

**Reading the locations.** Every file and line below was recorded at
the audit run. The line numbers have since moved, most of them by
tens or hundreds of lines, and one row changed file. The standard name
is the durable pointer: find a row's tag with
`grep -rn "@srrstats {STANDARD}" R/ tests/testthat/ vignettes/`. Every
tag in this audit is still present in the tree. Note also that the
sampling surface, the ODE seam and the two latent-state families left
the core package at 0.47.0 and 0.48.0, so a note below that names
`frm_sample()`, `check_laplace()`, `frm_ode()`, `hmm()` or `lca()` is
describing something that now lives under `extensions/`. The rows
where that changes the meaning are corrected in place.

Verdicts:

- **holds**: the claim is still literally accurate.
- **drifted**: the claim names something renamed, moved, gone, or
  counted wrong.
- **weakened**: the claim is broadly true but overstates or
  under-describes what the code does.

Scope: 118 tag entries. `srr` reports 126 standards covered (116
complied with, 10 N/A); the entry count differs because one tag can
claim several standards and `RE4.6` is tagged in two places.

## Totals

| Verdict | Count | Fixed in place | Logged |
| --- | --- | --- | --- |
| holds | 88 | n/a | n/a |
| drifted | 12 | 12 | 0 |
| weakened | 18 | 18 | 6 follow-ups |

Every non-holding tag now reads true. Six of them also left a
follow-up that is code, not prose; those are in "Logged" at the end.

## Verdict table

### `R/` general standards

| standard | location | verdict | note |
| --- | --- | --- | --- |
| G1.0 | `R/frmtmb-package.R:117` | weakened | Claimed all references live in the package-level `@references`. Four methods carry theirs on their own help page. Sentence added naming them. A separate gap is logged. |
| G1.1 | `R/frmtmb-package.R:127` | holds | "Algorithm provenance" says what the tag says it says. |
| G1.2 | `R/frmtmb-package.R:134` | drifted | Claimed `CONTRIBUTING.md` carries a life-cycle statement. It does not. Claim narrowed to the README. |
| G1.4 | `R/frmtmb-package.R:138` | drifted | `NAMESPACE` counts were 72/100. Actual: 134 `export()`, 232 `S3method()`. Updated and anchored to this version. |
| G1.4a | `R/frmtmb-package.R:143` | drifted | "about 290" internal functions; actual is above 500. Ten internal functions had no roxygen block at all, so "All ... are documented" was false. Both fixed: the count is corrected and the ten blocks were written. |
| G2.0 | `R/bf.R:38` | weakened | Named `se()` and `spde()` as routing tuning arguments through the scalar validator. Neither does. Replaced with `mm()` and `gr()`, which do. |
| G2.1 | `R/bf.R:44` | holds | Formula check, finite-numeric tuning check, and the character/factor multiplier refusal all verified. |
| G2.14c | `R/bf.R:50` | drifted | Claimed `mi(x, sdx)`. That spelling errors; the measurement-error form is `bf(x \| mi(sdx) ~ ...)`. Corrected. |
| G2.14c | `R/multiple.R:64` | holds | Rubin pooling, pooled `VarCorr`, `anova` and `hypothesis` methods all present. |
| G2.1 | `R/multiple.R:60` | holds | List-of-at-least-two assertion verified. |
| G2.0, G2.1 | `R/families.R:95` | holds | `stopifnot()` contract and `as_frmtmb_family()` dispatch verified. |
| G2.0, G2.1 | `R/priors.R:68` | holds | Length-one character assertion and per-distribution arity checks verified. |
| G2.3a | `R/priors.R:74` | holds | `match.arg()` over the five prior classes; LKJ checked in both directions. |
| G2.2 | `R/predict.R:2106` | holds | `single_response()` (nee `uni_resp()`) guards all six named methods and four more. |
| G2.3, G2.3a | `R/predict.R:927` | holds | `match.arg(type)`; unknown `resp` and `dpar` both list what is available. |
| G3.0 | `R/predict.R:932` | weakened | "Floating-point values are never compared for equality" is false as written: the zero- and one-inflation indicators and several degenerate-input guards use exact equality. Rewritten to claim the true property (no *computed* quantity is compared for equality) and to name the exceptions and why they are correct. |
| G2.4, G2.4a-c, G2.4e | `R/fit.R:89` | weakened | Said grouping level labels come from `as.character()`. They come from `levels()`; `as.character()` is applied to grouping *values* before matching. Reworded to what the code does. |
| G2.4d | `R/frame.R` (moved) | drifted | Was an `@srrstatsNA` claiming the package never converts an input to `factor`. It does, in two places (`R/frame.R` grouping levels, `R/sandwich.R` `cluster`), and one of them has a test. The standard is **met**, not N/A. The block moved out of `NA_standards` and is now tagged at `report_datetime_columns()` in `R/frame.R`, naming both sites. |
| G2.5 | `R/fit.R:98` | holds | `mo()` ordered-factor error text matches the source verbatim; unordered-ordinal warning and factor-response refusal verified. |
| G2.5 | `R/frame.R:816` | holds | `report_datetime_columns()` names column, unit and origin. |
| G2.6 | `R/fit.R:108` | holds | Every one-dimensional response reaches the objective as `as.numeric(as.vector(y))`. |
| G2.8 | `R/fit.R:112` | holds | `frmtmb_spec` and `frmtmb_frame` are the only two classes downstream; both reachable through `dry_run`. |
| G2.9 | `R/frame.R:822` | holds | The lossy-conversion diagnostic still fires on every assembly. |
| G2.13 | `R/fit.R:253` (was `:118`) | weakened | "assembly errors if any missing value remains" is false for `mi()` columns and for the response of a family that declares `keep_na`, which are exempt by design. The exemption is now stated. The tag is keyed on the generic `keep_na` flag, not on a family name; `hmm()` is one consumer of that flag and now lives in `frmtmb.latent`. |
| G2.14, G2.14a | `R/fit.R:123` | weakened | Same carve-out, plus a real behavior gap: `na.action = stats::na.fail` is silently downgraded on those models. Doc corrected; the behavior is logged for the input-validation lane. |
| G2.14b | `R/fit.R:129` | holds | One `message()` per fit; `na.action()` on the fit names the rows. |
| G2.15 | `R/fit.R:276` (was `:135`) | weakened | The "no defensive `na.rm`" claim is contradicted by the `mi()` and `keep_na` readers, which use `na.rm = TRUE` deliberately. The exception is now stated and justified. The reader is generic over the `keep_na` flag; there is no family-specific reader in core. |
| G2.16 | `R/fit.R:141` | holds | The response check is literally `any(!is.finite(y) & !is.na(y))`. |
| G5.0 | `R/frmtmb-package.R:147` | holds | All six named data sets are in live use. |
| G5.2a | `R/srr-stats-standards.R:48` | weakened | Opened with "Every condition message", but the described AST walk covers `stop()` only. Narrowed to `stop()`. The count is untouched: an independent recount at the end of this lane gives **651 literal `stop()` texts, 651 distinct, 0 duplicates**. |
| G5.10 | `R/frmtmb-package.R:153` | holds | All three environment variables present; 31 files use `skip_on_cran()`. |
| G5.12 | `R/frmtmb-package.R:160` | holds | `CONTRIBUTING.md` documents the variables, the requirements, and that nothing is downloaded. |
| RE1.4 | `R/fit.R:205` (was `R/interop.R:1478`) | holds | The tag names the two approximations, the regimes where each degrades, the three remedies inside this package, and the direct measurement in the companion sampling package. The measurement is `frmtmb.sample::check_laplace()`, which left core at 0.47.0; the tag points at the package rather than naming the function's return fields. The quadrature remedy is tested against `lme4::glmer(nAGQ = 25)` and GLMMadaptive. |
| RE4.13 | `R/interop.R:324` (was `R/interop.R:1833`) | holds | `getME()` Z dimnames come from the stored model frame's row names. |

### `R/` regression standards

| standard | location | verdict | note |
| --- | --- | --- | --- |
| RE1.0 | `R/fit.R:85` | holds | Formula interface only; no matrix entry point. |
| RE1.3 | `R/predict.R:937` | weakened | Claimed the named test asserts one naming scheme "across" `vcov()`, `confint()` and `fixef()`. It asserts the first two; `fixef()` is not compared anywhere. Claim narrowed to what is asserted; adding the missing assertion is logged. |
| RE2.0 | `R/fit.R:1014` | holds | `@param autoscale` states the qualification and the exclusions. |
| RE2.1 | `R/fit.R:145` | holds | `na.action` is a formal of `frm()`. |
| RE2.2 | `R/bf.R:58` | holds | The `na.pass` `mi()` branch is where the tag says. |
| RE2.3 | `R/fit.R:1024` | holds | Centering happens only when an intercept exists. |
| RE2.4, RE2.4a | `R/fit.R:148` | holds | QR rank drop with the stored null space. |
| RE2.4b | `R/confint.R:875` | weakened | Separation is found only when `diagnose()` is called, not at fit time, and only for binomial-type fits. Both limits now stated, with the reason a zero-residual continuous fit is not the same thing. |
| RE3.0 | `R/fit.R:156` | weakened | Two of the four warnings cannot fire under the default `se = FALSE`, and they are an `if`/`else if` pair, so "one per diagnostic" overstated. Rewritten to say when each fires and that the last two are exclusive. |
| RE3.1 | `R/fit.R:161` | holds | `warning()` conditions; the `diagnose()` fields match. |
| RE3.2 | `R/fit.R:1029` | holds | Documented defaults verified against `frmtmb_control()`. |
| RE3.3 | `R/fit.R:1034` | holds | All four thresholds are formals. |
| RE4.0 | `R/fit.R:167` | holds | Class plus the full `stats` method surface. |
| RE4.1 | `R/fit.R:170` | holds | `dry_run` returns spec and frame without fitting. |
| RE4.2 | `R/methods-fit.R:440` | holds | `coef()` is the fixef broadcast plus modes, with the no-random-effect fallback. |
| RE4.3 | `R/confint.R:333` | holds | Four methods; row-name identity asserted in `test-methods-audit.R`. |
| RE4.4 | `R/fit.R:175` | holds | `call` stored; `formula()` method registered. |
| RE4.5 | `R/fit.R:177` | holds | `nobs()` and `na.action()` both registered and tested. |
| RE4.6 | `R/methods-fit.R:340` | holds | Matches `vcov.frmtmb_fit`, including the non-finite guard. |
| RE4.6 | `R/sandwich.R:410` | holds | The cluster-robust claim uses `vcov(full = TRUE)` as the bread; consistent with the other RE4.6 claim, no contradiction. |
| RE4.7 | `R/confint.R:881` | holds | Every field named is in the `diagnose()` return list; all three accessors exported. |
| RE4.8 | `R/fit.R:179` | drifted | Named `offset` as an addition term in `fit$frame$aterm_values`. It is not one, and `trunc` is stored split. A reviewer can disprove this in one line. Corrected to the real set, with the offset placed on the linear predictor where it lives. |
| RE4.9 | `R/predict.R:943` | weakened | "asserted to equal `fitted()` for every family, with no exception" is contradicted by the package's own tests: `cox()` refuses both, and the fuzz invariant skips ordinal families. Rewritten to the true scope, naming the `cox()` case. |
| RE4.10 | `R/predict.R:2111` | holds | Four residual types; deviance section and DHARMa hook present. |
| RE4.11 | `R/confint.R:1145` | holds | `logLik()` carries `df` and `nobs`; the AIC family is registered. |
| RE4.12 | `R/families.R:102` | holds | Link, inverse link and `mu.eta` all present; insight methods registered. |
| RE4.14 | `R/predict.R:951` | holds | New-level marginal variance and GP kriging variance both returned. |
| RE4.16 | `R/predict.R:957` | holds | New groups accepted; the refusal names the offending levels. |
| RE4.17 | `R/fit.R:185` | holds | `print()` gives formula, family, fixef and VarCorr. |
| RE4.18 | `R/fit.R:188` | holds | `summary()` is genuinely distinct: z and p values, `ngrps`, BIC, rescor, autocor, smooth EDF, with a memoized `sdreport`. |
| RE6.0 | `R/conditional-effects.R:1345` | holds | `plot.frmtmb_fit` exists. |
| RE6.1 | `R/conditional-effects.R:1355` | holds | Registered as an S3 method in `NAMESPACE`. |
| RE6.2 | `R/conditional-effects.R:1350` | holds | Panel one is fitted values against Pearson residuals. |

### `tests/testthat/`

| standard | location | verdict | note |
| --- | --- | --- | --- |
| G2.7, G2.10, G2.11, G2.12 | `test-tabular-inputs.R:5` | holds | Tibble, data.table, attribute-carrying and list-column cases each have a matching test in the same file. |
| G5.2, G5.2b | `test-open-issues.R:6` | drifted | Counts were 193/16/8. Actual: 552 `expect_error()`, 28 `expect_warning()`, 26 `expect_message()`. Also "essentially every one" understated: all but one `expect_error()` carries a message regexp. Corrected. |
| G5.3 | `test-quadrature-defects.R:13` | holds | The `anyNA` sweep and every cross-referenced finiteness assertion exist. |
| G5.4 | `helper-reference.R:3` | weakened | "Every model class the package supports is checked against a package that implements the same likelihood" contradicted the very next tag, G5.4a, which covers the classes with no reference. Scoped to "every model class for which an existing implementation is available". |
| G5.4a | `helper-reference.R:17` | holds | 19 test files build a hand-written likelihood with `RTMB::MakeADFun()`. |
| G5.4b | `helper-reference.R:24` | holds | References are called live, never transcribed. |
| G5.5 | `helper-reference.R:28` | drifted | "55 of the 61 test files" is stale; it is 82 of 95. The qualitative tail still holds. Corrected. |
| G5.6, G5.6a | `helper-reference.R:35` | holds | A recovery test found for every listed area; tolerances at the call site. |
| G5.6b | `test-fuzz.R:17` | holds | The seed derivation, the truth constant and the binomial tail are all literally as described. |
| G5.7 | `test-perf.R:1` | weakened | The scaling assertion is exact, but the tag pointed at "the tape canary above" and the canary is below it. Changed to "in this file". |
| G5.8 | `test-edgecases.R:4` | holds | Every listed edge condition maps to a test in the file. |
| G5.8a | `test-edgecases.R:21` | holds | Both error messages verified verbatim against the source. |
| G5.8b | `test-edgecases.R:12` | weakened | Three of the four refusals live in `test-open-issues.R`, so a reviewer following the tag to its file finds only one. Pointer added. |
| G5.8c | `test-edgecases.R:26` | holds | All-`NA` and all-identical columns, including the boundary assertion. |
| G5.8d | `test-edgecases.R:32` | holds | p = 20, n = 12. |
| G5.9, G5.9a | `test-edgecases.R:37` | holds | Including the proportionality half. |
| G5.9b | `test-fuzz.R:28` | drifted | Two of the seven invariant names do not exist: `permutation` is `row_permutation`, and `vcov_summary` was never an invariant (the real ones are `vcov_dim`, `vcov_psd`, `summary_prints`). Corrected. |
| RE7.0, RE7.0a | `test-edgecases.R:62` | holds | The refusal message names the offending columns. |
| RE7.1, RE7.1a | `test-edgecases.R:44` | holds | 12 seeded pairs, timed at n = 2000. |
| RE7.2 | `test-edgecases.R:57` | holds | All five accessors plus the `na.action` case. |
| RE7.3 | `test-methods-audit.R:4` | weakened | Claimed `coef()` is exercised in this file. It is not; the `lme4::coef()` comparison is in `test-sugar.R`. The tag now attributes it, as it already does for the other borrowed assertions. |

### `vignettes/inputs.Rmd`

| standard | location | verdict | note |
| --- | --- | --- | --- |
| G1.3 | `inputs.Rmd:543` | holds | "Terminology" section present. |
| G1.5 | `inputs.Rmd:548` | holds | All eight scripts in the reproduction table exist in `dev/`. |
| G1.6 | `inputs.Rmd:551` | holds | The alternative-implementation comparison script exists. |
| G2.0a | `inputs.Rmd:555` | drifted | The argument table omitted `data2` and gave `dry_run` three values instead of four. Both fixed. |
| G2.1a | `inputs.Rmd:559` | weakened | Same omission; resolved by the `data2` row. |
| G2.3b | `inputs.Rmd:562` | drifted | The character-options table omitted five exported `match.arg()` arguments and gave `set_prior(class =)` four values instead of five. All six rows fixed. |
| G2.5 | `inputs.Rmd:567` | holds | "Dates and times" section covers unit, origin and centering. |
| RE1.1 | `inputs.Rmd:574` | holds | The stage diagram names each package in the chain. |
| RE1.2 | `inputs.Rmd:578` | holds | Accepted and rejected classes both listed. |
| RE1.3a | `inputs.Rmd:581` | holds | Preserved and not-preserved lists match the claim item for item. |
| RE5.0 | `inputs.Rmd:588` | holds | Stage table plus the n = 180 to n = 100000 scaling table. |

### `NA_standards`

| standard | verdict | note |
| --- | --- | --- |
| G3.1 | holds | `stats::cov()` appears nowhere in `R/` except inside this tag's own prose. |
| G3.1a | holds | Follows from G3.1. |
| G4.0 | holds | No `write*()`, `saveRDS()`, `sink()`, `save()` or `cat(file =)` call in `R/`; the three `saveRDS` hits are all comment or message text. No exported function takes a path. |
| G5.1 | holds | All six named data sets in use; `tests/` is not in `.Rbuildignore`, so `helper-reference.R` ships. |
| G5.4c | holds | Consistent with G5.4b: references are installable packages, called live. |
| G5.11 | weakened | "They simulate their own data from a seed" under-described the `NOT_CRAN` tier, which also uses data sets from `Suggests`. `CONTRIBUTING.md` already said this correctly. Tag brought into line. |
| G5.11a | holds | Follows from G5.11. |
| RE4.15 | holds | The `predict(se.fit = TRUE)` sentence is still literally true of `R/predict.R`. |
| RE6.3 | holds | Non-estimable rows return `NA` with a warning. |
| RE7.4 | holds | Follows from RE4.15. |
| G2.4d | drifted | Removed from this list; the standard is met, not N/A. See the `R/` table. |

## Logged: work this audit did not do

These are code, not prose, and each is larger than a tag edit.

1. **Primary references for three surfaces.** `G1.0` is satisfied for
   the package as a whole, but these implemented methods cite no
   primary reference anywhere: the LKJ correlation prior (Lewandowski,
   Kurowicka and Joe 2009), `mo()` (Buerkner and Charpentier 2020) and
   `tweedie()`. Add them to the package-level `@references` or to each
   function's own page. Two more surfaces were on this list and have
   since left the core package: `hmm()` (the forward algorithm) is in
   `frmtmb.latent`, and the NUTS path behind `frm_sample()` (Hoffman
   and Gelman 2014; Monnahan and Kristensen 2018 for tmbstan) is in
   `frmtmb.sample`. Each is now that package's gap to close, and
   neither can be closed from core's `@references`.

2. **`na.action = stats::na.fail` is silently downgraded.** In
   `R/frame.R`, a model with an `mi()` addition term or a family whose
   structure declares `keep_na` replaces the user's `na.action` with
   `stats::na.pass` for the whole model frame, re-implements row
   dropping over the non-`mi()` columns, and stamps the result
   `"omit"` unless the action was identically `na.exclude`. So
   `frm(..., na.action = stats::na.fail)` on such a model omits rows
   instead of erroring. Apply the requested action to the non-`mi()`
   subframe, or refuse explicitly. This is `G2.14a`, and the
   input-validation lane owns the file.

3. **RE1.3: assert the `fixef()` half.** `test-methods-audit.R` asserts
   the `vcov()`/`confint()` naming identity but never compares
   `names(fixef(fit))` with `rownames(vcov(fit))`. One
   `expect_identical()` would let the original, stronger wording come
   back.

4. **RE2.4b: separation is diagnosed on demand only.** Nothing reports
   complete separation at fit time; `check_convergence()` has no
   separation check. Consider raising it as a fit-time warning, which
   is what the standard reads most naturally as asking for.

5. ~~**`vignette("brms-migration")` is stale in one line.**~~
   **RESOLVED.** The complaint was that the "When you still want brms"
   section listed "`loo`-based model comparison" among the reasons to
   go back to brms, after `loo()` and `waic()` began working on
   `frm_sample()` draws. That section now opens by saying that with
   `frmtmb.sample` installed a posterior, brms's default priors and
   `loo()`-based comparison are all available, and it lists a narrower
   set of reasons. One line in that narrower list has since gone stale
   in turn: "discrete latent structure beyond observation-level
   mixtures" is now in `frmtmb.latent`, so it is no longer a reason to
   reach for brms.

6. **The counts in these tags will drift again.** Four hand-maintained
   counts had gone stale, two by a factor of three, and the `G5.2a`
   paragraph says out loud that it depends on a manual re-run. Either
   add a mechanical check over the tag text, or drop the numbers from
   the tags that do not need them. `G5.2a` does need its number.

## Verification

- `srr::srr_stats_pre_submit()`: **All applicable standards [v0.2.0]
  have been documented in this package (116 complied with; 10 N/A
  standards).**
- `srr::srr_report()` builds clean. Rendered copy checked in at
  `dev/srr-report.html`. It reports G 61/68 and RE 45/48 tagged, plus
  G 7 and RE 3 N/A, which closes both categories.
- `roxygen2::roxygenise()` after the audit leaves `man/`, `NAMESPACE`
  and `DESCRIPTION` unchanged, and adds no new roxygen warning.
- `G5.2a` recount by AST walk over `R/`: 651 `stop()` calls with
  literal text, 651 distinct texts, 0 duplicates. Unchanged by this
  lane, so the number in `R/srr-stats-standards.R` stays at 651.
