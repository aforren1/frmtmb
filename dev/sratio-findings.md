# Lane sratio: sratio()'s thresholds are unordered, as in brms

Round of 2026-09-25. Worktree `/home/user/wt/sratio`, branch
`lane/sratio` (the consolidated branch at a897706d), private library
`/opt/rlib/lane-sratio`.

## The defect

Found by lane thres (`dev/thres-findings.md`, "Defects seen, not
fixed"). brms 2.23.0 does not order sratio's thresholds:
`brms:::has_ordered_thres()` is TRUE for `cumulative()` and FALSE for
`sratio()` and `cratio()`, and brms's Stan program for sratio declares
`vector[nthres] Intercept` (read off `make_stancode()` here). frmtmb's
`fam_sratio()` held them as (first threshold, log increments), like
`cumulative()`. Where the unconstrained optimum has crossing
thresholds, frmtmb's estimate sat on its ordering boundary, with two
thresholds equal to about 1e-8, and differed from brms's mode without
a word. A stopping-ratio category probability is a product of hazards,
`F(tau_k - eta) * prod_{j<k} (1 - F(tau_j - eta))`, so it is positive
whatever order the thresholds take; nothing in the model needs the
constraint.

The defect was not rare. Three of the 19 sratio fits built from the
existing tests' models sat on the boundary (below), one of them brms's
own `inhaler` example with `cs(treat)`.

## What changed

`sratio()` now holds its thresholds as an unconstrained vector, exactly
as `cratio()` does. `tau_raw` of an sratio fit is the thresholds
themselves.

| file | change |
|---|---|
| `R/families.R` | `fam_sratio()`: the log-density reads `tau_raw` directly (the loop that built increasing thresholds from log increments is gone); `thres_finalizer()`, `ord_tau_init()`, `ord_sim()` and `post$ord_thresholds` get `ordered = FALSE`. Comments on `ord_threshold_map()` and `ord_tau_from_raw()`. The exported family help gains one paragraph: cumulative's thresholds are increasing, the other three families' are unconstrained, as in brms |
| `R/thres.R` | comment on `thres_tau()`; the grouped density, simulator, start values and map already follow the `ordered` flag the family passes |
| `R/priors.R` | class `"Intercept"` on an ordinal family: `ordered` is `identical(family, "cumulative")`, so sratio's entry is scale `"internal"` like cratio's, with no log-Jacobian. The wrong comment and the `@section Ordinal thresholds` text of `?set_prior` are corrected |
| `R/compat.R` | the `thres()` x `prior` row said "each level an ordered vector of its own"; it now says ordered for cumulative() and unconstrained for sratio(), cratio() and acat() |
| `man/set_prior.Rd`, `man/frmtmb-families.Rd` | roxygen 8.1.0 |
| `vignettes/case-studies.Rmd` | "`tau_raw` holds the thresholds on an internal increasing scale" was false for the sratio fit there; one sentence |
| `extensions/frmtmb.sample/R/methods-draws.R` | comment only |
| `dev/thres-findings.md` | its "Defects seen, not fixed" entry says the defect is fixed here |
| `tests/testthat/test-sratio-thresholds.R` (new), `extensions/frmtmb.sample/tests/testthat/test-sratio-draws.R` (new) | see Tests |
| `tests/testthat/helper-brms.R`, `test-brms-likelihood.R`, `test-brms-shapes.R`, `test-diagnostics-ux.R`, `test-review-v29.R`, `test-simulate-density.R`, `test-thres.R` | these hard-coded sratio as ordered; see below |
| `dev/sratio-invariance.R`, `dev/sratio-brms-stan.R`, `dev/sratio-saved-fit.R`, their logs, `dev/sratio-invariance-before.csv`, `dev/sratio-seen-failing-base.txt` | evidence |

### Every consumer of `tau_raw`, followed

`grep -rn "sratio\|ordered\|tau_raw\|ord_tau_from_raw\|ordthres"` over
`R/` and every extension's `R/`:

- **Family** (`fam_sratio()`): density, start values, threshold map,
  category probabilities (`ord_cat_probs()` takes the thresholds, not
  `tau_raw`) and simulator. Changed.
- **Grouped thresholds** (`R/thres.R`): `thres_tau()`,
  `thres_finalizer()`'s `ordered`, the grouped start values and map all
  take the flag from `fam_sratio()`; no family-name test exists there.
- **Priors**: the one family-name test (`c("cumulative", "sratio")`)
  was in `ordinal_threshold_entry()`. Changed. `"ordthres"` is now
  reached by cumulative alone; `draw_prior_entry()`'s refusal of an
  `"ordthres"` draw and `prior_logdens()`'s ordered branch need no
  change.
- **Read through the family's map, unchanged**: `fixef()`, `vcov()`,
  `summary()`, `variables()`, `hypothesis()` (`brms_fixef_rows()`,
  `hyp_put_ordinal()` via `ord_threshold_values()`), `coef()`,
  `fitted()`, `predict()`, the delta-method standard errors and
  `conditional_effects()` (all through the family's own log-density,
  `ord_probs_from_eta()`), `ordinal_ncat()`, `confint()` (its internal
  `tau_raw_k` rows are now the thresholds), the `marginaleffects` and
  `insight` coefficient vectors, and frmtmb.sample's
  `draws_fixef_ordinal()`, `posterior_epred()` and `posterior_predict()`.
- **Not about storage**: `coef_shift_thresholds()` keys brms's SIGN on
  `thres_minus_eta` (cumulative, sratio), which is still right; `cs()`
  counts `length(tau_raw)`; the multivariate per-response extras
  (`extra_tpl_name()`) name the block and do not read it; emmeans works
  on the latent predictor, which no threshold enters.

### Tests that hard-coded the old storage

`helper-brms.R`'s `brms_ord_thresholds()`, which feeds the gated
brms-likelihood row 12, and the reference computations in
`test-diagnostics-ux.R`, `test-thres.R` (two), `test-review-v29.R` and
`test-brms-shapes.R` read sratio's `tau_raw` as (first, log
increments). Each now reads it as the thresholds. `test-simulate-density.R`
passes `tau_raw` straight through, so its sratio entry changes from
`c(-1, 0, 0)` to `c(-1, 0, 1)`, the same thresholds as before.
`test-brms-shapes.R` asserted that sratio's reported thresholds DIFFER
from `tau_raw`; it now asserts they are the family map of `tau_raw`,
which is the identity for sratio.

## Validation

### Where the thresholds do not cross, the fit is unchanged

`dev/sratio-invariance.R`, logs `dev/sratio-invariance-before-log.txt`
and `dev/sratio-invariance-after-log.txt`. The 19 sratio models of the
existing tests, rebuilt with the tests' data constructors and seeds
(test-ordinal, test-ordinal-fitted, test-review-v29, test-brms-shapes,
test-brms-shapes-punch, test-brms-likelihood, test-brms-agreement's
`inhaler`, test-mv-gaps, test-thres), plus a random intercept and
weights on the test-ordinal data. The "before" arm ran on the
consolidated branch at a897706d installed into the lane library before
any edit, the "after" arm on the change.

On the 16 models whose thresholds do not cross:

- relative logLik difference at most 5.3e-11;
- the largest |estimate after - before| over all fixef() rows, divided
  by that row's standard error, at most 1.6e-4.

That is optimizer precision: the two arms stop at max |gradient|
between 4e-8 and 8e-4. Two notes: `thres_5` has an unidentified top threshold
(estimate near 19.9, standard error 5052) that moved 0.27, 5.4e-5 of
its standard error; `v29_43_cs` (`y ~ x + cs(x)`, not identified, see
"Defects seen") has no finite standard error after the change and a
standard error of 9.2e5 before, and its estimates moved at most 9.0e-6.

On the three that sat on the boundary, the fit moves, as it must:

| model | smallest gap before | logLik before | logLik after | gap after |
|---|---|---|---|---|
| inhaler, `rating ~ period + carry + cs(treat)` | 1.8e-8 | -455.114992284 | -451.256799338 | -1.87 |
| test-mv-gaps multivariate, response o2 | 6.3e-9 | -538.668768950 | -535.223385040 | -0.875 |
| thres(gr = g), seed 11 without the random effect | 5.6e-8 | -313.010364748 | -312.980708012 | -0.109 |

### Where they cross, frmtmb now reaches brms's mode

`dev/sratio-brms-stan.R`, log `dev/sratio-brms-stan-log.txt`. brms
2.23.0 writes the Stan program, rstan 2.32.7 (StanHeaders 2.32.10)
compiles it, and `log_prob(adjust_transform = FALSE)` and its gradient
are evaluated at frmtmb's estimate. brms's ungrouped `Intercept` is the
threshold minus `means_X' b` (it centers `X` and not `Xcs`); under
`thres(gr = )` it does not center.

Flat priors, against frmtmb's logLik at the maximum likelihood
estimate:

| case | brms log_prob | frmtmb logLik | difference | brms max abs gradient |
|---|---|---|---|---|
| grouped | -312.98070801181348 | -312.98070801181348 | 0 | 2.2e-4 |
| inhaler | -451.25679933758295 | -451.25679933758295 | 0 | 4.3e-5 |
| o2 | -126.85411925629103 | -126.85411925629103 | 0 | 7.5e-5 |

On the old build the grouped case read -313.0104 with a brms gradient
of +0.544 and -0.544 on the two thresholds held equal
(`dev/thres-sratio-order-log.txt`). The multivariate fit of
test-mv-gaps (cumulative o, sratio o2, gaussian y1) shares no
parameter, and its logLik equals the sum of the three univariate fits'
to 5.0e-12, relative; its o2 thresholds, 0.0904 and -0.7844, cross.

Classical reference: with `cs()` the stopping-ratio likelihood is
`K - 1` logistic regressions on shrinking subsets (Laara and Matthews
1985), whose intercepts are unconstrained. `test-sratio-thresholds.R`
asserts it on crossing data (relative logLik difference below
sqrt(eps), intercepts within 1e-3 standard errors), together with the
closed-form maximum of the intercept-only model, the link of each
observed hazard, for the logit, probit and cloglog links, and per
level under `thres(gr = )`. On the old build the glm check failed with a
relative logLik difference of 0.092, the closed form with an estimate
6.6 standard errors away.

### The MAP convention: no Jacobian, as in brms

With `normal(0, 2)` on class `"Intercept"`, against minus frmtmb's
penalized objective with NOTHING added (same script and log):

| family | case | brms log_prob | frmtmb | difference | brms max abs gradient |
|---|---|---|---|---|---|
| sratio | grouped | -328.24582770575796 | -328.24582770575796 | 0 | 1.1e-4 |
| sratio | inhaler | -456.81964153662528 | -456.81964153662528 | 0 | 9.1e-5 |
| cratio | grouped | -328.24582770575796 | -328.24582770575796 | 0 | 1.1e-4 |
| cratio | inhaler | -456.81964153662528 | -456.81964153662528 | 0 | 9.1e-5 |

brms's gradient vanishes at frmtmb's MAP, so frmtmb's MAP is brms's
mode, and neither carries a Jacobian: brms's unconstrained vector has
none, and frmtmb's entry is on scale `"internal"`. The inhaler rows also
exercise the centering offset. cratio is the control for the
convention, which it already followed. sratio and cratio agree to the
last digit because under a symmetric distribution function
`F(t - eta) = 1 - F(eta - t)`: the two are one model, an identity, not
a second check. (cumulative keeps its ordered map and its Jacobian:
`dev/thres-brms-stan-log.txt` has brms's gradient of -1 on each log
increment there.) `test-sratio-thresholds.R` asserts the MAP against a
one-dimensional `optimize()` per threshold and the objective at that
point against the log-likelihood plus the log prior; on the old build
the sratio objective there was off by 23.0.

### Fits saved by an earlier version

`dev/sratio-saved-fit.R`, log beside it. A fit carries its family's
closures, so an sratio fit made by the old build keeps reading its own
(first, log increments) storage. An inhaler fit saved by the base build
and read by this one prints the same `fixef()`, `fitted()`,
`variables()`, `hypothesis()` and `simulate()` (seeded) to every digit
printed. `update()` refits, and so gets the new storage. A `newparams =
list(tau_raw = )` or `start =` written for an old sratio fit is now
read as thresholds: that is the breaking part of the change.

## Speed

The change removes work from the objective: sratio's log-density no
longer builds the thresholds from log increments, a loop of `K - 1`
taped additions and exponentials per evaluation. Not timed, because
nothing was added.

## brms-suite ports

No port verdict changes. The ledger rows that mention sratio
(`priors:14`, `brmsfit-methods:837`, `families:35`, `standata:118`)
are about `threshold = "equidistant"`, residuals, a link name and a
response check, none about ordering. brms's `tests.distributions.R`
checks of `dsratio()` are not ported. The generator was not rerun, as
no verdict changed.

## Decisions and what was NOT done

- **cumulative stays ordered.** brms orders it too, and its category
  probabilities need it.
- **Start values unchanged.** `ord_tau_init()` still starts sratio at
  the link of the cumulative proportions; the old build started at the
  same thresholds (with increments floored at 0.05), so the optimizer
  starts where it did. The observed hazards would be the exact
  intercept-only answer; not needed for the maximum.
- **A family flag instead of the name test** in `priors.R` and in the
  test helper was considered and not done: the map is already on the
  family (`post$ord_thresholds`), and the prior code needs the scale
  of the map, which no family declares. One name, `"cumulative"`,
  replaces two.

## Defects seen, not fixed

- **`y ~ x + cs(x)` is not identified and is not refused.** The global
  slope and the `K - 1` category-specific slopes of the same variable
  span the same space (`test-review-v29.R`'s cs() block fits it). The
  standard errors were 9.2e5 on the old build and are NaN now (one
  negative eigenvalue of order 1e12 in the covariance); the estimates
  are the same to 9e-6. Not this lane's.
- **A prior draw on an unordered threshold vector would recycle one
  value.** `draw_prior_entry()` draws one number per entry, and a
  class `"Intercept"` entry on cratio, acat and now sratio covers the
  whole vector (`idx = seq_along(tau_raw)`), so `draw_prior_pars()`
  would write the same draw into every threshold. The comment in
  `test-simulate-ergonomics.R` says `frm_simulate()` stops earlier on
  an ordinal model, at the natural-scale `newparams` check; I did not
  construct a call that reaches the draw, so reachability is not
  settled here.

## Tests run

Every file one per process with `NOT_CRAN=true`, lane library unless
marked. RESULT lines are pasted from the runner.

    test-sratio-thresholds.R           pass=29 fail=0 err=0 skip=0
    test-sratio-thresholds.R (FRMTMB_LIB=base)
                                       pass=3 fail=22 err=2 skip=0
    test-ordinal.R                     pass=109 fail=0 err=0 skip=0
    test-ordinal-fitted.R              pass=129 fail=0 err=0 skip=0
    test-thres.R                       pass=67 fail=0 err=0 skip=0
    test-brms-shapes.R                 pass=72 fail=0 err=0 skip=0
    test-brms-shapes-punch.R           pass=55 fail=0 err=0 skip=0
    test-brms-shapes-punch2.R          pass=22 fail=0 err=0 skip=0
    test-review-v29.R                  pass=137 fail=0 err=0 skip=0
    test-simulate-density.R            pass=459 fail=0 err=0 skip=0
    test-diagnostics-ux.R              pass=127 fail=0 err=0 skip=0
    test-mv-gaps.R                     pass=46 fail=0 err=0 skip=0
    test-osa-inference.R               pass=34 fail=0 err=0 skip=0
    test-numerical-robustness.R        pass=709 fail=0 err=0 skip=0
    test-brms-families.R               pass=1033 fail=0 err=0 skip=0
    test-compat.R                      pass=639 fail=0 err=0 skip=0
    test-predfix.R                     pass=92 fail=0 err=0 skip=0
    test-simulate-newdata.R            pass=44 fail=0 err=0 skip=0
    test-v15.R                         pass=41 fail=0 err=0 skip=0
    test-v17.R                         pass=35 fail=0 err=0 skip=0
    test-prior-location-dpar.R         pass=55 fail=0 err=0 skip=0
    test-brms-formula-priors.R         pass=208 fail=0 err=0 skip=0
    test-case-studies.R                pass=30 fail=0 err=0 skip=0
    test-parity-integration.R          pass=8 fail=0 err=0 skip=0
    test-brms-agreement.R              pass=168 fail=0 err=0 skip=2
    test-simulate-ergonomics.R         pass=50 fail=0 err=0 skip=0
    test-bracket-access.R              pass=33 fail=0 err=0 skip=0
    test-custom-family.R               pass=126 fail=0 err=0 skip=0
    gated (FRMTMB_BRMS_FIT_TESTS=true, /opt/rlib/stan first):
    test-brms-likelihood.R, block "row 12: ordinal families, cumulative
      sratio cratio acat" alone (test_file(desc = ), Stan compiled)
                                       pass=20 fail=0 err=0 skip=0
    test-brms-suite-families.R         pass=84 fail=0 err=0 skip=0
    test-brms-suite-methods.R          pass=162 fail=0 err=0 skip=0
    test-brms-suite-priors.R           pass=36 fail=0 err=0 skip=0
    test-brms-suite-standata.R         pass=87 fail=0 err=0 skip=0
    frmtmb.sample (lane frmtmb and frmtmb.sample):
    test-sratio-draws.R                pass=2 fail=0 err=0 skip=0
    test-sratio-draws.R (FRMTMB_LIB=base)
                                       pass=0 fail=2 err=0 skip=0
    test-thres-draws.R                 pass=4 fail=0 err=0 skip=0
    test-brms-shapes-draws.R           pass=61 fail=0 err=0 skip=0
    test-draws-methods.R               pass=148 fail=0 err=0 skip=0
    test-predfix-new-levels.R          pass=62 fail=0 err=0 skip=0

The ungated brms-suite files report 0 of everything without the gate,
which is why they were run gated. The rest of test-brms-likelihood.R
was not run: its other blocks compile a Stan model each and do not
reach sratio. The vignettes were not built.

On the base build (`FRMTMB_LIB=base`) `test-sratio-thresholds.R` fails
in every block: five fail behaviorally (the second threshold is held
equal to the first, the closed form is 6.6 standard errors away, sratio
and cratio disagree by 0.082 relative in logLik, the glms by 0.092, the
MAP is 6.4 standard errors away and the objective 23.0 off), and two
error because base predates `thres()` and ordinal responses in
multivariate models. The messages are in
`dev/sratio-seen-failing-base.txt`. The frmtmb.sample test fails both
expectations on base.

## NEWS entry

- BREAKING: `sratio()`'s thresholds are unconstrained, as brms declares
  them, and no longer kept increasing. frmtmb held them as (first
  threshold, log increments), like `cumulative()`, so where the maximum
  has crossing thresholds the fit stopped on the ordering boundary,
  with two thresholds equal, and was not brms's mode. brms's own
  `inhaler` example with `cs(treat)` was one such fit: its log-likelihood
  rises from -455.115 to -451.257. brms's compiled log density now
  equals frmtmb's log-likelihood there with a vanishing gradient. Where
  the thresholds do not cross, a fit is unchanged to optimizer
  precision. The internal `tau_raw` of an sratio fit is now the
  thresholds themselves, so a `newparams` or `start` written for the
  old storage means something else. A class `"Intercept"` prior on
  sratio's thresholds carries no Jacobian any more, as in brms and as
  on `cratio()`.
