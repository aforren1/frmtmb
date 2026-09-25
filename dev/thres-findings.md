# Lane thres: brms's `thres()` addition term

Round of 2026-09-25. Worktree `/home/user/wt/thres`, branch
`lane/thres`, private library `/opt/rlib/lane-thres`.

## What changed

`y | thres(x, gr) ~ ...` is accepted on `cumulative()`, `sratio()`,
`cratio()` and `acat()` and means what it means in brms 2.23.0:

- `thres(K)` or `thres(x = K)` fits `K` thresholds, so the response
  has `K + 1` categories even when the top ones are never observed.
- `thres(gr = g)` gives each level of `g` a threshold vector of its
  own. The vectors are merged into one parameter vector in level order,
  as brms merges them (`merged_Intercept`, `Jthres`), and the linear
  predictor is shared.
- `thres(n, gr = g)` or `thres(x = n, gr = g)` sets a count per level;
  `n` must be constant within a level.

Files:

| file | change |
|---|---|
| `R/thres.R` (new) | the group codes, the count per level, the family finalizer, the grouped densities of the four families, the grouped simulator and category probabilities, newdata codes, brms's `[g,k]` labels, the fit-end warning |
| `R/parse.R` | `thres` joins `core_aterms`; `thres(x, gr)` is matched like `brms:::resp_thres(x, gr)`; `thres_gr` maps to `thres()` in `aterm_base()` and `aterm_spelling()`; a second `thres()` is a duplicate |
| `R/families.R` | the four ordinal families declare `accepts_aterms = c("weights", "thres")` and `family_finalize = thres_finalizer(...)`; `ord_tau_init()` takes `K` |
| `R/frame.R` | `thres(gr = )` is coded by `factor()` levels; `thres()` on a non-ordinal family is refused by name; `cs()` with grouped thresholds is refused (brms's sentence); the `cs()` coefficient count comes from the threshold count, not `max(y)`; `thres(x = )` beyond an ordered factor's levels is refused |
| `R/predict.R` | `ordinal_ncat()` is the largest group's; the category probabilities (`fitted()`, `frm_linpred(type = "response")`, the delta-method SEs, `conditional_effects()`) receive each row's group; `aterms_for_newdata()` codes `thres(gr = )` against the fitted levels and skips the count |
| `R/predict-brms.R` | `predict(newdata = )` hands the simulator each row's group |
| `R/priors.R` | class `Intercept` on grouped thresholds is one entry per level, each its own ordered vector, with no centering offset; `group = "<level>"` selects one level; `group =` on an ungrouped model is refused |
| `R/brms-shapes.R`, `R/confint.R` | threshold names carry the group: `Intercept[a,1]`, `b_Intercept[a,1]` |
| `R/compat.R` | `thres()` joins the aterm vocabulary with 22 rules |
| `R/fit.R`, `man/frm.Rd` | new section "Ordinal thresholds, thres()"; the G2.4 srrstats note names the `thres(gr = )` coding |
| `vignettes/brms-migration.Rmd`, `SPEC.md` | one bullet; `thres` joins the aterm list |
| `extensions/frmtmb.sample/R/methods-draws.R` | `posterior_predict(newdata = )` hands the simulator each row's group |
| `tests/testthat/test-thres.R` (new), `extensions/frmtmb.sample/tests/testthat/test-thres-draws.R` (new) | see Tests |
| `tests/testthat/helper-brms-suite.R` and its sample copy | the `standata()` shim reports `nthres`, `ngrthres`, `Jthres` |
| `dev/brmsport-verdicts.tsv`, `dev/brmsport-verdicts-manual.tsv`, `tests/testthat/test-brms-suite-standata.R` | 10 port rows flipped; see "brms-suite ports" |
| `dev/thres-validate.R`, `dev/thres-timing.R`, `dev/thres-brms-stancode.R` and their `-log.txt` | evidence |

### How it is built

A count alone changes only the length of the threshold vector, so it
reaches the existing densities through `extra_pars`. Groups change the
density. The response is needed to resolve them, so the ordinal
families use the existing `family_finalize()` slot: at frame assembly it
counts the thresholds per level, then REPLACES the family's `lpdf`,
`sim`, `extra_pars` and `post$ord_thresholds` with grouped versions,
and stores the layout on the family as `fam[["thres"]]`. Every later
stage reads the family, so prediction and simulation on new data need
nothing but the row's group code, which rides in `aterms$thres_gr` like
any per-row addition term (so `imp_expand()` and the other row
subsetters handle it unchanged). A model without `thres()` gets its
family back untouched.

The grouped densities are vectorized over rows; their loops run over
threshold positions (at most `max(nthres)`) or over the parameter
layout, as the ungrouped ones do. Each row reads its own slice by data
indices built once per tape. The cumulative interior term is evaluated
only on the rows with an interior category, so no row evaluates a
threshold pair outside its own slice (the ungrouped code clamps
instead, which a slice of length 1 cannot). A category past a row's own
count gets log probability `-Inf`, which makes the density a proper pmf
on `1..max(nthres) + 1` for every row and lets `fitted()` read category
probabilities out of it unchanged.

## Validation

### brms's density at a shared parameter point

`dev/thres-validate.R`, log `dev/thres-validate-log.txt`. Data seed 11,
n = 240, three levels with 4, 3 and 5 categories. frmtmb's objective at
a chosen point (slope 0.55, thresholds per level fixed by hand) against
the sum of brms's own R densities (`brms:::dcumulative`, `dsratio`,
`dcratio`, `dacat`, the functions `posterior_epred()` uses), each row
read through its own level's thresholds as the generated Stan code does
with `Jthres` (`dev/thres-brms-stancode-log.txt`):

| family | link | frmtmb log-lik | brms | relative difference |
|---|---|---|---|---|
| cumulative | logit | -320.971733615564 | -320.971733615564 | 5.3e-16 |
| cumulative | probit | -352.235618992248 | -352.235618992248 | 3.2e-16 |
| cumulative | cauchit | -320.226039178067 | -320.226039178067 | 1.8e-16 |
| sratio | logit | -359.947552754553 | -359.947552754553 | 3.2e-16 |
| sratio | probit | -405.636057698747 | -405.636057698747 | 2.8e-16 |
| cratio | logit | -359.947552754553 | -359.947552754553 | 3.2e-16 |
| cratio | probit | -405.636057698747 | -405.636057698747 | 2.8e-16 |
| acat | logit | -357.650855471611 | -357.650855471611 | 1.6e-16 |

cauchit takes the plain-difference branch of the cumulative density,
the others the log-space branch. sratio and cratio agree with each
other here because the logistic and the normal are symmetric; that is
an identity of the model, not a second check.

`thres(5)` without groups (seed 14, top two categories unobserved), the
same comparison: relative differences 1.6e-16, 0, 0 and 0 for the four
families.

`tests/testthat/test-thres.R` repeats the grouped comparison against a
second, independent implementation written from brms's Stan functions
inside the test, with the tolerance `64 * eps * sum(|terms|)`.

### Classical references at the ML optimum

- MASS::polr, seed 12, n = 600: with a slope per level
  (`y | thres(gr = g) ~ g:x`) the model factorizes into one polr fit per
  level. frmtmb log-lik -763.3918321075, sum of the three polr fits
  -763.3918321069, relative difference 7.6e-13; the nine thresholds
  agree to 3.2e-06 relative to the largest.
- ordinal::clm, seed 13, n = 600, two levels with 4 categories each:
  `clm(factor(y) ~ x, nominal = ~ g)` is a shared slope with a threshold
  vector per level. log-lik -740.0992040582 against -740.0992040564,
  relative 2.4e-12; slope 0.81822529 against 0.81822544; thresholds to
  3.3e-06.

The tests assert these as `|difference| / |log-lik| < sqrt(eps)` and
`|difference| / standard error < 1e-3`, both ratios to quantities the
run measures.

### Priors

A class `Intercept` prior on grouped thresholds is one density per
level, each on an ordered vector of its own, with its own log-Jacobian
and no centering offset (brms does not center the design under grouped
thresholds: `mu += X * b`, `b_Intercept_g = Intercept_g`). The test
checks the objective difference with and without `normal(0, 2)`
against `sum_g [sum log dnorm(tau_g) + sum log diff(tau_g)]` to
`64 * eps * |objective|`, and `group = "b"` against level b's term
alone.

## Downstream methods

| method | grouped thresholds | how checked |
|---|---|---|
| `summary()`, `fixef()`, `vcov()`, `confint()` | rows `Intercept[a,1]`... in merged order, then the slopes | test-thres.R; `confint()` keeps its internal `tau_raw_k` rows, as without groups |
| `variables()`, `hypothesis()` | `b_Intercept[a,1]`, brms's `rename_thres()` names | test-thres.R |
| `fitted()`, `predict()`, `frm_linpred(type = "response")` | `max(nthres) + 1` columns; a column past a row's own categories is 0, as in brms's `posterior_epred()`; rows sum to 1; newdata needs the grouping variable and known levels, refused by name otherwise | test-thres.R |
| `simulate()`, `predict()` draws, `frm_simulate()` | each row from its own level's thresholds, in sample and on newdata | test-thres.R; `frm_simulate(newparams = )` by hand |
| `conditional_effects()` | runs; the grid carries the grouping variable at its reference level, and `effects = "x:g"` shows each level | by hand (`/tmp/lanes/thres/thres-explore3.R`) |
| `emmeans` | the latent predictor, which no threshold enters | test-thres.R |
| random effects (Laplace, `quadrature = TRUE`, `importance = 400`, `REML = TRUE`), `coef()` | fit; `coef()` shifts every `Intercept[g,k]` by the group mode | test-thres.R (Laplace), by hand (the others) |
| `mo()` | fits; it changes only the latent predictor | test-thres.R |
| `weights()` | weights of 2 give the fit of the duplicated data (relative log-lik difference below sqrt(eps)) | test-thres.R |
| `cs()` | refused with `gr`, in brms's sentence; with `thres(x = K)` it takes K coefficients | test-thres.R |
| `residuals(type = "osa")` | refused with `gr`, naming dharma_residuals(); works with `thres(x = K)` | test-thres.R |
| `anova()`, `update()`, `confint(method = "profile")`, `dharma_residuals()`, `pp_check()` | run | by hand |
| `frm_sample()` draws: `fixef()`, `posterior_epred()`, `posterior_predict(newdata = )` | group names; zero columns; each row within its level | frmtmb.sample test-thres-draws.R |
| `prior_summary()` | lists the group-specific row | by hand |

## Decisions and what was NOT done

- **Counts follow brms exactly.** brms counts per level from the highest
  category observed there, and an ordered-factor response takes its
  categories from the levels that occur (unused levels are dropped with
  the model frame). Measured on brms 2.23.0: a factor with levels 1..4
  and 1..3 observed gives `nthres` 2, and 2 and 1 per level
  (`dev/thres-brms-stancode-log.txt`). frmtmb's counts are the same.
- **Unobserved categories under ML.** A threshold above a level's
  highest observed category is not identified: its estimate runs off
  toward infinity (`thres(5)` on data reaching 4 gave 19.5 with a
  standard error of 7184). brms never meets this because its default
  `student_t(3, 0, 2.5)` prior holds every threshold. frmtmb now warns
  at the end of the fit, naming the thresholds, unless a class
  `Intercept` prior covers them. Checked both ways in test-thres.R.
- **`thres(x = )` beyond an ordered factor's levels is refused.** brms
  hands back bare codes there; frmtmb returns simulated responses as the
  response's own factor, which cannot hold a category with no level.
  The refusal says to code the response as integers.
- **`cs()` with grouped thresholds is refused**, as brms refuses it.
- **`residuals(type = "osa")` with grouped thresholds is refused.** The
  one-step density picks the category arithmetically over one shared set
  of categories (`ord_cat_sel()`); a per-row set would need a second
  version of every OSA branch. `dharma_residuals()` covers the need.
- **`threshold = "equidistant"` and `"sum_to_zero"`** are not
  arguments of frmtmb's ordinal families at all (unchanged, and outside
  this lane: `priors:14` in the port verdicts already records it). So
  their combination with `thres()` does not arise.
- **Refits whose data lose a category.** `frm_bootstrap()` and other
  simulate-and-refit paths recount the thresholds on each simulated
  data set, as they do without `thres()`. A replicate in which a level
  never reaches its top category refits with one threshold fewer.
  Pinning the count with `thres(x = , gr = )` avoids it. Not changed.

## Defect found and fixed

`set_prior(..., class = "Intercept", group = "a")` on an ordinal model
WITHOUT grouped thresholds was applied to the whole threshold vector,
silently: `ordinal_threshold_entry()` never read `group`. It is now
refused, naming the reason. The test
`group = on an ungrouped threshold prior is refused` fails on the base
build with "Expected `frm(...)` to throw a error" (seen, RESULT below).

## Defects seen, not fixed

- `tests/testthat/test-pp-check-types.R` fails 4 expectations in block
  "every bayesplot ppc type does on a fit what it does in brms" on the
  base build too (`FRMTMB_LIB=base`: pass=168 fail=4), so it is not this
  lane's.
- At n = 6000 (`dev/thres-timing.R`), the ungrouped cumulative fit
  already warns "Large maximum absolute gradient at the optimum
  (0.001)" on the base build, and so do the grouped cumulative, cratio
  and acat fits (0.0012, 0.0059, 0.0040) with the objective near 8,700.
  The warning's threshold looks absolute; not investigated here.

## Speed

`dev/thres-timing.R`, log `dev/thres-timing-log.txt`, seed 21, n =
6000, three levels with 4, 3 and 5 categories.

- A model without `thres()` is untouched: on each of the four families
  the base and lane builds print the same objective, iteration count
  and estimates to 17 significant digits.
- The grouped model against the ungrouped one, interleaved in one
  process, blocks past 1.2 s, minimum of 5 rounds, with a control built
  from the same call twice:

| family | control B / A | grouped / A |
|---|---|---|
| cumulative | 1.15 | 0.92 |
| sratio | 0.86 | 1.13 |
| cratio | 1.01 | 0.98 |
| acat | 1.12 | 0.96 |

  Every grouped ratio lies inside the spread the control shows on
  identical work (0.86 to 1.15), so the instrument does not separate
  the two; the per-evaluation work is the same number of passes over
  the rows.

## brms-suite ports

Ten `standata` rows flip from "cannot transfer" to pass:
`standata:1063`, `:1065`, `:1072`, `:1073`, `:1074`, `:1075`, `:1078`,
`:1079`, `:1082`, `:1083` (block "standata handles grouped ordinal
thresholds correctly"). The shim in `helper-brms-suite.R` now reports
brms's `nthres` (per level), `ngrthres` and `Jthres`; the sample
extension's copy was synced, and its copy test passes.

The generator could NOT run: `dev/brmsport-blocks.R` needs
`dev/brms-suite/brms_2.23.0.tar.gz` (sha256-checked), which is
gitignored and absent from this machine, and CRAN is unreachable. So
`tests/testthat/test-brms-suite-standata.R` was edited by hand to what
the generator emits for a pass verdict with an empty reason (compare
`standata:1060`), the ten rows were set to `pass` with an empty reason
in `dev/brmsport-verdicts.tsv` and removed from
`dev/brmsport-verdicts-manual.tsv` (the ledger refuses a manual verdict
on an assertion that holds). `dev/brmsport-ledger.tsv` and
`dev/brmsport-log/` were not regenerated; they need the recorder and the
tarball. `standata:171` (the deprecated `cat()` spelling) stays
"cannot transfer".

## Tests run

Every file one per process with `NOT_CRAN=true`, lane library unless
marked:

    test-thres.R                      pass=67 fail=0 err=0 skip=0
    test-thres.R (FRMTMB_LIB=base)    pass=0 fail=1 err=10 skip=0
    test-brms-suite-standata.R (FRMTMB_BRMS_FIT_TESTS=true)
                                      pass=87 fail=0 err=0 skip=0
    test-brms-suite-standata.R (base, gated)
                                      pass=77 fail=10 err=0 skip=0
    test-compat.R                     pass=607 fail=0 err=0 skip=0
    test-ordinal.R                    pass=109 fail=0 err=0 skip=0
    test-ordinal-fitted.R             pass=129 fail=0 err=0 skip=0
    test-custom-family.R              pass=126 fail=0 err=0 skip=0
    test-bracket-access.R             pass=33 fail=0 err=0 skip=0
    test-cens-trunc.R                 pass=69 fail=0 err=0 skip=0
    test-brms-priors.R                pass=0 fail=0 err=0 skip=12 (gated)
    test-brms-likelihood.R            pass=20 fail=0 err=0 skip=33
    test-mo-terms.R                   pass=74 fail=0 err=0 skip=0
    test-osa-inference.R              pass=34 fail=0 err=0 skip=0
    test-simulate-density.R           pass=435 fail=0 err=0 skip=0
    test-ce-bands.R                   pass=167 fail=0 err=0 skip=0
    test-ce-facets.R                  pass=92 fail=0 err=0 skip=0
    test-brms-shapes-punch.R          pass=55 fail=0 err=0 skip=0
    test-brms-shapes-punch2.R         pass=22 fail=0 err=0 skip=0
    test-numerical-robustness.R       pass=679 fail=0 err=0 skip=0
    test-structure.R                  pass=135 fail=0 err=0 skip=0
    test-importance.R                 pass=226 fail=0 err=0 skip=0
    test-method-residue.R             pass=64 fail=0 err=0 skip=0
    test-pp-check-types.R             pass=168 fail=4 (same on base)
    test-brms-families.R              pass=986 fail=0 err=0 skip=0
    test-v17.R                        pass=35 fail=0 err=0 skip=0
    test-review-v29.R                 pass=137 fail=0 err=0 skip=0
    frmtmb.sample test-thres-draws.R  pass=4 fail=0 err=0 skip=0
    frmtmb.sample test-thres-draws.R (lane core, base sample)
                                      pass=0 fail=0 err=1 skip=0
    frmtmb.sample test-brms-suite-helper-copy.R (gated)
                                      pass=33 fail=0 err=0 skip=0
    frmtmb.sample test-draws-methods.R
                                      pass=148 fail=0 err=0 skip=0
    frmtmb.sample test-predfix-new-levels.R
                                      pass=62 fail=0 err=0 skip=0

On the base build test-thres.R fails in every block: ten error because
base refuses `thres()` at parse time, and the eleventh, the prior
`group =` block, fails behaviorally (the prior was accepted). With the
lane's core and the base frmtmb.sample, the sample test errors with the
lane's "rows reached the density without their group" refusal, the
pin for the `methods-draws.R` change; without that change the failure
is loud, not a wrong draw.

## NEWS entry

- The ordinal families take brms's `thres()` addition term.
  `y | thres(K) ~ x` fits `K` thresholds, and `y | thres(gr = g) ~ x`
  fits one threshold vector per level of `g`, with a count per level
  from the data or from `thres(n, gr = g)`. The log-likelihood agrees
  with brms's densities to about 5e-16, relative, on all four families,
  and with MASS::polr and ordinal::clm at the optimum. The thresholds
  have brms's names, `Intercept[g,k]` and `b_Intercept[g,k]`;
  `set_prior(class = "Intercept", group = "g")` reaches one level.
  `fitted()` and `predict()` give probability 0 to a category past a
  row's own, as brms does. `cs()` and `residuals(type = "osa")` are
  refused with `gr`. A threshold above the highest observed category
  is not identified without a prior, and the fit warns about it.
- A class `"Intercept"` prior with `group =` on an ordinal model
  without grouped thresholds is refused. It used to be applied to every
  threshold without a word.
