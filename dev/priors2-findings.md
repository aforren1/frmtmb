# PRIORS-CONSISTENCY lane (wt-priors2)

The decision this lane implements, verbatim: "re: sigma, I'd like to
maintain consistency with brms". That is option 3 of
`dev/priors-findings.md` ("flip the meaning of the existing
spelling"), which the review at `dev/reviews/2026-09-05-priors.md`
recommended against on cost grounds. The maintainer overrode it.

## brms's own acceptance rules, re-measured here

`brms::make_stancode()`, brms 2.23.0, no compilation. Each cell is the
`lprior` line brms writes, or its refusal.

| model | `class = <dpar>` | `class = "Intercept", dpar = <dpar>` |
| --- | --- | --- |
| `y ~ x`, gaussian | `normal_lpdf(sigma \| 0, 3)` | ERROR "do not correspond to any model parameter: Intercept_sigma" |
| `bf(y ~ x, sigma ~ 1)` | ERROR "... : sigma" | `normal_lpdf(Intercept_sigma \| 0, 3)` |
| `bf(y ~ x, sigma ~ x)` | ERROR "... : sigma" | `normal_lpdf(Intercept_sigma \| 0, 3)` |
| `cnt ~ x`, negbinomial | `normal_lpdf(shape \| 0, 3)` | ERROR "... : Intercept_shape" |
| `cnt ~ x`, zero_inflated_poisson | `normal_lpdf(zi \| 0, 3)` | (not measured) |
| `bf(cnt ~ x, zi ~ 1)`, zip | ERROR "... : zi" | `normal_lpdf(Intercept_zi \| 0, 3)` |
| `p ~ x`, Beta | `normal_lpdf(phi \| 0, 3)` | (not measured) |
| `bf(p ~ x, phi ~ 1)`, Beta | (not measured) | `normal_lpdf(Intercept_phi \| 0, 3)` |

Two things this pins that the review left implicit:

1. the rule generalizes past `sigma`: `shape`, `zi` and `phi` behave
   the same way, so it is a rule about dpars and not about one name;
2. **`sigma ~ 1` counts as a predictor.** brms refuses `class =
   "sigma"` on `bf(y ~ x, sigma ~ 1)` exactly as it does on
   `sigma ~ x`. The distinction is "did the user write a formula",
   not "does the design have more than an intercept".

## Can frmtmb make that distinction?

Yes, from the spec alone (`resolve_priorlist()` gets `spec` and
`frame`, never `bform`: `R/fit.R:722` and `R/par-template.R:157` call
it with a two-element shim). Measured field sets of
`fit$spec$responses[[1]]$dpars$sigma`:

| model | fields |
| --- | --- |
| `bf(y ~ x)` | `name link fixed re rhs smooth constant` |
| `bf(y ~ x, sigma ~ 1)` | `name link constant fixed re smooth mo miterms csterms gpterms carterms spdeterms acterms rhs` |
| `bf(y ~ x, sigma ~ x)` | same 14 as above |

`plain_dpar()` (`R/parse.R:1244`) builds the first; a written formula
goes through `parse_linpred()`, which adds the special-term slots. The
design matrix is `(Intercept)` only in the first two, so `ncol(X)` does
NOT separate them. The written-formula marker does.

## What landed

| # | change | file:line |
| --- | --- | --- |
| 1 | `set_prior()` takes a distributional parameter's own name as a class and marks it `natural` | `R/priors.R:404-405`, `:452-470` |
| 2 | the class vocabulary is open: frmtmb's ten own names, the brms classes refused BY NAME (now from the native path too), then any dpar name | `R/priors.R:479` (`frmtmb_prior_classes`), `:493` (`check_dpar_prior_class`) |
| 3 | the brms route is a pass-through: one rule, written once | `R/priors.R:838` (`brms_prior_route`), `:706-709` |
| 4 | the shape gate, brms's own acceptance rule, by name in both directions | `R/priors.R:1921` (`dpar_shape_refusal`), called at `:1646` |
| 5 | `dpar_has_predictor()` reads the spec's record of who wrote the formula | `R/priors.R:1896` |
| 6 | `set_prior()` refuses an unhonored `coef` on `sd`/`cor` (review residual R1) | `R/priors.R:444-445` |
| 7 | a natural spec prints and reports as the class it was WRITTEN with | `R/priors.R:1338` (`spec_spelling`), print at `:963`, `spec_target()` at `:1313` |
| 8 | `get_prior()` lists a predictor-free dpar under its own class, as brms does | `R/priors.R:1144-1148` |
| 9 | the sample route keys its defaults by the spelling the table lists | `R/sampling-api.R:335-338`, exported at `:231` |

No `scale =` or `natural =` argument was added to `set_prior()`: after
the flip the placement follows from the class, which is brms's own
rule, so there is nothing left for such an argument to override. The
`natural` field survives as an internal marker only, because
`frmtmb.sample:::natural_dpar_prior()` writes it directly.

## Measured, before and after

Same data, same seeds, same starts; `p2-before2.R` / `p2-after.R`.

| shape | spelling | before | after |
| --- | --- | --- | --- |
| `y ~ x`, gaussian, n = 200 (ML sigma 2.010268914) | `set_prior("normal(0, 0.5)", class = "sigma")` | error (class refused) | **sigma 1.942946001** |
| the same | `set_prior(..., class = "Intercept", dpar = "sigma")` | sigma 1.996511359 | **refused by name** |
| the same | `brms::prior(normal(0, 0.5), class = "sigma")` | sigma 1.942946001 | 1.942946001 (unmoved) |
| `bf(y ~ x, sigma ~ 1)` | `class = "Intercept", dpar = "sigma"` | sigma 1.996511359 | 1.996511359 (unmoved) |
| `bf(y ~ x, sigma ~ 1)` | `class = "sigma"` | (unreachable) | **refused by name** |

The two routes now agree to every digit on the same model, which is
the point: `set_prior(class = "sigma")` and
`brms::prior(..., class = "sigma")` are one code path.

`get_prior()` tables, compared against `brms::get_prior()` row for row:

| model | brms | frmtmb before | frmtmb after |
| --- | --- | --- | --- |
| `y ~ x` | b, b/x, Intercept, **sigma** | b, b/x, Intercept, **Intercept+dpar=sigma** | b, b/x, Intercept, **sigma** |
| `bf(y ~ x, sigma ~ x)` | ... Intercept+dpar, b+dpar, b/x+dpar | same | same (unmoved) |

## The eam sign flip the earlier review priced does NOT happen

The review's case against the flip rested on two costs, and one of
them is not a cost of the flip as the maintainer specified it.

`extensions/frmtmb.eam/tests/testthat/test-surface.R:133` writes
`set_prior("normal(0, 0.1)", class = "b", dpar = "mu")`. That is class
`"b"`, which is link-scale on the SLOPES in brms and here, so this
change does not touch it. Measured on the same fixture (`ddm_fit()`,
seed 404, n = 400):

| quantity | before | after |
| --- | --- | --- |
| `mu.cond` unpenalized | 0.9165381216 | 0.9165381216 |
| `mu.cond` with the class "b" prior | **+0.2654012893** | **+0.2654012893** |
| `shrink` ratio | 0.2895692858 | 0.2895692858 |

The review's `-0.3301564534` is what a flip that ALSO retargeted class
`"b"` on a dpar would produce. brms does not do that, so neither does
this. The `(0, 1)` ratio assertion is left exactly as it was, and it
still passes; the block gains a positive test of the new spelling
(`class = "bs"`, `normal(1, 0.05)`, pulling bs from 1.382781771 to
1.27555188) and of the refusal of the link spelling on the same model.

## The other test that had to move, and it is not the shrink one

`extensions/frmtmb.eam/tests/testthat/test-surface.R:126` asserted
`all(c("bs", "ndt") %in% gp$dpar)`. Those two dpars have no predictor
in that model, so `get_prior()` now lists them under their own CLASS
with an empty `dpar`, which is brms's own row. The assertion moves to
`gp$class`. This is in the eam lane's file and outside the one block
the brief named, because the brief could not have known which block
`get_prior()`'s row change would reach.

`extensions/frmtmb.sample/tests/testthat/test-prior-route.R:26` moved
for the same reason and in the same way: it read the sample route's
sigma default off `gp$class == "Intercept" & gp$dpar == "sigma"`, and
the row is now `gp$class == "sigma"`. The DEFAULT still attaches
(`student_t(3, 0, 2.5)`), because `registered_prior_defaults()` keys by
the spelling the table lists.

## The substantive-differences table, refreshed

Only the rows this lane moves are restated; every other row of
`dev/priors-findings.md`'s table stands as written there.

| brms row | frmtmb after this lane | substantive difference remaining | measured consequence |
| --- | --- | --- | --- |
| `class = "sigma"` and every other dpar class, from `brms::prior()` | a density on the parameter itself, through its inverse link | none | unmoved: sigma 1.942946001 on the gaussian shape before and after |
| the same class from frmtmb's OWN `set_prior()` | **the same density, the same code path** | **none. The lane's row 7 divergence is gone** | the two routes give one objective at 1e-12 on six shapes covering Intercept, sd, cor, sigma, shape, zi and phi (`test-brms-priors.R`) |
| `class = "Intercept", dpar = <name>` | the link-scale intercept of that dpar's predictor | none | brms writes `Intercept_sigma` for exactly the models frmtmb accepts it on |
| either spelling on the WRONG model shape | refused by name, naming the other | none | brms errors "do not correspond to any model parameter"; frmtmb errors naming the spelling that applies |
| `class = "b", dpar = <name>` | the dpar's slopes, link scale | none | unmoved: the eam `mu.cond` stays +0.2654012893 |
| `get_prior()`'s dpar row | its own class where the dpar has no predictor | none | matches `brms::get_prior()` row for row on `y ~ x`, `sigma ~ 1` and `sigma ~ x` |

### Every spelling brms honors that frmtmb still cannot

Swept 17 families through `brms::get_prior(y ~ x, family = )` and
compared the non-`b`, non-`Intercept` classes brms writes against the
dpars frmtmb's own family gives the same model:

`sigma`, `nu`, `alpha`, `beta`, `shape`, `phi`, `zi`, `hu`, `kappa`,
`ndt`, `quantile` - **every one is a dpar frmtmb has**, so every one is
now a class `set_prior()` takes. gaussian, student, lognormal,
skew_normal, exgaussian, Gamma, weibull, Beta, negbinomial, poisson,
bernoulli, zero_inflated_poisson, hurdle_gamma, cumulative, von_mises,
shifted_lognormal and asym_laplace all resolve with no class left over.
The families that carry a class frmtmb cannot name are the ones whose
CLASS was already refused for a structural reason (`simo`, `sds`,
`sdgp`, `lscale`, `sdcar`, `car`, mixture `theta`), and that list is
unchanged by this lane.

So the answer to the brief's question is: **none**. There is no dpar
spelling left whose brms meaning frmtmb cannot honor. What remains
refused is refused because frmtmb holds the QUANTITY somewhere else,
not because it spells the prior differently.

## What this cost, and what it did not

| item | cost |
| --- | --- |
| `print(prior_summary())` on a natural spec | the line changed: `class=Intercept dpar=sigma scale=natural` becomes `class=sigma scale=natural`. Deliberate: the old line printed a spelling the same model now refuses, and `print.frmtmb_priorlist()` promises what it prints can be pasted back into `set_prior()` |
| `frm_sample()`'s defaults | unchanged in substance and in sampling. `default_priors_for()` still builds `class = "Intercept", dpar = "sigma"` + `natural` for a sigma with no predictor and the plain link-scale spec for one with a predictor, which is brms's own distinction, and `registered_prior_defaults()` keys it onto the table's `sigma` row |
| `frm_sample()`'s disclosure MESSAGE | still labels that default `Intercept (sigma)  [natural scale]`. `announce_default_priors()` (`extensions/frmtmb.sample/R/sample.R:728`) reads `s$class` directly. Now inconsistent with the print and with `get_prior()`. Left to the wt-sample-ce lane, whose file it is; `spec_spelling()` is exported for it (`R/sampling-api.R:224`) and the fix is one line |

## Files touched

Core, mine by the brief:

- `R/priors.R` - the class vocabulary, the normalization, the shape
  gate, `get_prior()`'s dpar row, the printed spelling, the roxygen.
- `R/sampling-api.R` - the default-key spelling, plus `spec_spelling`
  added to the exported seam so the sampling package can label a
  natural default the way the table does.
- `tests/testthat/helper-brms-priors.R`, `test-brms-priors.R`,
  `test-prior-compat.R`, `test-get-prior-route.R`,
  `test-simulate-ergonomics.R`.
- `NEWS.md`, `vignettes/brms-migration.Rmd`, `SPEC.md`,
  `vignettes/inputs.Rmd`, `man/` via roxygen, `NAMESPACE`.

Outside the brief's list, each forced by the change and each minimal:

- `R/simulate-new.R` - two lines of ONE roxygen example, which wrote
  `class = "Intercept", dpar = "sigma"` on a model with no sigma
  formula and would now stop under `R CMD check`.
- `vignettes/frmtmb.Rmd` - one prior-predictive chunk, the same
  spelling, same reason. It now writes `exponential(1)` on
  `class = "sigma"`, a positive density on a positive parameter.
- `extensions/frmtmb.eam/tests/testthat/test-surface.R` - the
  `get_prior()` column assertion (see above), plus two assertions on
  the new spelling in the same block. The shrink assertion the brief
  named is untouched, because it did not move.
- `extensions/frmtmb.sample/tests/testthat/test-prior-route.R` - one
  line, the same `get_prior()` column.

## Verification

One `test_file()` per process, private library
`.../scratchpad/p2-lib` built from this worktree, my own copy of
`dev/stan-cache` (54 programs), `R_MAKEVARS_USER` at its
`makevars-cxx17.mk`. Every Stan program was a cache hit; nothing
compiled.

| run | result |
| --- | --- |
| `test-brms-priors.R` (gated) | **103 pass**, 0 fail (was 83 before the lane; +20 from the route-agreement block and the retargeted row 5) |
| `test-brms-likelihood.R` (gated) | **373 pass**, 0 fail, unmoved |
| `test-prior-compat.R` | **185 pass**, 0 fail (was 161) |
| `test-setprior.R` | 27, unmoved |
| `test-get-prior-route.R` | **37** (was 36; one row assertion became two) |
| `test-priors-autocor-classes.R` | 63, unmoved |
| `test-priors-bounds-grcov.R` | 49, unmoved |
| `test-simulate-ergonomics.R` | **50 pass** (was 45; two new assertions and the block that used to stop) |
| `test-sample-direct.R` (frmtmb.sample) | 136, unmoved |
| `test-message-uniqueness.R` | 6, unmoved |
| `test-bracket-access.R` | 8, unmoved |
| **full core suite**, one file per process | **117 of 117 files ran**, list `diff`ed against `ls tests/testthat/test-*.R` (identical, 0 duplicates), **6562 pass, 0 fail**, 93 skip, 2 errors |

The 2 errors are `test-influence-plot.R:68` and `:130`, both
`local_mocked_bindings()` without a `.package` argument, which needs
pkgload and so fails under a bare `test_file()` harness and not under
`test_check()`. `R CMD check`'s own `checking tests` runs the same file
green, which is the confirmation. The same artifact hits
`frmtmb.sample`'s `test-parallel-chains.R:54`. Neither touches priors.

### as-cran

`R CMD build` (vignettes built) then
`R CMD check --as-cran --no-multiarch`, `_R_CHECK_CRAN_INCOMING_=false`,
pandoc 3.8.3 on `PATH`.

| check | result |
| --- | --- |
| `checking tests` | **[31m] OK**, `FAIL 0 | WARN 1 | SKIP 98 | PASS 6554` |
| `checking examples` / `--run-donttest` | [76s] OK / [117s] OK |
| `checking re-building of vignette outputs` | **[427s] OK** |
| `checking files in 'vignettes'` / `package vignettes` | OK / OK |
| all seven Rd checks, `code/documentation mismatches`, `Rd contents` | OK |
| Status | 1 ERROR, 1 WARNING, 2 NOTEs |

The four remaining items are the environmental ones this machine always
has, each with a verified cause: `Rdlatex.log` says **"pdflatex is not
available"** (the manual WARNING, the ERROR, and the leftover
`frmtmb-manual.tex` NOTE) and **"package 'V8' unavailable"** (the HTML
NOTE). `Sys.which("pdflatex")` and `Sys.which("texi2dvi")` are both
empty here. Nothing that reads the package's code, documentation,
tests, examples or vignettes is anything but OK.

`_R_CHECK_FORCE_SUGGESTS_=false` was needed: `brokenstick` and
`frmtmb.spline` are Suggests and neither is installed on this machine.
Both uses are guarded (`skip_if_not_installed()`,
`requireNamespace()`), and main has since dropped `frmtmb.spline` from
core's Suggests (`612cde1`), so only `brokenstick` will remain.

The one WARNING inside the test run is
`test-prior-compat.R:505`, "Optimizer did not report convergence:
singular convergence (7)", on the `(x|g) + (z|h)` shape the previous
lane added. It is PRE-EXISTING: the same warning came out of that file
before any edit of mine.

### Merge

Main advanced to `b131fe1` during the run (`612cde1`, `de9d639`,
`b131fe1`). Those three commits touch `.github/`, `DESCRIPTION`,
`dev/build-docs.R`, `docs/**` and `vignettes/case-studies.Rmd`, and
this lane touches **none of them**, not even `NEWS.md`. There is no
overlap with main at all; the usual `NEWS.md` collision will only
appear against whichever sibling lane merges first.

## Where frmtmb's meaning still differs substantively from brms

The complete list after this lane, each row with what it costs. Only
the last two rows are about a prior SPELLING at all; the rest are
about the quantity underneath it.

| # | what | difference | measured consequence |
| --- | --- | --- | --- |
| 1 | a prior is a PENALTY | `frm()` is MAP, not a posterior, so the reported log likelihood, AIC and `anova()` are penalized quantities | the estimates land close where priors are weak (kidney `sd(patient)` 0.38 against brms's posterior mean 0.40) but they are not the same quantity |
| 2 | `class = "sd"` | brms writes the NORMALIZED half-t (density minus one `lccdf` per element); frmtmb evaluates the unfolded density plus the log-Jacobian | exactly `log(2)` per lower-bounded PRIORED element. A constant: it moves no mode. S1 `dT = 2 log 2` |
| 3 | `class = "cor"`, `"cortime"`, `"rescor"` | the same LKJ density on the same correlation, carried onto frmtmb's coordinate with that map's exact Jacobian, which is not Stan's coordinate | `(eta + (d-1)/2) log(1 - rho^2)`, measured -0.0021842 on S3: a coordinate change, not a different density |
| 4 | `class = "ar"`, `"ma"`, `"cosy"` | the same, for the stationarity map | ANALYTIC, not measured against Stan: brms declares `vector<lower=-1,upper=1>[Kar] ar` under `cov = TRUE` and frmtmb uses its own map |
| 5 | `simo` (`dirichlet`) | REFUSED. frmtmb holds a `mo()` simplex as free softmax coordinates and puts no density on it | a `mo()` table stops the call, naming the row |
| 6 | `sds`, `sdgp`, `sdcar` | REFUSED, presentational. The parameters exist as random-effect blocks and are reachable by hand (`class = "sd"`, `group = "s(x)"`), but `as_priorlist()` has no model in hand to read the label from | one named edit round per table |
| 7 | `lscale`, `car` | REFUSED. frmtmb keeps them in the raw internal covariance vector, so they are reachable only on the INTERNAL scale through `class = "theta"` | the density a user writes is not the density they get unless they rewrite it for the internal scale |
| 8 | mixture `theta`/`theta1`/`theta2` | REFUSED as a bare class. A modeled proportion with its own predictor is `class = "Intercept", dpar = "theta1"`, which is honored | brms writes no Stan statement for its reference component either |
| 9 | `horseshoe()`, `R2D2()`, `lasso()` | REFUSED. frmtmb has no hierarchical-shrinkage prior at all | the nearest honest substitute is a tight `normal()` on class `"b"`, which is a different estimator |
| 10 | `constant()`, `dirichlet()`, `uniform()` | REFUSED. `constant()` is brms's way of FIXING a parameter, whose frmtmb counterpart is a constant dpar in `bf()`; `dirichlet` needs a simplex target frmtmb does not have | those tables stop |
| 11 | a distributional parameter's class | **none any more.** Every class brms writes for a dpar is a class `set_prior()` takes, with brms's meaning and brms's model-shape rule | 17 families swept: `sigma`, `nu`, `alpha`, `beta`, `shape`, `phi`, `zi`, `hu`, `kappa`, `ndt`, `quantile`, all honored |
| 12 | `class = "theta"` | frmtmb-only. brms has no spelling for a raw internal covariance parameter | an addition, not a divergence |

Rows 2 to 4 are coordinate or normalization differences that no mode
can see. Rows 5 to 10 are refusals: frmtmb holds the quantity
somewhere else or does not hold it, and every one of them is refused
BY NAME from both routes now, which is the property this lane added to
the ones that were only refused on the translated route.

### Extension suites, against this lane's core

| package | files | result |
| --- | --- | --- |
| frmtmb.eam | 17 of 17 | **1292 pass, 0 fail, 0 error**, `test-surface.R` at **50** (was 47) |
| frmtmb.sample | 10 of 10 | **899 pass, 0 fail**, 1 error, `test-sample-direct.R` at **136** and `test-prior-route.R` at **9** |

The one frmtmb.sample error is `test-parallel-chains.R:54`, the
`local_mocked_bindings()`/pkgload harness artifact described above, not
a prior.

# Punch round, 2026-09-06

Against `dev/reviews/2026-09-05-priors2.md` (verdict PUNCH). The
review's own verification of the placement (4e-15 across five links)
and of the acceptance rule (60 of 60 cells on the two spellings this
round is about) is not repeated here; what follows is the nine items.

## 1. F1, blocker. `get_prior()` no longer names a class `set_prior()` refuses

`R/priors.R:1144-1157`. The own-class branch now skips any dpar whose
class `brms_prior_class_refusal()` rejects, which is the same gate
`set_prior()` uses, so the two cannot disagree by construction.

| model | class column before | after |
| --- | --- | --- |
| `bf(y ~ x) + mixture(gaussian(), gaussian())` | Intercept, b, b, sigma1, sigma2, **theta1** | Intercept, b, b, sigma1, sigma2 |

The row is DROPPED rather than sent back to `Intercept` + `dpar`,
because that spelling is refused on that model too (theta1 has no
predictor): falling back would have restored a quieter version of the
same lie. A mixture proportion frmtmb cannot prior is now absent from
the table, and `?set_prior` says so.

The regression test is `tests/testthat/test-get-prior-route.R:242-317`,
"every row get_prior() lists is one set_prior() accepts". It is
stronger than the review asked for: for **every row** of **seven
shapes** (ordinary, distributional, intercept-only dpar, negbinomial,
zero-inflated, Beta, mixture) it re-spells the row through
`set_prior()` AND resolves it against the fit, so a row that parses
but addresses nothing fails too. That is what makes it a test of the
promise rather than of one column.

## 2. F2, blocker. `class = "b", dpar = X` on an intercept-only predictor is refused

`R/priors.R:1976-1993`, the third case in `dpar_shape_refusal()`.
`dpar ~ 1` gives the parameter a predictor (so it takes the Intercept
spelling) but no population-level slopes (so class `"b"` addresses an
empty set). Measured before: the penalty was bit-identical to no prior
at all, with no message. Measured after: refused, naming
`class = "Intercept", dpar = X`.

Tests: `tests/testthat/test-prior-compat.R:703-716` for the refusal and
for `sigma ~ x` still being accepted, and `:718-737`, a new block that
asserts **brms refuses the same row**
(`brms::validate_prior()` answers "do not correspond to any model
parameter") and accepts it once sigma has slopes. The rule is brms's,
asserted against brms, not restated from the review.

## 3. F3. `frm_sample()` announces the spelling the table lists

`extensions/frmtmb.sample/R/sample.R:724-735`: `lab` is built from
`spec_spelling(s)`. Taken rather than deferred, as the review asked;
the file is this lane's now that wt-sample-ce has merged.

```
before:   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
after:    sigma              student_t(3, 0, 2.5)  [natural scale]
```

which is the word `get_prior(route = "sample")` prints and the word
`set_prior()` accepts.

## 4. F7. `coef` and `group` on a natural class are refused

`R/priors.R:726-763`, `unhonored_coef_refusal()` gains the natural
case, and both call sites pass `group` (`:444`, `:693-694`). The call
already saw the class the user WROTE, so no move was needed. `resp` is
still accepted, because it is the one narrowing such a class takes.
Test: `tests/testthat/test-prior-compat.R:772-786`.

## 5. F4. Which mode, in the migration vignette

`vignettes/brms-migration.Rmd:80-95`, next to the `log(2)` sentence.
Re-measured here rather than copied, because the review's data is not
in the tree; gaussian `y ~ x`, n = 150, `normal(2.5, 0.4)` on sigma:

| quantity | mine | the review's |
| --- | --- | --- |
| frmtmb penalized sigma | **2.133344** | 2.170587 |
| mode of `L * p` in sigma coordinates | **2.126533** | 2.163707 |
| mode in log-sigma coordinates | **2.133358** | 2.170589 |
| unpenalized ML sigma | 2.091048 | - |

Same phenomenon, same gap to three digits (0.0068 against 0.0069), on
different draws. The vignette carries MY pair, because a number in the
repository should be one the repository can reproduce; the review's
pair is recorded here as the independent confirmation.

## 6. F8 and F6, documented in the same pass

`vignettes/brms-migration.Rmd:97-105`. A resp-less natural class
broadcasts across a multivariate model where brms refuses the row:
**kept, not refused**, because it is frmtmb's own convention for every
class-wide prior (`class = "sd"` with no `group` covers every block,
`class = "b"` covers every slope), and making one class family the
exception would be a worse inconsistency than the difference with
brms. Documented instead, with `resp` named as the way to narrow.
`cumulative()` has no `disc` dpar at all, so the acceptance rule is
never reached there; also documented.

## 7. F5. `?frm_simulate` says which scale a dpar column is on

`R/simulate-new.R:497-506`.

## 8. NEWS

`NEWS.md:30-56`: the `get_prior()` bullet gains the dropped mixture
row, and three bullets are added for the `b`/dpar refusal, the
`coef`/`group` refusal, and the changed `print()` / `spec_target()`
text.

## 9. Re-measured: as-cran on a tarball built WITH vignettes

`R CMD build` (no `--no-build-vignettes`), then
`R CMD check --as-cran --no-multiarch`, `_R_CHECK_CRAN_INCOMING_=false`,
`_R_CHECK_FORCE_SUGGESTS_=false`, pandoc on `PATH`, and **no
`--library=`** this time. Status: **1 ERROR, 1 WARNING, 2 NOTEs**, and
every one of the four is named with its cause below.

The review's two vignette WARNINGs are gone, which confirms they were
its own `--no-build-vignettes` artifact: `checking files in
'vignettes'`, `checking package vignettes` and `checking installed
files from 'inst/doc'` are all OK, and `frmtmb.Rcheck/frmtmb/doc/`
holds the eight built `.html` files.

| item | text | measured cause |
| --- | --- | --- |
| `checking PDF version of manual` WARNING | "LaTeX errors when creating PDF version" | `frmtmb.Rcheck/Rdlatex.log` says **"pdflatex is not available"**, twice. `Sys.which("pdflatex")` and `Sys.which("texi2dvi")` are both `""`, and `R_LATEXCMD` is unset |
| `checking PDF version of manual without index` ERROR | "Re-running with no redirection of stdout/stderr" | the same retry, the same missing pdflatex |
| `checking HTML version of manual` NOTE | "Skipping checking math rendering: package 'V8' unavailable" | `requireNamespace("V8")` is FALSE on this machine |
| `checking for non-standard things in the check directory` NOTE | `'frmtmb-manual.tex'` | the leftover of the two items above: `Rd2pdf` writes the `.tex`, `texi2dvi` then fails, so nothing cleans it up. There is **no `frmtmb-manual.pdf`** beside it |

One INFO, not counted in the status: `checking package dependencies`
reports "Package suggested but not available for checking:
'frmtmb.spline'". It is a Suggests of core at this base and is not
installed here; main dropped it from core's Suggests at `612cde1`, so
it disappears on merge. `brokenstick`, the other Suggests, IS
installed and is not reported.

### The reviewer's doubt about pdflatex, settled

The review inferred that pdflatex must be available because
`frmtmb-Ex.pdf` was produced. That file is not a LaTeX product:

```
$ grep -a -o "Producer[^)]*)" frmtmb.Rcheck/frmtmb-Ex.pdf
Producer (R 4.6.1)
$ grep -a -o "Creator[^)]*)"  frmtmb.Rcheck/frmtmb-Ex.pdf
Creator (R)
```

`frmtmb-Ex.pdf` is R's own `pdf()` graphics device output from running
the examples, which needs no TeX at all. The LaTeX product would be
`frmtmb-manual.pdf`, and it does not exist in either check directory,
while `frmtmb-manual.tex` does. Together with `Rdlatex.log`'s explicit
"pdflatex is not available" and the two empty `Sys.which()` results,
that is three independent measurements agreeing. The attribution
stands; the inference from `frmtmb-Ex.pdf` does not.

### One correction to my own earlier report

My first as-cran attempt died on "Packages suggested but not
available: 'brokenstick', 'frmtmb.spline'", and I recorded brokenstick
as uninstalled. That was wrong, and the cause was my own
`--library=<private lib>` on `R CMD check`: it also narrows where
check LOOKS for suggested packages, hiding the user library where
brokenstick lives. The same flag was what killed 11 files of the first
full-suite run mid-flight, because check installs the package into
that library. Dropped here, as instructed.
