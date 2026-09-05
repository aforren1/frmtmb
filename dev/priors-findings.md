# PRIORS lane: making a brms prior specification mean in frmtmb what it
# means in brms

Worktree `C:/Users/adf44/source/r/frmtmb-wt-priors`, branch `wt-priors`,
off `2210aa1` (frmtmb 0.51.0). Written incrementally while the work
happens; the table the maintainer decides from is at the end.

The instruction: follow brms as closely as possible, within reason,
WITHOUT silently changing what an existing frmtmb spelling means, and
report every place where matching brms would be a substantive change
with its measured consequence.

Evidence read first, end to end: `dev/reviews/2026-09-05-brms-priors.md`
(1060 lines, including the punch re-check and its residual R3),
`dev/brms-priors-findings.md` (917 lines), the "Follow-on, priors"
section of `dev/brms-likelihood-tests.md`,
`tests/testthat/test-brms-priors.R` (61 assertions) and
`tests/testthat/helper-brms-priors.R`, `R/priors.R`, `R/fit.R`'s prior
application, `R/sampling-api.R`.

## Environment

- R 4.6.1, private library
  `.../scratchpad/pr-lib`, user library
  `C:/Users/adf44/AppData/Local/R/win-library/4.6` as a READ fallback
  only (brms 2.23.0, rstan 2.32.7, StanHeaders 2.39.1 live there).
- Stan cache: `C:/Users/adf44/source/r/frmtmb/dev/stan-cache` COPIED to
  `C:/Users/adf44/source/r/frmtmb-wt-priors/dev/stan-cache` (53 programs
  plus `makevars-cxx17.mk`); `FRMTMB_STAN_CACHE` points at the copy and
  `R_MAKEVARS_USER` at the copied makevars, whose `CXX17FLAGS` ends in
  `-std=gnu++17`.
- The brms tiers run under `FRMTMB_BRMS_FIT_TESTS=true` and
  `NOT_CRAN=true`, one `test_file()` per process.

## What landed

Four decisions, in the order the sequencing constraint required. Every
file:line below is in this worktree.

### D2. A brms distributional class is ROUTED, not refused

`brms_prior_route()` (`R/priors.R:682-697`) replaces the flat gate.
Nine brms class names are frmtmb's own and pass through
(`brms_direct_prior_classes`, `R/priors.R:593-594`: `b`, `Intercept`,
`sd`, `cor`, `ar`, `ma`, `cosy`, `cortime`, `rescor`. The last five
were refused before, although `set_prior()` spells them identically.
`brms_prior_class_refusal()` (`R/priors.R:604-680`) keeps refusing
`theta*`, `simo`, `sds`, `sdgp`, `lscale`, `sdcar` and `car`, each
naming where the quantity actually lives. Everything else is a distributional parameter
and is routed to `class = "Intercept"`, `dpar = <class>` with the
internal `natural` flag, written ONLY when TRUE
(`as_priorlist()`, `R/priors.R:565-583`).

The placement itself is decided where the link is known, in
`coef_placement()` (`R/priors.R:1655-1672`), which reads the target's
own linear predictor: a log link reuses the exact `scale = "sd"` path
class `"sd"` already had, an identity link needs no map, and every
other link goes through the link object's `linkinv` and `mu_eta`
(`prior_logdens()`, `R/priors.R:1694-1730`). That generality is not
speculative: brms's default on `nu` is `gamma(2, 0.1)` and `nu`'s link
is `logm1`.

**R3 of the punch re-check is fixed and measured.** `internal_bound()`
(`R/priors.R:1674-1692`) carries a bound through the link before it
boxes an internal parameter.

| spec | sigma | max abs gradient |
| --- | --- | --- |
| review's naive D2a routing (bound passed through) | 1.000000000 | 117 |
| this lane, S5 nonlinear | **0.1354697363** | **3.52e-05** |
| the lane's hand-built natural spec (target) | 0.135469736 | 3.5e-05 |

brms writes `lb = 0` on every dispersion default; on a log-linked sigma
that is `log(0) = -Inf`, not a floor of 1.

### D1. A row applies what its prior string says

The `source`-keyed drop and its message are gone
(`as_priorlist()`, `R/priors.R:549-560`). `source` records who BUILT a
row, not who wrote the density in it, and brms does not update it when
a user edits the `prior` cell in place. An empty prior string still
applies nothing; a flat row of a class frmtmb cannot name is skipped
rather than refused, because brms echoes a parameter's declared bounds
onto such rows and a refusal there would be an error the caller did not
ask for.

### D3. class "Intercept" is the intercept brms constrains

`lp_center_offset()` (`R/priors.R:1592-1609`) reads the column means of
the PARAMETRIC design columns, which is brms's own rule, verified off
its generated code rather than assumed: `Xs` (a smooth's unpenalized
part) and `Xmo` sit outside `Xc` and outside `means_X`, and frmtmb
appends those columns after `n_param_cols`. The offset rides on the
entry and is evaluated on the tape by `entry_offset()`
(`R/priors.R:2221-2228`). The map between the two parameterizations is
unit triangular, so it carries no Jacobian, which the measured S7
residual confirms. NO change was needed in `R/frame.R`: the linear
predictor already carries `X` and `n_param_cols`.

### D4. The refusals that misfired, and the rows with no spelling

- the `theta*` hint no longer names `dpar = "theta2"`, a spelling that
  then failed with "Prior target not found"
  (`R/priors.R:608-621`).
- `logistic` and `gamma` parse (`R/priors.R:746-753`), with
  `prior_logistic()` and `prior_gamma()` exported beside the other
  constructors (`R/priors.R:1824`, `:1838`) and AD-safe densities at
  `R/priors.R:1750-1758`.
- brms's ordinal `class = "Intercept"` reaches the thresholds
  (`ordinal_threshold_entry()`, `R/priors.R:1413-1441`), on the
  threshold vector computed from `tau_raw` with the ordered map's
  log-Jacobian and brms's own centering
  (`ordinal_center_offset()`, `R/priors.R:1627-1653`;
  `prior_logdens()`'s `ordthres` arm, `R/priors.R:1696-1710`).
  `get_prior()` now lists the slot (`R/priors.R:965-974`), and the
  resolver's "Prior target not found" names thresholds when that is
  what was asked for (`R/priors.R:1300-1318`).

## The measurements

Every shape below is the tier's, warm from the copied cache.
`frm - stan` is frmtmb's log prior minus the Stan program's at
frmtmb's estimates; `dT` is `log_prob(adjust_transform = TRUE)` minus
frmtmb's penalized objective; `gT` is the AT=TRUE gradient.

| shape | before: frm - stan | after: frm - stan | after: dT | after: gT |
| --- | --- | --- | --- | --- |
| S1 `Reaction ~ Days + (1\|Subject)` | +2.69708527 | +5.64778636 | **2 log 2** | joint |
| S2 `y ~ x + z, sigma ~ x` | -0.00331008 | **0.00000000** | -5.7e-14 | 5.7e-04 |
| S3 `Reaction ~ Days + (Days\|Subject)` | +3.39555038 | +6.16426331 | 2.08016961 | joint |
| S4 `cumulative(y ~ x)` | 0 (no prior carried) | +0.40481353 | **-5.7e-14** | **5.6e-04** |
| S5 `y ~ a exp(-b x)` | 0 (no prior carried) | -2.69215419 | **log 2** | **3.5e-05** |
| S6 `mixture`, `theta1 ~ x` | 0.00000000 (2 rows) | -1.44629885 (5 rows) | 2.74798120 | 1.00 |
| S7 `Reaction ~ Days` | -0.22234198 | +3.16732816 | **log 2** | **6.8e-06** |

On every shape whose optimum is a real mode (S2, S4, S5, S7),
`gT` now vanishes: **frmtmb maximizes exactly the density brms
samples**, up to one `log(2)` per lower-bounded parameter, which is a
constant. S1 and S3 are joint random-effect shapes where Stan's optimum
is not the marginal mode (the flat tier's caveat); their z-block
gradients are 1.7e-14 and 8.4e-15.

### The headline: S7's slope bias is gone

| quantity | no prior | before | after | brms mode |
| --- | --- | --- | --- | --- |
| slope `Days` | 10.46724855 | 10.38302871 | **10.46728412** | 10.46728315 |
| shift, in SEs of the slope | - | **0.0684** | **0.0000289** | 9.8e-06 |
| raw intercept | - | - | 251.35915273 | 251.35914902 |
| centered intercept | - | 298.6627025 | **298.46193126** | 298.461943 |
| sigma | 47.4498179 | 47.48875616 | **47.48791879** | 47.48791992 |

`SE(Days) = 1.231295522`. The previous columns are the measurement
lane's and the reviewer's, reproduced independently.

### The ordinal thresholds, which had no spelling at all

| quantity | value |
| --- | --- |
| frmtmb `tau_raw` | (-1.0609848097, 0.4048135282) |
| thresholds, `(tau_1, tau_1 + cumsum(exp(rest)))` | (-1.0609848, 0.4380378) |
| minus `mean(x) * slope` | (-1.0718448, 0.4271778) |
| Stan's `Intercept` vector | (-1.0718449982, 0.4271779503) |
| `dT` | -5.68e-14 |

Nothing separates them, and it is worth being precise about why,
because the two reasons are different. There is no RENORMALIZER because
an ordinal program has no lower-bounded parameter, so brms writes no
`lccdf` line. There IS a JACOBIAN, and it cancels: Stan's `ordered`
transform contributes `lpT - lpF = 0.4048135282` on this fit, and
frmtmb's `ordthres` entry adds the same `sum(raw[-1])` because it
applies the same map. The cancellation is what makes `dT` zero rather
than a constant, which is a stronger result than the absence of a
renormalizer alone would give.

## The table: every brms prior behavior against frmtmb after this lane

"Substantive difference remaining" is what the maintainer decides on.
Everything else is now the same prior in both packages.

| brms row | frmtmb after this lane | substantive difference remaining | measured consequence |
| --- | --- | --- | --- |
| `class = "b"` (and with `coef`, `dpar`, `nlpar`, `resp`) | same class, same link scale, same coefficients | none | S2's residual is 0 to 1e-13 |
| `class = "Intercept"` | evaluated at the intercept at the predictor MEANS, brms's own | none | S7 slope shift 0.0684 SE before, **0.0000289 SE** after; brms's own is 9.8e-06 SE |
| `class = "Intercept"`, `dpar =` | the same, per sub-formula, as brms's `means_X_sigma` | none | S2 `frm - stan` **exactly 0** |
| `class = "sd"` (whole block, `group =`) | natural sd with the log-Jacobian | brms adds `log(2)` per element that frmtmb does not | a constant; moves no mode. S1 `dT = 2 log 2` |
| `class = "sd"` or `"cor"` narrowed by `coef =` | REFUSED by name, naming the whole-block spelling | none any more | brms writes `exponential_lpdf(sd_1[2] | 1)` and keeps its default on `sd_1[1]`; frmtmb resolves per BLOCK and never read `coef`, so all of `coef = "Intercept"`, `coef = "x"` and no coef gave one bit-identical objective (the reviewer measured 1609.2605567405; reproduced here as -439.6363853199 on another shape). Refusing it is the punch round's C1 |
| `class = "cor"` (`lkj`) | same density on the correlation, carried onto frmtmb's coordinate with that map's exact Jacobian | frmtmb's unconstrained coordinate is not Stan's | `(eta + (d-1)/2) log(1 - rho^2)`, measured -0.0021842 on S3: a coordinate change rather than a different density |
| `class = "sigma"` and every other dpar class (`shape`, `phi`, `nu`, `kappa`, `sigma1`, ...) | ROUTED to `class = "Intercept"`, `dpar =`, marked `natural`: the density is on the parameter itself through its inverse link | none in the translation | S5 sigma **0.1354697363** against brms's mode 0.135469737; `gT = 3.5e-05` |
| the same parameter in **frmtmb's OWN spelling**, `class = "Intercept", dpar = "sigma"` | still LOG sigma; deliberately unchanged | **NARROW.** frmtmb matches brms's meaning for EACH spelling, and brms accepts only one of them on any given model, so no ported script is ambiguous | see below |
| `class = "ar"`, `"ma"`, `"cosy"`, `"cortime"`, `"rescor"` | translated (they were refused before although `set_prior()` spells them identically) | frmtmb's unconstrained coordinate is not Stan's, as for `cor` | ANALYTIC, not measured against Stan: brms declares `vector<lower=-1,upper=1>[Kar] ar` under `cov = TRUE` and frmtmb uses its own stationarity map, so the two differ by that coordinate change |
| ordinal thresholds (`class = "Intercept"` on `cumulative`/`sratio`/`cratio`/`acat`) | the threshold vector at the predictor means, with the ordered map's log-Jacobian | none | S4 `dT` = **-5.7e-14**, `gT` = 5.6e-04. Was: no prior carried at all |
| `class = "simo"` (`dirichlet`) | refused by name | **YES** | frmtmb holds a `mo()` simplex as free softmax coordinates and puts no density on it. A `mo()` model's `get_prior()` table now STOPS the call instead of being dropped in silence |
| `class = "sds"` (a smooth's wiggliness SD) | refused by name, naming `class = "sd", group = "s(x)"` | **YES, presentational** | the parameters exist and are reachable by hand; the automatic translation is not, because `as_priorlist()` has no model to read the smooth's label from |
| `class = "sdgp"`, `"lscale"` (a `gp()`) | refused by name, naming `class = "sd", group = "gp(x)"` and `class = "theta"` | **YES, presentational** for `sdgp`; `lscale` reaches only the INTERNAL scale | same cause |
| `class = "sdcar"`, `"car"` (a `car()`) | refused by name, naming `class = "sd", group = ` and `class = "theta"` | **YES, presentational** for `sdcar`; `car` reaches only the INTERNAL scale | measured: frmtmb lists `sd`/`g` plus `theta_1`, `theta_2` for the same block |
| brms's mixture `theta`/`theta1`/`theta2` (bare class) | refused by name; the modeled proportion's `class = "Intercept", dpar = "theta1"` is honored | none that costs a density | brms writes no Stan statement for the reference component either. The old hint named a spelling that then failed; it does not now |
| flat default (empty `prior` string) | applies nothing | none | brms's own rule. `as_priorlist(prior("", class = "b"))` is `NULL` |
| a flat row carrying brms's declared bound | skipped where frmtmb cannot name the class | none | brms echoes a parameter's declaration onto such rows; refusing one would stop a call the caller did not write |
| `lb`/`ub` on a translated row | carried through the placement's own map | none | brms's `lb = 0` on a log-linked sigma is `log(0) = -Inf`. Naive routing pinned sigma at **1.000000000** with gradient **117**; it is now 0.1354697363 with gradient 3.5e-05 |
| `lb`/`ub` on `class = "b"` | unchanged, link scale on both sides | none | `lower[x] = 0` |
| `beta(1, 1)` on `zi` and `hu` | translated onto `zi`/`hu` ITSELF through the logit link | none | brms's default for EVERY zero-inflated and EVERY hurdle family. It stopped those tables until this punch round added `prior_beta()`. Measured: the entry is `betad[1]`, `scale = "natural"`, link `logit`, and brms's `lb = 0`/`ub = 1` become -Inf/+Inf on the logit scale |
| `inv_gamma(a, b)` on `shape` | translated onto `shape` itself | none | brms's default on a negbinomial `shape`. Also its default on a `gp()` length-scale, where the `lscale` CLASS is refused for its own reasons |
| `horseshoe()`, `R2D2()`, `lasso()` on `class = "b"` | refused, naming the supported densities | **YES** | brms's shrinkage priors, which is what a porting user reaches for on a wide model. frmtmb has no hierarchical-shrinkage prior at all; the nearest honest substitute is a tight `normal()` on class `"b"`, which is not the same estimator. `R2D2()` used to get the generic "Cannot parse prior" because its name is upper-case and its argument list empty; it now reaches the density message |
| other unsupported densities (`uniform`, `constant`, `dirichlet`) | refused, naming the nine supported | **YES** | `constant(1)` is brms's way of FIXING a parameter and has no prior spelling here at all; frmtmb's counterpart is a constant dpar in `bf()`. `dirichlet` needs a simplex target, which frmtmb has none of |
| brms's `tag`, `check` | refused | none | they name things inside a Stan program frmtmb does not build |
| the `source` column | not read | none | this is what made a user's edited row survive |

### The remaining divergence, and why it is narrower than it looks

`set_prior("student_t(3, 0, 2.5)", class = "Intercept", dpar = "sigma")`
means a density on LOG sigma and keeps that meaning. Making it mean
sigma (the review's D2a-prime) would change what an existing frmtmb
script says, and the review measured the cost:

- one hard failure, `tests/testthat/test-simulate-ergonomics.R:120`,
  where a deliberate link-scale prior `normal(log(0.6), 1e-9)` has
  every natural-scale draw negative;
- a SILENT sign flip in a fitted ddm coefficient, `mu.cond` from
  `+0.2654012893` to `-0.3301564534`, which the old assertion did not
  catch.

So it is left alone, and the two spellings sit side by side with
different meanings. `?set_prior`'s "Translating a brms prior" section
and the migration vignette say so, and `test-brms-priors.R`'s row-5
block measures the gap rather than leaving it to prose: the link
spelling puts sigma at 0.1351045271 against brms's 0.1354697363,
capturing about 35 percent of the intended shift on that fit and
effectively none of it on `sleepstudy`, where brms's data-derived
scale is large.

**The review narrowed this, and the correction matters.** The two brms
spellings are mutually exclusive BY MODEL SHAPE, measured off
`make_stancode`:

| model | `class = "sigma"` | `class = "Intercept", dpar = "sigma"` |
| --- | --- | --- |
| `y ~ x`, no sigma predictor | `student_t_lpdf(sigma \| ...)` | brms ERRORS |
| `bf(y ~ x, sigma ~ 1)` | brms ERRORS | `student_t_lpdf(Intercept_sigma \| ...)`, i.e. log sigma |

frmtmb routes `class = "sigma"` to the natural scale and keeps
`class = "Intercept", dpar = "sigma"` on the link scale, which is
**exactly brms's meaning for each spelling, in both directions**. So:

- there is no brms script whose prior is ambiguous after porting;
- the only way to meet the divergence at all is to write frmtmb's
  link spelling on a model with no sigma predictor, which is a spelling
  brms itself rejects;
- and on such a model the link spelling barely moves the parameter (the
  reviewer measured 0.962056457 with no prior, 0.962076230 with the
  link spelling, 0.964013095 with brms's natural one).

That is a documentation question rather than a correctness one, and it
weakens the case for urgency on option 3 below.

**The question for the maintainer.** Three ways to close it, none
taken here:

1. leave it (this lane's choice): two spellings, one documented
   divergence, and no existing script changes meaning.
2. expose `natural = TRUE` on `set_prior()`, so the natural placement
   has a frmtmb spelling as well as a brms one. Additive, and no
   existing meaning changes; costs one argument, its documentation and
   a test.
3. flip the meaning of the existing spelling (D2a-prime). Costs the one
   hard failure and the silent ddm sign flip above, and buys little
   given the mutual exclusivity just measured.

Option 2 is recorded as a low-priority additive follow-up; the review
agrees it is worth doing eventually and is not worth blocking on.

### What is left out, and why

- **`inv_gamma`, `beta`, `dirichlet`, `uniform`, `constant`.** Each is
  one `switch` arm except `dirichlet` (no simplex target) and
  `constant` (not a density at all). `inv_gamma` is the one that still
  costs a real row: it is brms's default on `shape` and on a `gp()`
  length-scale. `inv_gamma` and `beta` were added in the punch round
  (see the section below), which took the stop set from 6 to 3 of 20
  shapes; `dirichlet` and `constant` stay out for the reasons above.
- **`sds`/`sdgp` automatic routing.** `as_priorlist()` runs at the
  argument boundary with no model in hand, so it cannot know a smooth's
  group label. Routing would mean deferring the translation to
  `resolve_priorlist()`, which is a larger change than this lane.
- **brms's mixture ordered intercepts.** brms declares the component
  intercepts as an `ordered` vector; frmtmb does not order them. That
  is a parameterization difference in the likelihood rather than in the
  prior, and the mixture family belongs to another lane. It is why S6's
  `gT` is 1.00 where every other mode-bearing shape's vanishes.
- **`ar`/`ma`/`cosy`/`cortime`/`rescor` identity against Stan.** Routed
  on the strength of the class name being frmtmb's own and the
  placement being the one `?set_prior` already documents. Measuring the
  coordinate change costs new Stan compiles on a shape no tier covers.

### Files touched

`R/priors.R`, `R/simulate-new.R`, `NEWS.md`,
`vignettes/brms-migration.Rmd`, `dev/brms-likelihood-tests.md`,
`man/set_prior.Rd`, `man/get_prior.Rd`, `man/frmtmb-priors.Rd` and
`NAMESPACE` (those four via roxygen),
`tests/testthat/test-brms-priors.R`,
`tests/testthat/helper-brms-priors.R`,
`tests/testthat/test-prior-compat.R`,
`tests/testthat/test-simulate-ergonomics.R`, and this file. No change
was needed in `R/frame.R`, `R/fit.R` or `R/sampling-api.R`: the linear
predictor already carries `X` and `n_param_cols`, which is all the
centering needed.

`SPEC.md:429-431` and `vignettes/inputs.Rmd:194` carried the pre-lane
translation rule and are corrected in this merge (one hunk each; the ce
lane edits `inputs.Rmd` elsewhere). `README.md:59-62` is NOT in that
list and is deliberately unchanged: it never carried that sentence. Its
actual claim, "a prior object brms itself built is translated", was
false on main for a `get_prior()` table and this lane makes it true.

## Verification

One `test_file()` per process, private library `.../pr-lib`, Stan cache
copied into `dev/stan-cache` and warm throughout. Counts audited by
name, not by total.

| run | result |
| --- | --- |
| `test-brms-priors.R` (gated, warm) | **83 passed**, 0 failed, 0 error, 0 skipped, 11 blocks, 33 s, no new compiles. Was 61 assertions, all of them pins |
| `test-brms-likelihood.R` (the log-density tier, gated, warm) | **351 passed**, 0 failed, 0 error, 58 s. The expected count exactly: **no fit moved** |
| `test-prior-compat.R` | 121 passed, 0 failed (was 105; three pins flipped, four blocks added) |
| `test-setprior.R` | 27 passed, 0 failed |
| `test-get-prior-route.R` | 36 passed, 0 failed |
| `test-priors-autocor-classes.R` | 63 passed, 0 failed |
| `test-priors-bounds-grcov.R` | 49 passed, 0 failed |
| `test-simulate-ergonomics.R` | 45 passed, 0 failed (was 41; one block added) |
| `test-message-uniqueness.R` | 6 passed, 0 failed |
| `test-bracket-access.R` | 8 passed, 0 failed |
| `frmtmb.sample` `test-sample-direct.R` | 135 passed, 0 failed. `expect_null(sg$natural)` at `:220` still holds |
| `frmtmb.ddm` `test-surface.R` | 47 passed, 0 failed. No sign flip: the dpar prior there is class `"b"` and is untouched |
| **full core suite**, one file per process | **109 of 109 files, 6148 passed, 0 failed, 0 error, 87 skipped.** Every file present by name, none duplicated. Run twice: an earlier pass gave 6146, the two extra being the `sdcar`/`car` assertions added after it started |
| `roxygen2::roxygenise()` | idempotent: a second pass writes nothing |
| `R CMD check --as-cran` | **Status: OK.** No ERROR, no WARNING, no NOTE. `checking tests ... [406s] OK`, `checking re-building of vignette outputs ... [238s] OK`, `checking examples ... [74s] OK`, `--run-donttest [59s] OK`. Run with `_R_CHECK_CRAN_INCOMING_=false` and pandoc on PATH |

## Punch round, 2026-09-05

The twelve items of `dev/review-priors.md`'s punch list, worked and
measured. The verdict was GO WITH FIXES; item 12 was the reviewer's own
and is kept unchanged.

### 1 (must). `coef =` on `class = "sd"` was silently ignored

CONFIRMED independently before fixing, on my own shape rather than the
reviewer's: `y ~ x + z + (x | g) + (z | h)`, class `"sd"` with
`group = "g"`.

| row | resolved entries | penalized objective |
| --- | --- | --- |
| `coef = "Intercept"` | `theta[1]/sd theta[2]/sd` | -439.6363853199 |
| `coef = "x"` | `theta[1]/sd theta[2]/sd` | -439.6363853199 |
| no `coef` | `theta[1]/sd theta[2]/sd` | -439.6363853199 |

Bit-identical, as the reviewer measured at 1609.2605567405 on his.
brms writes `exponential_lpdf(sd_1[2] | 1)` and keeps its default on
`sd_1[1]`, so honoring the row without its `coef` is a WIDER prior than
the one written, which is the failure D1 exists to remove.

**Fixed by refusing, not by implementing.** `unhonored_coef_refusal()`
(`R/priors.R:642-666`) refuses a translated `sd` or `cor` row carrying
a `coef`, naming the whole-block spelling and the `class = "theta"`
escape hatch; `as_priorlist()` calls it at `R/priors.R:594`. Refusing
rather than implementing per-coefficient `sd` priors keeps the change
inside this lane: the resolver's block model is not this round's
question, and a refusal cannot silently mean the wrong thing.

**The other half of the question, which the fix must not break.** `coef`
on class `"b"` and `group` on `"sd"`/`"cor"` DO narrow. Measured on a
design chosen to discriminate, two slopes and two correlated blocks,
because with one of each every spelling picks the same parameters and a
test would pass while proving nothing:

| spelling | entries | objective |
| --- | --- | --- |
| `b` (no coef) | `beta[2] beta[3]` | -447.6365220974 |
| `b`, `coef = "x"` | `beta[2]` | -444.7695348189 |
| `b`, `coef = "z"` | `beta[3]` | -439.0628495457 |
| `cor` (no group) | `theta[3] theta[6]` | -436.9803411611 |
| `cor`, `group = "g"` | `theta[3]` | -437.2284285786 |
| `cor`, `group = "h"` | `theta[6]` | -436.0504330882 |
| `sd` (no group) | four thetas | -444.3299757510 |
| `sd`, `group = "g"` | `theta[1] theta[2]` | -439.6363853199 |

Pinned at `tests/testthat/test-prior-compat.R:460-527` (the refusal and
the two honored narrowings, in two blocks).

### 2 (must). Two rows the decision table was missing

Both added to the table above: `beta(1, 1)` on `zi` and `hu`, which is
brms's default for EVERY zero-inflated and EVERY hurdle family, and
`horseshoe()` / `R2D2()` / `lasso()` on `class = "b"`. The density row
was also split, because three of the six stops the reviewer measured
were DISTRIBUTION refusals rather than class refusals and the old single
row hid that.

### 3 and 9 (must, was nice to have). `inv_gamma` and `beta` parse

`prior_inv_gamma()` and `prior_beta()` are exported beside the other
constructors (`R/priors.R:1951`, `:1964`), parsed at `R/priors.R:830`
and `:834`, AD-safe densities at `R/priors.R:1847-1859`, samplers at
`R/simulate-new.R:217-222`.

Values and gradients against R's own, tolerance 1e-15:

| density | x | value deviation | gradient deviation |
| --- | --- | --- | --- |
| `inv_gamma(0.4, 0.3)` | 0.15 | 4.441e-16 | 0 |
| `inv_gamma(0.4, 0.3)` | 0.8 | 0 | 0 |
| `inv_gamma(0.4, 0.3)` | 3.2 | 0 | 5.551e-17 |
| `beta(2, 3)` | 0.02, 0.5, 0.97 | 0 | 0 |
| `beta(1, 1)` | 0.02, 0.5, 0.97 | 0 | 0 |

Worst deviation **4.441e-16**, every one under 1e-15. The reference for
the inverse gamma is `dgamma(1/x, shape, rate = scale, log = TRUE) -
2 log(x)`, which is the density by its own change of variables rather
than a restatement of the implementation. Pinned at
`test-prior-compat.R:550-585`.

**The stop set, re-measured over the reviewer's 20 shapes.**

| | before | after |
| --- | --- | --- |
| stops | **6 of 20** | **3 of 20** |
| which | `s()`, `gp()`, `mo()`, negbinomial, zero_inflated_poisson, hurdle_gamma | `s()`, `gp()`, `mo()` |

All three remaining stops are CLASS refusals, structures frmtmb keeps
elsewhere. No shape stops for want of a density any more.

**The `zi` row end to end**, because "it parses" is not the same as "it
lands right": the entry is `betad[1]`, `scale = "natural"`, link
`logit`, kind `beta`, and brms's `lb = 0`/`ub = 1` become `-Inf`/`+Inf`
on the logit scale. That is the non-log branch of the bound transform
exercised in BOTH directions, which is also where the reviewer's
`isTRUE(b == Inf)` repair earns its place. zi moves 0.2614650030 to
0.2635928578. Pinned at `test-prior-compat.R:586-612`.

### 4 (must). Every refused row named in one message

`as_priorlist()` collects refusals and reports them together
(`R/priors.R:552-608`). Measured on `y ~ gp(x)`, the shape the reviewer
found costs two edit rounds:

```
A brms prior table has 2 rows with no faithful frmtmb spelling:
  row 3 (inv_gamma(1.494197, 0.056607), class = "lscale"): ...
  row 4 (student_t(3, 0, 2.5), class = "sdgp"): ...
```

One edit round now, not two. A single bad row still reads as "1 row
with". Pinned at `test-prior-compat.R:529-548`.

### 5 to 8, 10, 11 (should)

| item | what changed |
| --- | --- |
| 5 | `vignettes/brms-migration.Rmd:91-105`: `sdcar` and `car` added to the refused set, `car()` added to the remedy sentence, the `coef` refusal and the 3-of-20 stop set stated |
| 6 | `SPEC.md:426-440` and `vignettes/inputs.Rmd:194` now state the post-lane translation rule. `inputs.Rmd` is ONE hunk on that one line, because the ce lane edits the file elsewhere. `README.md` deliberately unchanged, see below |
| 7 | `helper-brms-priors.R:62-65` comment rewritten; the dead `dropped_as_default` field deleted |
| 8 | `parse_prior_dist()`'s name pattern widened to `[A-Za-z_][A-Za-z_0-9]*` (`R/priors.R:785-793`), so `R2D2()` reaches the unsupported-density message like `horseshoe(1)` and `lasso(1)`. It needed both halves: an upper-case name and an empty argument list |
| 10 | `?frm_simulate` carries the `pars` note (`R/simulate-new.R:487-497`) |
| 11 | row 7 of the table restated with Note B, below |

### 9 (should). The three findings corrections

- **README dropped from the stale-doc list.** Verified: `grep` for
  "faithful" and for "refused by name" in `README.md` returns nothing.
  Its actual sentence is "a prior object brms itself built is
  translated", which was false on main for a `get_prior()` table and
  which this lane makes true. It must NOT change.
- **The ordinal paragraph** now separates the missing renormalizer from
  the Jacobian that CANCELS: Stan's `ordered` transform contributes
  0.4048135282 and frmtmb's `ordthres` entry adds the same, which is
  why `dT` is zero rather than a constant.
- **Row 4 (class `sd`)** was incomplete: it described the whole-block
  case and said nothing about `coef`. Now two rows.
- **Row 19 (unsupported densities)** was incomplete and is now four
  rows, separating what translates from what does not.
- **Row 7 restated.** The divergence is real but narrower than "the one
  to decide" suggested: brms REJECTS `class = "Intercept",
  dpar = "sigma"` on a model with no sigma predictor and REJECTS
  `class = "sigma"` on a model with one, so the two spellings never
  collide in a real brms script, frmtmb matches brms's meaning for each,
  and no ported script is ambiguous.

### Verification after the punch round

Re-run in full, one `test_file()` per process, counts audited by name.

| run | result |
| --- | --- |
| `test-brms-priors.R` (gated) | 83 passed, 0 failed |
| `test-brms-likelihood.R` (gated) | **351 passed**, 0 failed. Still the expected count: no fit moved |
| `test-prior-compat.R` | **161 passed**, 0 failed, 24 blocks (was 121 / 19; five blocks added for C1, the combined message and the two densities) |
| `test-setprior.R` | 27 passed, 0 failed |
| `test-get-prior-route.R` | 36 passed, 0 failed |
| `test-priors-autocor-classes.R` | 63 passed, 0 failed |
| `test-priors-bounds-grcov.R` | 49 passed, 0 failed |
| `test-simulate-ergonomics.R` | 45 passed, 0 failed |
| `test-message-uniqueness.R` | 6 passed, 0 failed |
| `test-bracket-access.R` | 8 passed, 0 failed |
| `frmtmb.ddm` `test-surface.R` | 47 passed, 0 failed |
| **`frmtmb.sample`, all 10 files** | **886 passed, 0 failed**, 2 skipped (`test-sample-direct.R` 135, `test-reparam.R` 260, `test-sampling-ported.R` 207, `test-draws-methods.R` 97, `test-simulators.R` 50, `test-draws-spellings.R` 48, `test-loo.R` 70, `test-prior-route.R` 9, `test-message-uniqueness.R` 6, `test-parallel-chains.R` 4) |
| **full core suite**, one file per process | **109 of 109 files, 6188 passed, 0 failed, 0 error, 87 skipped.** None missing by name |
| `roxygen2::roxygenise()` | idempotent |
| `R CMD check --as-cran` | **Status: OK** |

**One regression of my own, found by as-cran and fixed.** The first
post-punch check returned `Status: 1 WARNING`, a duplicated `scale`
argument entry in `frmtmb-priors.Rd`: `prior_inv_gamma()`'s `scale` was
documented both in the shared `location,scale,df` entry and in a new one
of its own. Folded into the shared entry (`R/priors.R:1877-1879`); the
rerun is `Status: OK` with the Rd usage check OK. Worth recording,
because it is exactly the kind of thing the per-file test runs cannot
see.

### Not taken

`natural = TRUE` on `set_prior()` (the review's option 2). Recorded as
a low-priority additive follow-up: it changes no existing meaning and
the review agrees it is not worth blocking on.

Per-coefficient `sd` priors. Refusing the row is this round's contract;
implementing brms's `sd_1[2]` narrowing means teaching
`resolve_priorlist()`'s block model to address one standard deviation,
which is a change to the resolver rather than to the translation.

## Status

DONE. Nothing was committed, and the main checkout was not touched.
