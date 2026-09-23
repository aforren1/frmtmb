# Seven brms divergences, closed: lane `wt-correct`

Worktree `C:\Users\adf44\source\r\frmtmb-wt-correct`, base `2d70b33`
(frmtmb 0.61.0, frmtmb.sample 0.9.0). The BASE arm of every measurement
reads the round's reference build `C:/Users/adf44/source/r/rellib-r3`;
the LANE arm reads `C:/Users/adf44/source/r/correct-lib` first. brms
2.23.0, bayesplot 1.16.0, mgcv 1.9.4, rstan 2.32.7, R 4.6.1.

**Answer first.** All seven items are done, each against a brms
measurement on the same data. `residuals()` refuses brms's two residual
types on every polytomous family; `pp_check()` on a fit resolves
`group` and `x` against the model frame and behaves as brms does on
every one of the 49 types bayesplot exposes; a class `"sd"` prior
reaches only its own response, distributional parameter and nonlinear
parameter; frmtmb.sample's defaults give a multivariate model brms's
rows (`rescor` `lkj(1)` included) and take the mean offset off the
intercept location; a `mi()` coefficient is `bsp_`; a written `theta`
formula makes the unwritten component the reference; and `?frm` says
what mgcv's REML integrates that frmtmb's does not, with the gap
measured. Nothing else moved: 121 of 126 recorded quantities are
`identical()` across the two builds and the five that are not are an
environment pointer (4) plus the mixture family object the refactor
changes (1). The ported brms bin-1 ledger goes from 243 to 250 passes.

## How to reproduce every number here

| what | script | output under `dev/correct-log/` |
| --- | --- | --- |
| library order, both arms | `dev/correct-prelude.R` | |
| the data every probe uses | `dev/correct-data.R` | |
| what brms does, per item | `dev/correct-brms-probe.R <part>` | `brms-<part>.txt` |
| what frmtmb does, per arm | `dev/correct-probe.R <arm> <part>` | `probe-<arm>-<part>.txt` |
| `pp_check()` on draws | `dev/correct-probe-draws.R <arm>` | `probe-lane-draws.txt` |
| REML against mgcv, 5 seeds | `dev/correct-reml-mgcv.R <arm>` | `reml-mgcv-lane.txt` |
| nothing else moved | `dev/correct-bitwise.R base\|lane\|compare` | `bitwise-*.txt`, `bitwise-compare.txt` |
| one test file, one process | `dev/correct-runtest.R <arm> <pkg> <file>` | `seefail-base-*.txt` |
| install, lane library only | `dev/correct-install.R core sample` | `install-*.txt` |
| the whole ungated suite | `dev/correct-run-suite.ps1` | `suite.log` |
| the gated tier | `dev/correct-run-gated.ps1` | `gated.log` |
| the ported ledger | `dev/brmsport-record.sh`, `dev/brmsport-ledger.R`, `dev/brmsport-gen.R` | `brmsport-*.txt` |
| `R CMD check --as-cran`, both packages | `dev/correct-run-check.ps1` | `check.log` |
| the library matches the sources | `dev/correct-libmatch.R` | |
| the tier counts in this file | `dev/correct-runs.R` | |
| punch round 1, every fix | `dev/correct-punch1.R` | `punch1-lane.txt` |
| punch round 1, the same calls on base | `dev/correct-punch1-base.R` | `punch1-base.txt` |
| the examples NOTE, both arms interleaved | `dev/correct-punch2-extime.R` | `punch2-extime.tsv` |
| which test files leave a `Rplots.pdf` | `dev/correct-punch4-rplots.ps1` | `punch4-rplots.txt` |
| what brms does with `x` on a type that has none | `dev/correct-punch3-brms-xarg.R` | `punch3-brms-xarg.txt` |
| the `simulate(newdata = )` probe | `dev/correct-newdata-feasibility.R` | `newdata-feasibility.txt` |
| the seventh library loss | `dev/correct-loss7-evidence.R`, `dev/correct-libcheck.R` | `loss7-evidence.txt`, `restore-loss7.txt` |

`parts` are `resid`, `ppcheck`, `ppcheck2`, `priors`, `priors2`, `mi`,
`theta`. The brms probe caches its fits in `dev/correct-brmsfits/`,
which is gitignored.

## 1. `residuals()` on a polytomous fit

**brms.** `.predictive_error()` refuses before it looks at the type,
for every `is_polytomous()` family: ordinal, categorical, multinomial,
dirichlet and logistic-normal (`brms/R/predictive_error.R:120`).
Measured (`brms-resid.txt`), on the seed-102 ordinal data and the
seed-103 multinomial data:

```
cumulative   residuals(type = "ordinary")  ERROR: Predictive errors are not defined for ordinal or categorical models.
cumulative   residuals(type = "pearson")   ERROR: (the same)
cumulative   predictive_error(posterior_predict | posterior_epred)  ERROR: (the same)
cumulative   pp_check(type = "error_binned")  ERROR: Type 'error_binned' is not available for polytomous models.
cumulative   pp_check(type = "error_hist")    OK
sratio, categorical: the same
multinomial  residuals and predictive_error: the same; BOTH pp_check
             types ERROR: 'pp_check' is not implemented for this family.
gaussian     every one OK
```

**frmtmb, base.** `residuals()` answered for `cumulative` and `sratio`
under `"response"`, `"ordinary"` and `"pearson"`
(`probe-base-resid.txt`), returning the score residual
`y - sum_k k P(Y = k)`. `categorical` was already refused for every
type. `multinomial` was refused by accident, through "family
'multinomial' declares no mean".

**Changed.** `residuals.frmtmb_fit()` refuses `"response"`,
`"ordinary"` and `"pearson"` for an ordinal or multinomial family with
a `frmtmb_error` naming the family; `categorical` keeps its own
refusal, which covers every type including `"osa"`. `"osa"` still
answers on an ordinal fit: it has no brms counterpart and uses only the
order. `"deviance"` was and is refused. `pp_check(type =
"error_binned")` is refused on a polytomous fit, as brms refuses it.
On draws, `predictive_error()` and therefore `residuals()` refuse the
same families.

**Seen failing on base**: `tests/testthat/test-ordinal-fitted.R` 95
pass, 28 fail, 3 errors on base; 129 pass, 0 fail on the lane
(`punch3-base-test-ordinal-fitted.txt`,
`punch3-lane-test-ordinal-fitted.txt`). The recheck's table restates
every one of these counts and says what moved them. The gated block "residuals
are refused on the shapes brms refuses them on" in
`test-brms-methods.R` asserts brms's refusal and frmtmb's on the same
six shapes.

## 2. `pp_check()` on a fit

**brms.** `pp_check.brmsfit()` looks `group` and `x` up in
`names(model.frame(object))` and hands bayesplot the COLUMN, refuses a
type `bayesplot::available_ppc()` does not list, and requires `group`
for a type whose function has a `group` formal. Measured on all 49
types of bayesplot 1.16.0, on a gaussian, a poisson and a bernoulli fit
(`brms-ppcheck.txt`, `brms-ppcheck2.txt`).

**frmtmb, base** (`probe-base-ppcheck.txt`), on the gaussian design:
15 of the 49 types drew a plot, 12 died on "length(group) must be equal
to the number of observations" and 7 on "is.numeric(x) is not TRUE",
because `group = "g"` and `x = "x"` reached bayesplot as strings; 7 more
are the `loo_*` types and 8 are bayesplot's own refusals of continuous
`y` (bars, rootogram, calibration) or the missing ggfortify;
`type = "violin"` died on "object 'ppc_violin' not found";
`pp_check(fit, "stat_grouped")` with no group died inside ggplot2's
faceting; a group naming a column the model does
not use gave bayesplot's length error rather than brms's "could not be
found in the data".

**Changed.** The fit method now takes brms's arguments in brms's order
(`type`, `ndraws`, `prefix`, `group`, `x`, `newdata`, `resp`, with
`re_formula` after the dots, as the draws method already did), resolves
`group` and `x` against `model.frame(fit)`, refuses a type that is not
a ppc (or ppd) type with the list of valid ones, requires `group` where
bayesplot does, refuses `error_binned` on a polytomous fit, refuses the
`loo_*` types (they weight posterior draws, which an ML fit has none
of), and passes a categorical fit's simulated categories as the 1..K
codes `y` is on, so `bars` and `error_hist` work there as in brms.
`newdata`, `draw_ids`, `nsamples` and `subset` are refused by name:
`simulate()` takes no `newdata`, and the old method swallowed it in its
dots and plotted the fitted rows instead. `resp` is accepted and
IGNORED on a one-response fit, which is what brms does with it there
(punch round 1, minor 2).

`tests/testthat/test-pp-check-types.R` walks every type
`available_ppc()` lists against the verdict brms gave on the same
design, for three responses (147 pairs), and fails rather than skips if
bayesplot grows a type the table does not know. Base: 78 pass, 85 fail,
10 errors. Lane: 183 pass, 0 fail
(`punch3-base-test-pp-check-types.txt`,
`punch3-lane-test-pp-check-types.txt`).

**Not done, on draws.** `pp_check()` on a draws object answers every
type except the seven `loo_*` ones, which die on bayesplot's "One of
'lw' and 'psis_object' must be specified" where brms computes
`loo(save_psis = TRUE)` for them (`probe-lane-draws.txt`). Passing a
PSIS object needs the loo draws and the `posterior_predict()` draws to
be the same subset, which `loo.frmtmb_draws(ndraws =)` and
`posterior_predict(ndraws =)` do not guarantee; filed rather than
guessed.

## 3. The scope of a class `"sd"` prior

**brms.** `prior_re()` writes one class `"sd"` row per PREFIX
(response, distributional parameter, nonlinear parameter), and
`stan_prior()` subsets the prior by that prefix, so an empty field
means the empty one. Measured with `stancode()` (`brms-priors.txt`,
`brms-priors2.txt`):

```
bf(yb ~ x + (1 | g), phi ~ (1 | g)), Beta(), sd group = "g":
  lprior += normal_lpdf(sd_1 | 0, 5)          the mu block
  lprior += student_t_lpdf(sd_2 | 3, 0, 2.5)  phi keeps its default
the same with class sd and no group: the same two lines
mvbind(y1, y2), sd group = "g", no resp:  ERROR: The following priors do
  not correspond to any model parameter: sd_g ~ normal(0, 5)
nonlinear a ~ 1 + (1 | g), class sd, no nlpar:  the same refusal
nonlinear, class sd, nlpar = "a":  lprior += normal_lpdf(sd_1 | 0, 5)
(1 | q | g) in mu and sigma, one block:
  sd group = "g":    normal on the WHOLE sd_1 vector
  sd dpar = "sigma": normal on sd_1[2], student_t on sd_1[1]
```

The last case is brms's own exception: for a block that spans several
prefixes `stan_prior()` reads each prefix's rows with the empty field
allowed too (`brms/R/stan-prior.R:46-73`), and the more specific row
wins by `stan_base_prior()`'s ranking, group first.

**frmtmb, base** (`probe-base-priors.txt`): `sd group = "g"` and class
`sd` both reached theta1 AND theta2, the phi block included; in the
multivariate model a resp-less `sd group = "g"` was accepted and
broadcast to both responses.

**Changed.** `sd_spec_reach()` resolves a class `"sd"` specification
per standard deviation, by the prefix of the column it belongs to
(`block_col_prefix()`, `block_sd_prefix()`), with brms's spanning-block
exception. A specification whose prefix no standard deviation has is
refused, and the message lists the prefixes the model does have.
`prior_table()` lists one class-wide `sd` row per prefix, as brms does,
and `fill_prior_table()` fills only the rows a specification reaches.
Class `"cor"` is unchanged: brms's cor rows carry no prefix (measured:
`cor dpar = "phi"` is refused by brms, and a class-wide `cor` prior
reaches both blocks).

**Seen failing on base**: `tests/testthat/test-prior-sd-scope.R`, 13
pass / 11 fail on base, 24 pass / 0 fail on the lane.

## 4. Two frmtmb.sample defaults

**brms** (`brms-priors.txt`), on brms's own `tests.priors.R` data:

```
bf(mvbind(y1, y2) ~ x + (x | ID1 | g)) + set_rescor(TRUE):
  lkj(1) cor; lkj(1) rescor; student_t(3, 0.3, 2.5) Intercept y1;
  student_t(3, 2, 2.5) Intercept y2; student_t(3, 0, 2.5) sd y1, sd y2,
  sigma y1, sigma y2; b rows flat
y ~ 1 + offset(off), y in {1, 3}, off = 10:
  student_t(3, -8, 2.5) Intercept
y ~ x + offset(o), seed 106:  student_t(3, 1.78834576782975, 2.5)
cnt ~ x + offset(log(expo)), poisson, seed 107:
  student_t(3, -0.900791310755016, 2.5)
sigma ~ z + offset(o2):  sigma's Intercept stays student_t(3, 0, 2.5)
```

**frmtmb.sample, base**: the multivariate model got NO defaults at all
(every row `(flat)`, `rescor` included), and the offset models got
`student_t(3, 2, 2.5)`, `student_t(3, 4.8, 2.5)` and
`student_t(3, 2.2, 2.5)`: the offset was not subtracted.

**Changed.** `default_priors_for()` loops over the responses, reading
each one's location and scale off that response and writing `resp`, and
adds `lkj(1)` on `rescor` when the model has one. The intercept
location is `round(median(y*), 1) - mean(offset)` for the predictor the
intercept belongs to, and only where the response was transformed at
all, which is brms's rule read line for line. A location is now printed
with `as.character()`, not `format()`, so the 15 significant digits
brms's own Stan code carries survive: `format()` gave `1.788346`, a
prior 2.3e-7 away from brms's. The `sd` default is written once per
prefix, which item 3 requires and which leaves every density on every
standard deviation where it was.

**Seen failing on base**:
`extensions/frmtmb.sample/tests/testthat/test-default-priors-brms.R`, 4
pass / 7 fail / 2 errors on base, 15 pass / 0 fail on the lane
(`punch3-base-test-default-priors-brms.txt`). The
last block SAMPLES a multivariate model and checks `prior_summary()`
against the `get_prior(route = "sample")` table row by row.

**What this moves, by design.** Every multivariate model sampled
without `prior = "flat"` (it had no priors and now has brms's); every
model with an `offset()` in its location formula; and any model sampled
with a user class-wide `sd` prior where a distributional parameter also
has a block, because that prior no longer reaches the second block. A
univariate model with no offset and no user `sd` prior is unmoved: the
negative log prior at the estimate is `identical()` across the arms on
five designs (`bitwise-compare.txt`).

## 5. `mi()` is `bsp_`

**brms** (`brms-mi.txt`, seed 108): `variables()` gives `bsp_y_mixm`
and `bsp_y_mixm:z`; `fixef()` rows are `y_Intercept`, `xm_Intercept`,
`y_z`, `xm_z`, `y_mixm` (and `y_mixm:z` after it), so the special-term
coefficients come after every ordinary one and a `mi(x) * z` main
effect comes before its interaction; the prior row stays class `"b"`,
coef `"mixm"`, resp `"y"`; and `hypothesis(fit, "y_mixm > 0")` is
refused, because brms's default `class = "b"` writes `b_y_mixm`, which
no longer exists.

**frmtmb, base**: `b_y_mixm`, `b_y_mixm:z`, `fixef()` rows in design
order (`y_mixm` before `xm_z`, `y_mixm:z` before `y_mixm`), and
`hypothesis(fit, "y_mixm > 0")` answered.

**Changed.** `brms_coef_table()` gives a `mi()` column the `bsp_`
prefix, beside `mo()`, which is the one vector brms holds both in.
`brms_fixef_rows()` sorts the `bsp` rows by interaction order, so a
main effect precedes its interaction as in brms; `mo()` columns are
already in that order in the design, so nothing moves for them.

**Seen failing on base**: `tests/testthat/test-mi-bsp-names.R`, 3 pass
/ 7 fail / 1 error on base, 11 pass / 0 fail on the lane
(`punch3-base-test-mi-bsp-names.txt`,
`punch3-lane-test-mi-bsp-names.txt`). The seventh base failure is minor
4's order assertion.

## 6. A written `theta` formula's reference component

**brms** (`brms-theta.txt`): `stan_mixture()` predicts the mixing
weights only when `length(theta_pred) == nmix - 1`, and the component
with no formula is the reference:

```
K = 2, theta1 ~ x:   vector[N] theta2 = rep_vector(0.0, N)
K = 2, theta2 ~ x:   vector[N] theta1 = rep_vector(0.0, N)
K = 3, theta2 + theta3 ~ x:  theta1 = rep_vector(0.0, N)
K = 3, theta1 ~ x:   ERROR: Can only predict all but one mixing proportion.
K = 2, theta1 + theta2 ~ x:  the same refusal
```

**frmtmb, base**: the family's dpars were always `theta1 ...
theta{K-1}`, so `theta2 ~ x` on two components was "dpar(s) not
available for family", and `theta1 ~ x` on three components fitted with
an intercept-only `theta2` where brms refuses.

**Changed.** `mixture()` builds through `mixture_build(comps, groups,
ref)` and carries a rebuilder; `mixture_theta_reference()` reads the
written theta formulas at parse time, refuses a count other than
`K - 1` with brms's sentence, and rebuilds the family against the
unwritten component. With no theta formula nothing changes: the family
is the same object it was, against the last component, which is where
`brms_coef_table()` reads brms's simplex `theta1 ... thetaK` off it.

The two parameterizations are the same model, which the fit shows:
seed-109 data, `theta2 ~ x` gives `theta2_Intercept -0.2198151`,
`theta2_x 0.7886597`, and `theta1 ~ x` on the same rows gives exactly
the negatives, with the same log-likelihood. brms's posterior mean on
the same data is `-0.2196` and `0.7992`, and `theta2` at `x = 0` is
0.4453 here against brms's 0.4455.

**The gated evidence** is `test-brms-likelihood.R` row 17b: brms's own
`log_prob()` at frmtmb's estimate equals frmtmb's log-likelihood to
`1e-6 * |logLik|`, with brms's gradient there below `1e-3`, for
`bf(y ~ 1, theta2 ~ x)`. That is the same check row 17 makes for
`theta1 ~ x`.

**Seen failing on base**: `tests/testthat/test-mixture-reference.R`, 0
pass / 2 fail / 2 errors on base, 12 pass / 0 fail on the lane
(`punch3-base-test-mixture-reference.txt`). Base cannot go further than
those two errors: it refuses `theta2 ~ x` at model build, so the fit
the rest of the block needs does not exist there.

## 7. REML against mgcv on a location-scale smooth

Documentation plus a measurement, as the item asks: no behavior
changed. `?frm` now says, under `REML`, that mgcv's `method = "REML"`
for a location-scale family integrates the coefficients of every linear
predictor, so the two criteria differ for a smooth in `sigma`.

`dev/correct-reml-mgcv.R`, 5 seeds (44 to 48), n = 500,
`bf(y ~ s(x), sigma ~ s(z))` against
`gam(list(y ~ s(x), ~ s(z)), family = gaulss(b = 0))`, reporting
`log(1 / sd^2) - log(sp)` for each smooth and the largest curve
difference as a fraction of that curve's range
(`reml-mgcv-lane.txt`):

| criterion | mu smoothing parameter | sigma smoothing parameter | mu curve | sigma curve |
|---|---|---|---|---|
| ML | 0.129 to 0.262 | 0.138 to 0.165 | 0.31% to 0.95% | 1.10% to 2.61% |
| REML | 0.005 to 0.011 | 0.142 to 0.170 | 0.029% to 0.055% | 1.30% to 2.93% |

The smoothing-parameter columns are absolute differences; under ML both
are positive and under REML the mu one is negative, which is why the
test compares magnitudes.

Read: REML closes the mu smooth's gap (22 to 28 times smaller) and
leaves the sigma smooth's exactly where ML left it (ratio 1.026 to
1.037), which is what "frmtmb integrates the mu coefficients and mgcv
integrates both predictors'" predicts. The first design tried put the
sigma smooth on its boundary in BOTH packages (mgcv `sp` 6.3e4, edf
1.0006; frmtmb sd 2.1e-5), where the log-sp difference of 10.4 measures
nothing; the recorded design has a wigglier sigma
(`exp(0.8 cos(5 z) - 0.5)`) and both smooths interior.

The test in `tests/testthat/test-smooths.R` pins the RATIOS above, not
the numbers, and passes on base and lane alike, because nothing about
REML changed.

## Punch round 1

The reviewer's one major and eight minors, each with what it now
measures (`dev/correct-punch1.R`, log `punch1-lane.txt`), and the
feasibility question answered rather than built.

**M1, the stale ledger row. FIXED.** `brmsfit-methods:837` still read
`defect / accepts-refused` with the reason "residuals() on the ordinal
fit4 RETURNS the response code minus the expected score", which item 1
made false: the row's recorded message is now this package's refusal,
and brms's assertion (`expect_error(residuals(fit4), "Predictive errors
are not defined")`) fails only on the WORDING. Under the user's rule of
2026-09-17 that is an own-words pass, so the row moved from
`dev/brmsport-verdicts-manual-fit.tsv` to `dev/brmsport-verdicts-own.tsv`
with the pattern `not defined for the .sratio. family, as in brms`. The
tier was re-recorded on the punch build and the ledger rebuilt:
**pass 243 to 250, defect 71 to 64**. The reviewer's reading of 675,
682 and 685 is confirmed: each fails for the reason its verdict now
states, and 682's reason was rewritten in this lane to name the newdata
refusal above it rather than the group-name defect it used to describe.

**Minor 1, the multinomial prose. FIXED on the second pass.** brms
refuses `pp_check()` outright for multinomial ("'pp_check' is not
implemented for this family."), so neither `error_binned` nor
`error_hist` reaches bayesplot there; only the ordinal and categorical
fits answer `error_hist`. Section 1's block says so now. **The first
pass claimed this and did not do it**: the line it describes still read
"sratio, categorical, multinomial: the same", and multinomial now has a
row of its own. The behavior agrees on both sides, which is measured
rather than assumed: on the lane build a multinomial fit refuses
`error_binned` with the polytomous sentence and `error_hist` with
"pp_check() on a fit supports vector responses" (`punch1-lane.txt`,
last block).

**Minor 2, `resp` on a one-response fit. MATCHED.** brms's
`validate_resp()` returns `NULL` when the model has one response, so
brms plots whatever `resp` names, a name the model does not have
included. The lane refused. `resp` is now accepted and ignored on a
univariate fit: `pp_check(fit, resp = "w")` and the bare call give
`identical()` plot data at one seed. A multivariate fit is still
refused by `single_response()` before `resp` could select. The DRAWS
method had the same stricter refusal, and minor 3's rule (one case, one
answer) applies to it too, so it now ignores `resp` on a one-response
model as well and keeps its refusal for a multivariate one.

**Minor 3, two messages for one case. FIXED.** `pp_check(type =
"error_binned")` on an ordinal model now gives the same sentence on
draws as on a fit: "Type 'error_binned' is not available for polytomous
models: the 'cumulative' response is a category, ...". Measured on both
in `punch1-lane.txt`.

**Minor 4, `variables()` order. FIXED, by sorting the terms rather than
the names.** `mi()` terms now pass through
`mo_terms_in_brms_order()`, the sorter `mo()` has always used, so the
design columns are in brms's terms() order and everything downstream
follows: `variables()` reads `bsp_y_mixm` before `bsp_y_mixm:z`, and
`fixef()`/`vcov()`/`summary()` are unchanged.

**Whether the 2026-09-17 decision covers this, worked out rather than
assumed.** That decision is "`variables()` keeps frmtmb's order"
(`dev/round-handoff.md:106`), and what it stands on is
`dev/brmsnames-findings.md`, "Found and NOT fixed" item 1: brms stores
every LEVEL of one group-level coefficient before the next, and brms
orders population-level coefficients BY CLASS, while frmtmb.sample
keeps the template order that every positional read of `x$draws`
relies on. Neither half is this case. `bsp_y_mixm` and `bsp_y_mixm:z`
are two terms of the SAME class, and the fix does not re-sort
`variables()` at all: it sorts the design columns, so the parameter
vector itself holds the main effect first and `variables()` reports it
in the order it has always reported, which is what the decision
protects. Nothing positional moves, because the two columns move
together everywhere.

**Minor 5, smoothing standard deviations. SAID.** On
`bf(y ~ s(x), sigma ~ s(z))` a bare class `"sd"` prior reached both
smooths on base and reaches the `mu` one only now (theta 1 of 1, 2),
and `default_prior()` grows the `dpar = "sigma"` rows; the sampling
route gives each smooth prefix its own `student_t(3, 0, 2.5)`, so no
smoothing standard deviation loses a default. Both NEWS entries name
smooths now.

**Minor 6, the quoted type. FIXED.** `residuals(fit, type = "ordinary")`
says `residuals(type = "ordinary")`; `"response"` and `"pearson"` say
their own.

**Minor 7, the heading. NARROWED.** The reproduction is a poisson fit;
the gaussian one of the same shape answers. "Found and NOT fixed" item
1 says so.

**Minor 8, the missing warning. MATCHED.** brms hands bayesplot
`data[[group]]` whatever the type, so a type with no `group` argument
gets the COLUMN in its dots and bayesplot warns "The following
arguments were unrecognized and ignored: group"; a name that does not
resolve gives `NULL`, which drops out and warns nothing. Both halves
are copied and asserted, and the test's comment no longer claims
parity of the whole call where only the plot agreed.

### The second pass

The first pass was stopped part way, so its tiers ran against a build
its last fixes had not reached, and it wrote one claim it had not
carried out (minor 1 above). What the second pass ran, all of it on
**RTMB 2.0** and against `C:/Users/adf44/source/r/correct2-lib`, which
`dev/correct-libmatch.R` reports as holding the worktree sources
exactly (1118 of 1118 and 209 of 209):

| tier | why it had to be rerun | result |
|---|---|---|
| the ported tier (`dev/brmsport-record.sh`) | `dev/brmsport-gen.R` rewrote all 14 files AFTER the last recorded run | 15 files, 631 pass, 0 fail, 0 error (`punch2-record.txt`) |
| the ledger (`dev/brmsport-ledger.R`) | to confirm M1's totals against that rerun | byte-identical to the first pass, 250 pass and 64 defect, and it stopped on nothing (`punch2-ledger.txt`) |
| `dev/brmsport-gen.R` | the ledger is its only input, and that input did not move | 14 files, 494 ids, 0 without a verdict |
| the bitwise battery, BOTH arms | the recorded arms predate the punch fixes | 121 of 126, with the control above |
| the bitwise lane arm again | the recheck changed `R/conditional-effects.R` | 121 of 126, the same five |
| the whole ungated suite | the first pass's run was cut off inside `frmtmb.eam` | in the generated block below |
| `R CMD check --as-cran` | the passing check ran from 18:41 to 19:50, BEFORE the punch edits to `R/predict.R`, `R/conditional-effects.R`, `R/frame.R` and `methods-draws.R` | in the generated block below |
| the gated tier | the same reason: it ran at 17:02, before the punch edits | in the generated block below |
| `R CMD check` a second time | two assertions were added after the first recheck built its tarball | in the generated block below |

Two pins were added, because the first pass changed behavior that no
assertion held. Minor 6's refusal wording is now asserted per type, on
four ordinal families and on multinomial, which is 12 more failures on
base. Minor 2's and minor 8's four assertions moved into a test_that
block of their OWN: in the block that held them they never ran against
0.61.0, because that block's `type = "violin"` assertion errors there
and testthat ends the block at the error. Split out, they fail on base,
which is the point of writing them.

**The examples-timing NOTE, and why it is the box.** frmtmb checked
with 2 NOTEs on two runs and with 1 on a third, all of this same tree.
The one that comes and goes is the examples-timing NOTE, and it is not
the lane's code.
`dev/correct-punch2-extime.R` runs the `?residuals.frmtmb_fit` example
verbatim, one arm per process, alternating base and lane over six
rounds so a load excursion hits both, and carries a fixed arithmetic
CONTROL that reaches no frmtmb code. Minimum of six
(`punch2-extime.tsv`):

| minimum of six | base | lane | lane/base |
|---|---|---|---|
| the example | 8.51 s | 6.94 s | 0.816 |
| the control | 1.24 s | 1.23 s | 0.992 |

**What this measures is the LEVEL, not a difference between the arms.**
The raw times in one interleaved sequence run 6.94, 8.63, 8.80, 8.98,
12.98, 13.50, 15.00, 15.11, 15.64, 16.50, 18.26 seconds on identical
code, so the 5 second threshold is crossed by load: **the RELEASED
build's own example takes 8.51 to 18.26 seconds on this box**, well
over the threshold, and raises the NOTE by itself. Another lane's suite
was running throughout, with the processor at 71 percent.

**The lane/base ratio itself does not reproduce, and that is recorded
rather than deleted.** The 0.816 above reads as the lane being faster.
An independent set of nine interleaved rounds in the recheck read
**1.072**, the opposite direction, with its control at 1.048 against
0.992 here. So at these counts the instrument does not resolve a
difference between the arms; alternation removes staleness and ordering
bias, which is what it is for, and it does not remove variance. **No
claim about which arm is faster is made, and none is needed**: what the
NOTE turns on is the level, and both sets agree that the level is over
the threshold on both arms.

**A third check then agreed, from the other direction.** The examples
phase of this tree has now been checked three times and taken 108 s
(NOTE), 96 s (NOTE) and **39 s (OK)**. The final check is the 39 s one,
so the generated block below records `Status: 1 NOTE`, and that NOTE is
the V8 one every release records. Two instruments, the interleaved
timing and the repeated check, support the same conclusion, which is
about the level and not about the arms: identical code crosses and
recrosses the threshold as the box's load moves.

`dev/correct-punch1.R` also grew the two measurements the first pass
did not take: `resp` on a one-response DRAWS object, which is minor 2's
other half and was edited after that pass's probe ran, and what a
multinomial fit answers, which is minor 1.

**The base arm of minors 2, 6 and 8, constructed rather than counted.**
`tests/testthat/test-pp-check-types.R` cannot show what 0.61.0 does
with these three calls: its "pp_check() refuses as brms refuses" block
ERRORS on base at the `type = "violin"` assertion, which ends the block
before the three calls are reached, so on base they never run. A count
from that file would say nothing about them. `dev/correct-punch1-base.R`
runs them on the reference build instead (`punch1-base.txt`), beside
the lane's `punch1-lane.txt`:

| call | 0.61.0 | the lane | brms |
|---|---|---|---|
| `residuals(fit, type = "ordinary")` on `cumulative()` | answers, 0.9457 -0.2 0.312 | refuses, quoting `type = "ordinary"` | refuses |
| `pp_check(fit, resp = "w")`, one response | `resp` is not a formal; it falls into the dots and bayesplot warns "unrecognized and ignored: resp" | accepted and ignored, silent | accepted and ignored, silent |
| `pp_check(fit, type = "dens_overlay", group = "g")` | warns "unrecognized and ignored: group" | the same | the same |
| `pp_check(fit, type = "dens_overlay", group = "nosuch")` | warns, because the STRING reaches bayesplot | silent, because a name the model does not use is not looked up for a type with no `group` | silent |

The base arm of the three test files the punch round changed was
re-taken twice: once as the first pass left them, and again after the
two pins above, which is where the pins prove themselves.

| file | base, before the pins | base, after |
|---|---|---|
| `test-pp-check-types.R` | 74 pass / 80 fail / 4 errors | 76 / **82** / 4 |
| `test-ordinal-fitted.R` | 95 / 16 / 1 | 95 / **28** / 1 |
| `test-mi-bsp-names.R` | 2 / 7 / 1 | not changed here |

The two extra failures in the first row are minor 2's silent `resp` and
minor 8's silent unused `group`, which is what the block split bought:
neither ran on base before it. The twelve in the second are minor 6's
wording, one per family and type. The logs are
`punch2-seefail-base-*.txt` for the first column and `punch2-base-*.txt`
for the second. **The recheck below split more blocks and moved these
counts again; its table is the current one.**

### The recheck of punch round 1

The recheck found no code defect and two RECORD defects. Both are
closed here, and the four minors with them.

**A1. The suite headline did not measure this tree, and the record said
it did.** `dev/correct-log/suite.log` finished at 06:24:10;
`test-ordinal-fitted.R` was written at 06:24:42 and
`test-pp-check-types.R` at 06:26:14, both after it. The log's rows read
`pass=114` and `pass=180`, where the tree at the time gave 129 and 181.
So "274 files / 16,194 pass" was a number from before the last edit,
which is not a number. The suite was rerun; the new total is in the
generated block below.

**The same question asked of every other generated number**, by
comparing each log's mtime against the mtime of what it covers. This is
the check that would have caught A1, so it is written down rather than
done once:

| number | its log | covers | verdict |
|---|---|---|---|
| suite | `suite.log` | all 274 test files | **was stale, rerun** |
| gated | `gated.log` | 38 gated files, none of them edited in either pass | was sound, rerun anyway because `R/` changed |
| check | `check.log` | the built tarball | rerun, because `R/conditional-effects.R` changed at the recheck |
| in-check counts | `<pkg>-testthat.Rout` | the same run | **had no source at all, see A2** |
| ledger | `brmsport-ledger.tsv` | the 15 recorded tier files | sound: recorded 05:23:15, ledger 05:23:25, and `brmsport-gen.R` leaves all 15 byte-identical |
| bitwise | `bitwise-compare.txt` | both arms | sound: the lane arm was re-recorded at 09:32:46 after the recheck's `R/` change, the base arm stands at 05:24:37 because the reference build did not move, and the compare at 09:32:47 is after both |
| the examples timing | `punch2-extime.tsv` | both installed builds | sound: the arms are interleaved in one sequence, so neither can be stale against the other |
| minor 6's see-fail | `punch2-base-test-ordinal-fitted.txt` 06:25:46 | the file at 06:24:42 | sound |
| minor 2's see-fail | `punch2-base-test-pp-check-types.txt` 06:27:03 | the file at 06:26:14 | sound |

**A2. The in-check counts had no surviving source.** Closed in the
"Runs" section below, and the driver now keeps `testthat.Rout`.

**Minor 1, the `x` half of the unrecognized-argument rule. MATCHED, and
brms was measured rather than read.** brms is ASYMMETRIC here, and only
the `group` half had been copied. `dev/correct-punch3-brms-xarg.R`
measures all four cases on brms 2.23.0
(`punch3-brms-xarg.txt`):

| call, on a type with neither formal | brms |
|---|---|
| `group = "g"`, a name the data carries | WARNS |
| `group = "nosuch"` | silent |
| `x = "x"`, a name the data carries | WARNS |
| `x = "nosuch"` | **WARNS** |

brms writes `ppc_args$group <- data[[group]]` and
`ppc_args$x <- as.numeric(data[[x]])`, so an unknown `group` is NULL and
drops out of the argument list while an unknown `x` is `numeric(0)`,
which stays. Controls, on a type that DOES take the argument: both
names are refused with "could not be found in the data". The lane's
`use` vector now reads `x = !is.null(x)`, which reproduces all four,
and both halves are asserted.

**The pin for it needed a mutant, because base cannot fail it.** 0.61.0
has no `x` formal at all, so `x` falls into its dots and bayesplot warns
in BOTH cases; the assertion passes there by accident. The build it
distinguishes is this lane's own, before the one-line fix. That build
was rebuilt as a mutant (the one line reverted, installed to a throwaway
library, `punch3-mut-install.txt`) and the file run against it
(`punch3-mut-xarg.txt`): **182 pass, 1 fail**, and the one failure is
`x = "nosuch"` at line 238. Exactly one assertion moves, and it is the
right one. The mutant library was deleted after the run.

**Minor 2, the dead tails, swept rather than patched.** testthat ENDS a
block at an error, so every assertion behind one is skipped and a count
from that block says nothing about them. Four files had such tails and
are now split one behavior per block:

| file | what hid them | now |
|---|---|---|
| `test-pp-check-types.R` | `type = "violin"` errors on base, hiding 8 | one block per refusal, 8 blocks from a table |
| `test-pp-check-types.R` | the first grouped `pp_check()` call errors, hiding 3 | three blocks, one per type |
| `test-ordinal-fitted.R` | the first of three types errors, hiding the rest | one block per type |
| `test-mi-bsp-names.R` | `hypothesis()` errors with a message of its own | three blocks: the names, `hypothesis()`, the prior row |
| `test-default-priors-brms.R` | the missing `sd phi` row errors on lookup, hiding 5 | two blocks: the per-prefix rows, and what a user `sd` prior reaches |

**Where splitting would buy nothing, said rather than done.**
`test-mixture-reference.R` cannot be proved on base BY CONSTRUCTION:
`theta2 ~ x` is refused at model build there, "dpar(s) not available for
family 'mixture(gaussian, gaussian)': theta2"
(`punch3-base-test-mixture-reference.txt`), so the fit the block needs
does not exist and splitting produces the same error more times. Its
base arm stays 0 pass / 2 fail / 2 errors.

The recheck grouped `test-default-priors-brms.R` with it as unprovable
by construction. **That is not what it measures.** Its models build on
base; what fails is downstream. Splitting its `sd` block therefore
bought real coverage, 2 more passes and 1 more failure on base. Only
its last block is unreachable, for a different reason: on base a
multivariate model gets no priors at all, so `prior_summary()` has no
rows and `ps$prior` is NULL, and every assertion after that reads a
table that does not exist. That is the defect itself, not a test
problem.

**Every see-fail count, restated against the current tree.** Each file
one process, both arms, on RTMB 2.0:

| file | base | lane |
|---|---|---|
| `test-pp-check-types.R` | 78 pass / 85 fail / 10 errors | 183 / 0 / 0 |
| `test-ordinal-fitted.R` | 95 / 28 / 3 | 129 / 0 / 0 |
| `test-mi-bsp-names.R` | 3 / 7 / 1 | 11 / 0 / 0 |
| `test-mixture-reference.R` | 0 / 2 / 2 (unprovable past the build) | 12 / 0 / 0 |
| `test-default-priors-brms.R` | 4 / 7 / 2 | 15 / 0 / 0 |

The logs are `punch3-base-*.txt` and `punch3-lane-*.txt`. The error
counts rose where blocks were split, which is the point: each error is
now one behavior rather than a lid over several.

**Minor 4, the shipped plot file.**
`extensions/frmtmb.sample/tests/testthat/Rplots.pdf` is left by a test
that draws to the default device. `.gitignore` hides it, so it is
invisible in `git status`, and `R CMD build` would have shipped it. The
ignore rule is fixed in both packages: core's `^Rplots\.pdf$` anchored
at the package root and so never matched `tests/testthat/`, and
frmtmb.sample had no rule at all. Both now carry `Rplots\.pdf$`.

**Deleting the file is not the fix**: it came back during the reruns,
because a test recreates it. What keeps it out of the tarball is the
ignore rule, so the rule was measured rather than trusted, with the
file PRESENT in `tests/testthat/` both times:

| `.Rbuildignore` | `tests/testthat/Rplots.pdf` in the tarball |
|---|---|
| with `Rplots\.pdf$` | **absent**, against 37 files kept under `tests/testthat/` |
| with the line deleted | **present** |

The second row is the absent-guard case this project requires, and it
is what says the rule is doing the work rather than something else.

**Which files actually draw, measured, because the first list was a
grep and the grep was wrong.** `dev/correct-punch4-rplots.ps1` deletes
`Rplots.pdf`, runs each of the 31 frmtmb.sample test files in its own
process with the GATED environment on, which is the superset, and
records whether the file came back (`punch4-rplots.txt`, 31 of 31
produced a RESULT line):

| file | leaves `Rplots.pdf` |
|---|---|
| `test-brms-suite-methods.R` | **yes** |
| `test-draws-methods.R` | **yes** |
| the other 29 | no |

**Two of 31, not the six this file named a moment ago.** The grep that
produced that list had five false positives and one false negative:

- `test-conditional-effects-draws.R`, `test-draws-spellings.R`,
  `test-generic-collision.R`, `test-loo.R` and `test-sampling-ported.R`
  all run CLEAN. A drawing call that is written is not a drawing call
  that runs: `test-loo.R`'s only one sits inside `expect_error()`, so
  it throws before it can draw.
- `test-draws-methods.R` DOES draw and the grep excluded it, because
  the file contains a `grDevices::pdf(NULL)` and the grep was
  whole-file. That guard is at line 504 and covers the
  `conditional_effects()` block only; what draws is the `pairs()` call
  at line 361, outside it. A per-file grep cannot see that a guard is
  partial.

So the device-guard fix, when it is taken, is two files and one of them
is a PARTIAL guard rather than a missing one. It is not taken here
because `test-brms-suite-methods.R` is GENERATED from the ledger by
`dev/brmsport-gen.R` and must not be hand-edited, so half of it needs
the generator changed, which is outside this round.

### The feasibility question: `simulate(newdata = )`

Reported, not built (`dev/correct-newdata-feasibility.R`, log
`newdata-feasibility.txt`).

**What it would take: about 25 lines in `simulate.frmtmb_fit()`, and no
new machinery.** frmtmb.sample's `posterior_predict.frmtmb_draws()`
already draws at `newdata`, and every piece it uses is EXPORTED by
core's sampling API: `frm_linpred(newdata =)` for the distributional
parameters, `aterms_for_newdata()` for `trunc()` bounds, `has_trunc()`,
`sim_is_structured()` for the refusal, and
`sim_context()`/`sim_draw()`/`fit_extras()` for the draw itself. The
probe assembles exactly that OUTSIDE the package and draws 10
replicates at a 12-row `newdata`.

**Whether the ported rows would pass: yes, both.**
`expect_ggplot(pp_check(fit1, newdata = fit1$data[1:10, ]))`
(`brmsfit-methods:675`) and the `violin_grouped` call with `group` read
off the newdata (`:682`) both build a ggplot in the probe. `pp_check()`
would also have to read `group`/`x` off `newdata` when one is given, as
brms's `current_data()` does; that is two lines in the method.

**What it would cost.** Three things, and the third is the one to
measure before shipping:
1. A structured family must refuse: a group-level mixture's draw walks
   the fitted groups, and the probe's refusal is the sentence the draws
   method already uses ("structured draw: the structure indexes the
   fitted rows").
2. `simulate()` gains a formal, which is additive to
   `stats::simulate()`, and three callers see it: `frm_bootstrap()`,
   `pp_check()` and `dharma_residuals()`. Only `pp_check()` would pass
   it.
3. `re_formula = NA` is the one piece that is not a copy. On the fitted
   rows `simulate(re_formula = NA)` draws FRESH group effects per
   replicate (`draw_b()`), while `frm_linpred(re_formula = NA)` drops
   them; at newdata the two must be composed, by drawing `b` and
   evaluating the newdata design at that draw. The probe does compose
   them, and the composition is right: over 2000 replicates the mean
   within-group correlation of two rows is 0.105 against the
   `sd_g^2 / (sd_g^2 + sigma^2)` of 0.116 that a shared draw implies,
   where dropping the effects gives -0.013. **At 200 replicates this
   same probe read 0.027 and I wrote down that the composition did not
   work**; the standard error there is 0.07, so the instrument could
   not tell 0 from 0.116. The dissolved reading is recorded in the
   script, because the next reader will want to know which number to
   trust.

So the recommendation is that it IS worth building, in a lane that can
own `simulate()`'s argument surface: it removes a refusal this lane had
to invent, returns `:675` to a plain pass and `:682` with it, and the
only new thinking is the `re_formula = NA` composition above, which
wants its own test at a replicate count that can see it.

## Nothing else moved

`dev/correct-bitwise.R` fits six designs (gaussian, poisson, a
distributional model, a mixed model, a mixed poisson, a distributional
mixed model) under ML, `REML = TRUE` and
`frmtmb_control(profile = TRUE)`, and five more the lane's code paths
touch (a mixture with no theta formula, a mixture with `theta1 ~ x`, a
`mo()` fit, an ordinal fit, a MAP fit with a class `sd` prior), and
records `vcov()`, `fixef()`, `logLik()`, `summary()`, `predict()` at
seed 7, `variables()`, `residuals()` (or `residuals(type = "osa")` and
`fitted()` on the ordinal one), and the negative log prior of
frmtmb.sample's defaults at the estimate on five designs.

**121 of 126 quantities are `identical()`** across base and lane
(`bitwise-compare.txt`). The five that are not are all `summary()`:
three differ only in the `.Environment` of the stored formula, which is
a pointer and never survives two processes (`all.equal()` TRUE), and
the mixture pair differ in that pointer and in the `family` object,
which is the refactor: the family now carries `mix$ref` and
`mix_rebuild`. Their `vcov()`, `fixef()`, `logLik()`, `variables()` and
`residuals()` are identical.

**Punch round 1 moved none of it, and the pointer is a process
artifact.** Both arms were re-recorded after the punch fixes, on RTMB
2.0, and four comparisons are in
`dev/correct-log/punch2-bitwise-deltas.txt`:

| comparison | identical | what differs |
|---|---|---|
| base against lane, before the punch | 121 of 126 | the five above |
| base against lane, after the punch | 121 of 126 | the same five |
| lane before against lane after | 121 of 126 | the formula pointer only |
| base before against base after | 121 of 126 | the formula pointer only |

The last row is the CONTROL: the reference build did not change, so its
five differences are what two R processes always produce, and the third
row says the punch fixes produce nothing beyond that.

Two accessors ERROR identically on both builds and are recorded as
their message: `vcov()` and `summary()` on a poisson fit with no random
effects under `REML = TRUE` or `profile = TRUE` die in
`vcov_estimated()` with "length of 'dimnames' [1] not equal to array
extent". Pre-existing, found by this harness, filed below.

## Found and NOT fixed

1. **`vcov()` dies on a POISSON fit with no random effects under
   `REML = TRUE` or `profile = TRUE`.** `frm(bf(cnt ~ x) + poisson(),
   REML = TRUE)` then `vcov(fit)`: "length of 'dimnames' [1] not equal
   to array extent", from `vcov_estimated()`. The gaussian fit of the
   same shape ANSWERS, so the heading is the poisson one and how much
   wider it goes was not measured. Identical on base, so pre-existing,
   and loud. Reproduce with `dev/correct-bitwise.R`.
2. **`pp_check()`'s `loo_*` types on draws** die on bayesplot's "One of
   'lw' and 'psis_object' must be specified" where brms computes the
   PSIS object itself (section 2).
3. **A multivariate class `"b"`, `"Intercept"` or `"sigma"` prior with
   no `resp` is accepted and broadcast**, where brms refuses it exactly
   as it refuses the `sd` one ("The following priors do not correspond
   to any model parameter: b ~ normal(0, 5)", measured in
   `brms-priors.txt`). Only the `sd` half was in this lane's brief. The
   same rule and the same helper would close it; it breaks more call
   sites, so it needs the user's word rather than a lane's.
4. **Class `"cor"` accepts a `dpar`** (`cor dpar = "phi"` reaches the
   phi block), where brms refuses it: brms's cor rows carry no prefix.
   Measured in `brms-priors2.txt`. Left alone because narrowing a
   correlation prior to one predictor is a capability, not a wrong
   answer, and the item named `sd`.
5. **`residuals(type = "osa")` on a multinomial fit** dies with TMB's
   "'observation.name' must be in data component", a raw internal
   message. brms has no counterpart (it refuses every residual there),
   so this lane's refusal does not reach it.
6. **`pp_check(type = "error_binned")` on a fit bins the SIMULATED
   responses**, where brms bins `posterior_epred()` draws. The plot is
   noisier than brms's, not wrong; the frequentist analog of an epred
   draw under `re_formula = NA` is a separate question.
7. **The `error_binned` refusal calls a multinomial response "a
   category".** It is a set of counts OVER categories, which the
   `residuals()` refusal beside it says correctly. One word in one
   message, found by punch round 1's own probe (`punch1-lane.txt`).
   Left alone because a message change needs its assertion rerun and
   this round's brief named the eight items only.

## What needs the user

Nothing blocks. Two decisions are recorded rather than taken:

1. **Item 3 breaks more than its own item.** In a multivariate model
   brms refuses every prefix-less prior, not only `sd`. Finding 3 above
   is that other half.
2. **The version numbers.** Both NEWS sections are under
   "(development version)"; the `sd` scope, the `mi()` name, the
   mixture reference, the `residuals()` refusal and the `pp_check()`
   argument order are BREAKING, and frmtmb.sample's defaults move
   multivariate and offset sampling.

**For release day.** The examples-timing NOTE on this box is load
dependent, so the release check may report 1 NOTE or 2 on IDENTICAL
code. This tree checked at 108 s (NOTE), 96 s (NOTE) and 39 s (OK).
Compare the example timings against the level this file records before
reading a second NOTE as a regression.

## The seventh library loss, in the middle of this lane

`R CMD check --as-cran` was running when the lane's session hit an
account limit at about 17:03 on 2026-09-22. Two minutes later 223 of
the 410 directories in the user library were empty. The lane library
and the round's reference build were untouched, and so was `pinlib`.
The evidence is in `dev/correct-log/loss7-evidence.txt` and the reading
is recorded in `dev/machine-library.md` ("The seventh loss"). The
restore recovered 223 of 223 and was verified by fitting.

What a reader of the first pass's numbers needed: **the suite and the
gated tier ran BEFORE the loss**, on RTMB 2.0 as the user installed it,
while its `R CMD check` ran after the restore, on RTMB 1.9, because the
restore takes CRAN's Windows binary. The handoff records that the two
agree on the objective to 2 ulp, so both toolchains have run these
tests.

**This no longer splits the numbers.** The user put RTMB 2.0 back
before punch round 1's second pass, and every tier in the generated
block below was rerun on RTMB 2.0, the check included.

## Runs

The tiers below were run on a lane library holding the punch-round
build (`C:/Users/adf44/source/r/correct2-lib`; the first pass used
`correct-lib`, and `CORRECT_LIB` selects either). `dev/correct-libmatch.R`
checks the library against the worktree sources function by function:
frmtmb 1118 identical, 0 differ; frmtmb.sample 209 identical, 0 differ.
Pointed at the reference build (`CORRECT_LIBMATCH_LIB=.../rellib-r3`)
the same script reports what the lane changed and nothing else, which is
what says it can see a stale install: frmtmb 1099 identical, 11 differ
(`brms_coef_table`, `brms_fixef_rows`, `pp_check.frmtmb_fit`,
`assemble_frame`, `mixture`, `parse_one_response`,
`residuals.frmtmb_fit`, `fill_prior_table`, `prior_table`,
`resolve_priorlist`, and the compat rules table), 8 new;
frmtmb.sample 201 identical, 6 differ (`pp_check.frmtmb_draws`,
`predictive_error.frmtmb_draws`, `default_prior_scale`,
`natural_dpar_prior`, `default_priors_for`, `default_prior_notes`), 2
new. `assemble_frame` and `pp_check.frmtmb_draws` joined that list in
punch round 1, at minors 4 and 2.

<!-- BEGIN GENERATED: dev/correct-runs.R -->
```
== generated by dev/correct-runs.R ==
suite    files 274  pass 16212  fail 0  error 0  skip 163  SUITE ran 274 of 274
gated    files 38  pass 3157  fail 0  error 0  skip 0  GATED ran 38 of 38
check    frmtmb-wt-correct: Status: 1 NOTE | frmtmb.sample: Status: OK
ledger0  | **total** | **494** | **243** | **71** | **35** | **145** |
ledger1  | **total** | **494** | **250** | **64** | **35** | **145** |
bitwise  identical: 121 of 126   differing: 5 
```
<!-- END GENERATED -->

The check's own test run is the whole suite of each package in ONE
process. **A passing "checking tests ... OK" does not echo the testthat
summary into the check log**, so those counts are not in `check.log` and
never were; punch round 1 quoted them and cited that file, which was
wrong about what the file holds. They live in
`<pkg>.Rcheck/tests/testthat.Rout`, and the low-disk rule deletes that
tree as soon as the Status line is read, so the two are in direct
tension. `dev/correct-run-check.ps1` now copies that file, a few KB, to
`dev/correct-log/<pkg>-testthat.Rout` before the tree can go. It DELETES
that destination first: a `-Force` copy onto a fixed name would leave
the previous run's file beside a log line saying this run produced none,
and a reader who quotes the file rather than the line would get the
wrong run's counts, which is the same failure again. The clearing is
verified in both directions rather than assumed: with a previous run's
file in place and no source to copy, the destination is GONE and the
log carries "NO testthat.Rout"; with a source, the destination holds
THIS run's content.

frmtmb FAIL 0, PASS 10700, SKIP 156
(`dev/correct-log/frmtmb-testthat.Rout`); frmtmb.sample FAIL 0, PASS
1788, SKIP 11 (`frmtmb.sample-testthat.Rout`). The check log itself is
`check-run4.log`, and what it holds is the phase lines and the two
Status lines, which is all it ever held.

frmtmb's ONE NOTE is this machine's environmental "V8 unavailable" on
the HTML manual, which every release records. The examples-timing NOTE
that two earlier checks of this same tree raised is not in this one; it
is load, measured above.


### Test files this lane added or changed

| file | why |
|---|---|
| `tests/testthat/test-ordinal-fitted.R` | the ordinal residual block is now the refusal; a multinomial block beside it, one per type since the recheck; punch round 1 pins the type the refusal quotes |
| `tests/testthat/test-pp-check-types.R` | new: every `available_ppc()` type against brms's verdict, on three responses; punch round 1 added a block for `resp` and an unused `group` or `x`, and the recheck split the refusals and the grouped types one behavior per block |
| `tests/testthat/test-prior-sd-scope.R` | new: which standard deviations a class `"sd"` prior reaches |
| `tests/testthat/test-mi-bsp-names.R` | new: `bsp_` and brms's `fixef()` order; punch round 1 pins the main effect before its interaction; the recheck split `hypothesis()` and the prior row out |
| `tests/testthat/test-mixture-reference.R` | new: the reference component follows the formula |
| `tests/testthat/test-smooths.R` | the REML-against-mgcv measurement |
| `tests/testthat/test-brms-methods.R` | gated: brms and frmtmb refuse residuals on the same shapes |
| `tests/testthat/test-brms-likelihood.R` | gated: row 17b, brms's density at frmtmb's `theta2 ~ x` estimate |
| `tests/testthat/test-prior-compat.R` | a block that asserted the OLD `sd` scope (`group = "g"` reaching two nonlinear parameters) now asserts brms's refusal |
| `tests/testthat/helper-brms-methods.R` | the exclusion's reason: both packages refuse, so there is no value to compare |
| `extensions/frmtmb.sample/tests/testthat/test-default-priors-brms.R` | new: the multivariate rows, the offset location, the per-prefix `sd` default, split at the recheck from what a user `sd` prior reaches |
| `extensions/frmtmb.sample/tests/testthat/test-draws-methods.R` | a category response refuses a predictive error; the multinomial guard meets that refusal first |
| `extensions/frmtmb.sample/tests/testthat/test-sampling-ported.R` | the loss model's `sd` default now carries its nonlinear parameter |
| the generated `test-brms-suite-*.R` | regenerated from the rebuilt ledger |
