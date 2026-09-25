# Lane mv: three multivariate gaps (2026-09-25)

Base: frmtmb 0.63.0 (commit 91cbf355). Library: `/opt/rlib/lane-mv`
(frmtmb and frmtmb.sample). Evidence files are in `dev/mv-log/`.

## What changed

### 1. Any number of `bf()` formulas summed with `+`

Cause: `+.frmtmb_formula` and `+.frmtmb_mvformula` were two different
S3 methods. When the operands of `+` have different methods, R warns
"Incompatible methods" and falls back to the internal `+`, so the
third `bf()` in a sum stopped with "non-numeric argument to binary
operator".

Fix (`R/bf.R`): both classes now carry the parent class
`frmtmb_bform`, and one method, `+.frmtmb_bform`, serves both. It
dispatches on the left operand to `plus_bf()` or `plus_mvbf()`. This
works on every R version the package supports (R >= 4.1), without
`chooseOpsMethod()`. A non-formula on the left (`gaussian() + bf()`)
is refused by name.

Also covered, each with a test:

- `bf() + bf() + bf() + bf()`, and `bf() + (bf() + bf())`, which
  flattens. brms refuses the parenthesized form ("Expecting a single
  value when fixing parameter 'forms'"); frmtmb accepts it.
- `set_rescor()` anywhere in the chain.
- A family after each `bf()`, and a family at the end.
- `lf()` and `nlf()` take brms's `resp =`. Without it, an `lf()` added
  to a multivariate formula is refused with the spelling that works
  (`+ lf(sigma ~ z, resp = "y3")`) and the list of responses. An
  unknown `resp`, and a `resp` that does not match the single `bf()`
  the `lf()` is added to, are refused by name. Before, the only
  advice was "put it directly after the bf()", which cannot reach the
  second of three responses.

`summary()` and `print()` of a multivariate fit printed `Family:` with
nothing after it and only the first formula. They now print brms's
`Family: MV(gaussian, gaussian, gaussian)` and one formula line per
response (`R/methods-fit.R`). A three-response gaussian rescor fit
reports a 3 x 3 `rescor_matrix()` and brms's names `rescor__y1__y2`,
`rescor__y1__y3`, `rescor__y2__y3`.

A DIVERGENCE KEPT, not introduced: a family added to a multivariate
formula fills only the responses without a family. brms's
`plus_mvbrmsformula()` applies `+ family` to every form, which
REPLACES each family, so in brms `bf(o ~ x) + cumulative() + bf(y1 ~
x) + gaussian()` makes `o` gaussian (measured on brms 2.23.0,
`/tmp/lanes/mv/mv-brms-plus.R`). frmtmb's rule predates this lane and
is what the gap's own example expects. It is now documented in
`?"+.frmtmb_bform"` (topic `plus-bform`) and in the brms-migration
vignette.

### 2. `student()` with `rescor = TRUE`

brms fits `multi_student_t(nu, Mu, Sigma)` with ONE `nu` and
`Sigma = D C D`, D the per-response sigmas, per row when sigma is
distributional (`make_stancode()`, read by hand:
`/tmp/lanes/mv/mv-brms-student.txt`). frmtmb now fits the same
density.

- Parse (`R/parse.R`, `rescor_share_nu()`): all responses must be
  student. The first response carries the one `nu` as its dpar,
  flagged `shared`; the other responses carry none.
- Objective (`R/objective.R`): the standardized residual matrix goes
  through `mvt_std_loglik()`, which is the multivariate-t block the
  autocor student case already used, factored out of
  `autocor_loglik()` in `R/autocor.R` so both call one function. The
  log-sigma Jacobian is the gaussian rescor one.
- Naming: `nu` has no response in it, as in brms: `variables()` lists
  `nu`, the internal coefficient is `nu_(Intercept)`, `summary()`
  prints `nu`, and `get_prior()` lists class `nu` with `resp` empty.
- Priors: `set_prior(class = "nu")` reaches the one `nu`. With `resp`
  it is refused, as `class = "rescor"` is ("takes no resp").
- Post-fit: every path that evaluates a response's dpars fills the
  shared `nu` into the responses that do not carry it
  (`with_shared_nu()` in `eval_dpars()`, `dpars_natural()`,
  `predict_dpar_values()`), because each marginal is a Student-t with
  that `nu`. `predict()` draws the responses jointly as a
  multivariate t: one normal draw through the correlation factor,
  divided by `sqrt(W / nu)` with one `W ~ chi^2(nu)` per row shared by
  the responses. `fitted()` is the n x 4 x K array. `simulate()` and
  `residuals()` refuse, exactly as for gaussian rescor. Nothing was
  widened.
- `log_lik()` (frmtmb.sample): its rescor branch computed the
  GAUSSIAN joint density unconditionally. It now calls a new core
  export `rescor_row_loglik()` (`R/predict-brms.R`, listed in
  `?frmtmb-sampling-api`), a plain-R implementation of brms's
  `multi_student_t_lpdf` term by term (and the gaussian one). Without
  this change a student rescor model's pointwise log-likelihood would
  have been silently gaussian: 38.63 log units off in total on the
  test fixture (`dev/mv-log/sample-pin.txt`).

Decisions:

- A predicted `nu` (`nu ~ z`) or a fixed one (`nu = 5`) is refused
  with brms's words, "Cannot predict or fix 'nu' in this model",
  measured on brms (`make_stancode()` refuses both). A multivariate t
  has one shape parameter per row, so a row-varying `nu` has no
  counterpart.
- A distributional sigma works (brms accepts it: it builds `Sigma[n]`
  per row). Validated below with `sigma ~ z` on one of three responses.
- A mix of gaussian and student responses is refused. brms refuses it
  too ("estimating 'rescor' is only possible in multivariate gaussian
  or student models", measured with `mvbf(bf(y1 ~ x, family =
  student()), bf(y2 ~ x, family = gaussian()), rescor = TRUE)`).

### 3. Ordinal families inside a multivariate model

The refusal lived in frame assembly (`R/frame.R`): every family that
declares `extra_pars` was refused in a multivariate model, because
the extras went into ONE global parameter list under the family's own
names (`tau_raw`), so two responses would collide and every consumer
read `fit$estimates$tau_raw` without knowing the response. "Extra
parameters" are whatever a family's `extra_pars` returns: the ordinal
thresholds (`cumulative`, `sratio`, `cratio`, `acat`: `tau_raw`), the
Cox baseline (`cox`: `sbhaz_raw`), and the class covariances of
`mixture_mvn()`; `mixture()` already refuses components with extras,
and an extension family may declare its own.

Fix: in a multivariate frame each response's family extras are stored
under a name that carries the response (`o_tau_raw`), with a map in
`frame$extra_map`. `fit_extras(fit, resp)` (new `resp` argument) and
the internal `resp_extras()` give a response its own block back under
its family's names, which is what its density and simulator read. A
univariate frame has no map and every path is unchanged (bitwise
check below). Response-aware now: the objective, `ordinal_ncat()`,
`ord_probs()`, `ord_prob_se()` (it perturbs each template entry once,
so a block is not counted twice), `sim_context()`, conditional
effects, `predict()`'s category proportions, the brms threshold rows
(`brms_extra_fixef()`, `hyp_put_ordinal()`, `ord_extra_comps()`),
`coef()`'s threshold shift, and the threshold prior
(`ordinal_threshold_entry()`, `has_ordinal_thresholds()`, the
`default_prior()` row). The `mo()` and `cs()` extras already had
unique global names and are untouched.

Lifted for `type = "ordinal"` families only. Every other family with
extras keeps a named refusal that now says why and that the ordinal
families are supported (`check_mv_extra_family()` in `R/predict.R`):
the Cox baseline and the mvn-mixture covariances are read by post-fit
code that has not been checked with a response in hand
(`cox_baseline()`, `mixture_posterior()`), and an extension family has
declared nothing about its extras.

Naming, as brms: `b_o_Intercept[1]`, `fixef()` rows `o_Intercept[1]`,
`vcov(full = TRUE)` rows `o_tau_raw_1`, `hypothesis(fit,
"o_Intercept[2] > o_Intercept[1]")`, and `set_prior(class =
"Intercept", resp = "o")` for the thresholds (a class `Intercept`
prior without `resp` on a multivariate model was already refused by
`resp_missing_refusal()`, brms's rule). `fitted()` of all responses
stacks the ordinal response as its `P(Y = k)` layers, the shape brms
2.23.0 gives (`fitted_mv()`, unchanged); `predict()` of an ordinal
response is asked one response at a time, as before for categorical
responses. frmtmb.sample's `log_lik()` now gives each response its own
block (`dev/mv-log/sample-pin.txt`: with the lane core and the old
frmtmb.sample the cumulative density received no thresholds and the
row-lpdf backstop stopped).

## Validation

All numbers from the scripts named, seeds inside them, on this
machine; outputs in `dev/mv-log/`.

### Student rescor, `dev/mv-validate-student.R`

Three responses, n = 400, `sigma ~ z` on the second, generated as a
multivariate t with 5 degrees of freedom.

| check | frmtmb | reference | difference |
|---|---|---|---|
| objective at a perturbed point vs `mvtnorm::dmvt` summed | -1808.490955043689 | -1808.490955043688 | relative 5.0e-16 |
| `logLik()` vs `sum(rescor_row_loglik())` | -1786.502695046127 | -1786.502695046127 | 0 at 16 digits |
| ML logLik vs hand-written RTMB objective, optimized separately | -1786.5026950461 | -1786.5026950653 | frmtmb higher by 1.9e-8 |
| nu at the optimum | 4.817534 | 4.817581 | |
| rescor (3 pairs) | 0.450770 0.181481 -0.350295 | 0.450767 0.181482 -0.350296 | |

Against brms 2.23.0 itself (`dev/mv-brms-student-loglik.R`, output
`dev/mv-log/brms-student-loglik.txt`): brms run with `algorithm =
"fixed_param"` from an init equal to a perturbed frmtmb parameter
point (two responses, n = 150), then `brms::log_lik()` summed, which
is brms's own R implementation of the multivariate t: -444.134618990114
against frmtmb's objective -444.134618990114, relative difference
2.6e-16.

The gaussian limit, as a check on the head and tail terms: at the
gaussian fit's parameters, log(nu - 1) = 5, 10, 20 gives -1873.10,
-1886.1795, -1886.2859416 against the gaussian -1886.2859465.

### Ordinal in a multivariate model, `dev/mv-validate-ordinal.R`

n = 600 in 60 groups; `o` cumulative (4 categories), `o2` sratio (3),
`y1` gaussian.

- IDENTITY, not a measurement: with no shared random effect the joint
  density factorizes and every parameter belongs to one factor, so the
  multivariate objective is the sum of the univariate ones. At a
  shared perturbed parameter point: mv -2360.95232649737, sum
  -2360.95232649737, residual 0.000e+00. At the three separate ML
  optima the sum differs from the mv optimum by 3.1e-8 (-2343.29696111157
  against -2343.29696108029), which is optimizer tolerance: the mv
  fit's `nlminb` stopped at max gradient 0.0029 (see Open questions).
- Against `ordinal::clm` (univariate cumulative logit): logLik
  -765.4534156202 (clm) against -765.4534156227 (frmtmb); thresholds
  -1.049105 0.041768 1.073592 (clm) against -1.049106 0.041767
  1.073588 from the MULTIVARIATE fit, slope 0.861733 against 0.861732.
- Shared `(1 | p | g)` across `o` and `y1`, against a hand-written
  RTMB objective (cumulative logit written from indicator sums,
  bivariate normal group effects, Laplace over the same effects):
  at a perturbed shared point -1610.540025385376 in both (relative
  0.000e+00); at the ML optima -1607.6213753016 (frmtmb) against
  -1607.6213753241 (reference); group-effect correlation 0.713075
  against 0.713069 (true 0.6).

- Against brms 2.23.0 (`dev/mv-brms-ordinal-loglik.R`, output
  `dev/mv-log/brms-ordinal-loglik.txt`), the same fixed-point
  construction on `bf(o ~ x, family = cumulative()) + bf(y1 ~ x)`:
  `brms::log_lik()` summed -400.061086992231 against frmtmb
  -400.061086992230 (relative 1.4e-16), and brms's `b_o_Intercept`
  draw equals frmtmb's thresholds (-1.032468, -0.2007506, 1.21134),
  which checks the threshold map and the centering together.

### Existing models are unchanged bit for bit, `dev/mv-bitwise.R`

`fn` and `gr` in hexadecimal at a fixed perturbed point, base against
lane, for univariate cumulative with `(1 | g)`, acat, a two-response
gaussian model with `(1 | g)`, gaussian rescor, and student with
`ar(cov = TRUE)` (the refactored multivariate-t block): `diff` is
empty (`dev/mv-log/bitwise-base.txt`, `bitwise-lane.txt`). No work was
added to any existing objective, so no timing was taken; the new
branches run only for the new models.

## Tests

New `tests/testthat/test-mv-gaps.R`, 8 blocks, 46 expectations
(agreement with `mvtnorm::dmvt`, the sum identity, the hand-written
|ID| objective, fitted() against the logit probabilities written out
from `fixef()`, every new refusal, naming). Against the base build all
8 blocks fail on behavior, not on missing symbols
(`dev/mv-log/base-pin.txt`): "non-numeric argument to binary
operator", "rescor = TRUE requires all responses to be gaussian (got:
student)", "Families with extra parameters ('cumulative') are not
supported in multivariate fits yet", and the old lf() message.

New block in `extensions/frmtmb.sample/tests/testthat/test-loo.R`
("the pointwise density of student rescor and mv ordinal fits"),
seen failing with the lane core and the old frmtmb.sample
(`dev/mv-log/sample-pin.txt`).

Tolerances are multiples of `.Machine$double.eps` times the magnitude
the run measured; none is absolute.

## Not done, and why

- `simulate()` and `residuals()` stay refused for every multivariate
  fit, rescor or not, as the brief asks.
- Cox, `mixture_mvn()` and extension families with extras stay
  refused inside a multivariate model (named refusal, compat rule
  `mvbf x cox`).
- brms's `mvbf() + family` replacement semantics were not adopted
  (above).
- frmtmb.sample's sampling-route defaults leave `nu` and ordinal
  thresholds flat, where brms uses `gamma(2, 0.1)` and
  `student_t(3, 0, 2.5)`. That is the package's documented existing
  behavior for univariate fits too (`R/sample.R`, "WHAT IS NOT
  MATCHED"), not something this lane introduced.
- `frm_linpred(resp = "y2", dpar = "nu")` on a student rescor model is
  refused; the message now names the response that carries the shared
  `nu`.
- The `thres()` and `0 + Intercept` lanes: nothing implemented. The
  namespacing is generic over whatever names `extra_pars` returns, so
  a `thres(gr = )` that adds threshold blocks is namespaced the same
  way; the threshold consumers look up `tau_raw` through
  `extra_tpl_name()`, which is the one place a new name would need
  adding.

## Open questions

- Multivariate models with ordinal responses often end with the
  "Large maximum absolute gradient" warning (0.0029 and 0.0017 on the
  validation fits, grad_tol 1e-3). The univariate cumulative fit of
  the same data already stops at 8.4e-4, and a joint fit sums the
  responses' gradients, so this looks like `nlminb`'s relative
  tolerance against a larger objective rather than a defect; the
  shared-point identity is exact. Not changed here.
- The generated brms-suite ports: no verdict row refers to these three
  gaps (`dev/brmsport-verdicts.tsv` searched; the brms tests that use
  a student rescor model are `stancode()` string checks, which have no
  port), so no port flips and the generator was not run. The gated
  tier (`FRMTMB_BRMS_FIT_TESTS=true`) was not run: no existing model's
  objective changed (bitwise check above), and the two brms
  comparisons above cover the new models directly.

## Files

R: `bf.R`, `parse.R`, `frame.R`, `objective.R`, `autocor.R`,
`predict.R`, `predict-brms.R`, `families.R`, `brms-names.R`,
`brms-shapes.R`, `confint.R`, `priors.R`, `methods-fit.R`,
`conditional-effects.R`, `insight.R`, `structure.R`, `sampling-api.R`,
`compat.R`. Generated: `NAMESPACE`, `man/lf.Rd`, `man/nlf.Rd`,
`man/plus-bform.Rd` (new), `man/frmtmb-extension-api.Rd`,
`man/frmtmb-sampling-api.Rd`. Docs: `_pkgdown.yml` (new topic),
`vignettes/brms-migration.Rmd`, `vignettes/case-studies.Rmd`,
`SPEC.md`. Extension: `extensions/frmtmb.sample/R/loo.R`, its
`tests/testthat/test-loo.R`. Tests: `tests/testthat/test-mv-gaps.R`.
Dev: `dev/mv-validate-student.R`, `dev/mv-validate-ordinal.R`,
`dev/mv-brms-student-loglik.R`, `dev/mv-brms-ordinal-loglik.R`,
`dev/mv-bitwise.R`, `dev/mv-log/`, this file.

frmtmb.sample now calls `rescor_row_loglik()` and `fit_extras(fit,
resp)`, which exist from this change on; its `DESCRIPTION` floor
(`frmtmb (>= 0.63.0)`) should move with the release that carries it.
DESCRIPTION was not edited, per the round rules.

## NEWS entry

- Any number of `bf()` formulas can be summed with `+`:
  `bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x)` used to stop with
  "non-numeric argument to binary operator", because the formula and
  the multivariate formula had two different `+` methods. `lf()` and
  `nlf()` take brms's `resp =` to name the response they belong to in
  a multivariate formula. `summary()` of a multivariate fit prints
  brms's `Family: MV(...)` line and every response's formula.
- `student()` with `rescor = TRUE` fits brms's multivariate Student-t:
  one `nu` shared by all responses (named `nu`, with no response),
  per-response sigma, distributional sigma allowed. A formula or a
  constant for `nu`, and a mix of gaussian and student responses, are
  refused as brms refuses them. Verified against `mvtnorm::dmvt` and
  against brms's own `log_lik()`.
  `predict()` draws the responses jointly from the multivariate t.
- Ordinal responses (`cumulative()`, `sratio()`, `cratio()`, `acat()`)
  can be part of a multivariate model, with their own thresholds
  (`b_o_Intercept[1]`, `set_prior(class = "Intercept", resp = "o")`)
  and `|ID|` group effects shared with other responses. Other families
  with extra parameters (`cox()`, `mixture_mvn()`) stay refused there,
  now with the reason.
- `fit_extras()` takes `resp`, and the new sampling-API export
  `rescor_row_loglik()` is the joint row density of a rescor fit.
- frmtmb.sample: `log_lik()` on a student rescor model used the
  gaussian joint density; it now uses the multivariate t. On a
  multivariate model with an ordinal response it gives each response
  its own thresholds.
