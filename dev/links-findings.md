# wt-links: brms link parity in the AD-safe link registry

Worktree `C:/Users/adf44/source/r/frmtmb-wt-links`, branch `wt-links`.
Private library `<scratchpad>/lk-lib`, Stan cache `<scratchpad>/lk-stan-cache`.
brms 2.23.0, rstan 2.32.7, StanHeaders 2.39.1, RTMB 1.9, R 4.6.1.

## Scope

Core's registry held identity, log, logit, cloglog, inverse, logm1,
tan_half and power12. The brief started at probit, probit_approx and
cauchit; the coordinator extended it to every brms 2.23.0 link name core
lacked. Measured from brms's own `.family_*` tables that is nine names:

    probit  probit_approx  cauchit  softit
    softplus  squareplus  sqrt  log1p  1/mu^2

All nine landed. After this branch, `setdiff(brms link names, core link
names)` is **empty**.

`1/mu^2` is NOT core's `power12`. `power12` is the tweedie power link, a
logit onto (1, 2) (`R/links.R:168`); `1/mu^2` is the inverse Gaussian
canonical link. No alias was possible; it is its own entry
(`R/links.R:227`).

## Two corrections to the brief

1. **softit's log-odds.** The coordinator gave it as
   `log(softplus(eta)) - log1p(softplus(eta))`. That is `log(mu)`, not
   the log-odds. With `y = softplus(eta)` and `mu = y / (1 + y)`,
   `1 - mu = 1 / (1 + y)`, so `mu / (1 - mu) = y` exactly and the
   log-odds is `log(y)`. Verified to 1e-15 at eta in {-3, -0.5, 0, 1, 4}
   (`dev/lk-measure3.R`). The registry uses `log(y)`.

2. **What `log_eta` is for.** The brief framed the robust fields purely
   as accuracy repairs, and by that test softplus and sqrt would have
   had none. Measurement said otherwise: the field is also the SWITCH
   that selects a family's robust log-density. See "The softplus
   finding".

## What landed

| file:line | change |
|---|---|
| `R/links.R:1-23` | registry header rewritten: the roster follows brms; the robust fields are present only where the plain path measurably saturates AND an exact form exists; each entry states its own measured boundary |
| `R/links.R:64` | `probit`, with `logit_eta` |
| `R/links.R:86` | `probit_approx`, with `logit_eta` |
| `R/links.R:107` | `cauchit`, no robust field, with the reason |
| `R/links.R:114` | `softit`, with `logit_eta` |
| `R/links.R:152` | `log1p` |
| `R/links.R:180` | `softplus`, with `log_eta` |
| `R/links.R:196` | `squareplus`, with `log_eta` |
| `R/links.R:210` | `sqrt`, with `log_eta` |
| `R/links.R:227` | `` `1/mu^2` ``, no robust field |
| `tests/testthat/test-numerical-robustness.R:365-` | 13 new `test_that()` blocks |
| `tests/testthat/test-brms-likelihood.R` | row 22, the link roster |
| `vignettes/brms-migration.Rmd:40-70` | the roster in the port table, plus a `### Links` subsection with the two porting cautions |
| `NEWS.md:31-62` | the development bullet |

`R/families.R` was not touched.

## brms agreement, per link

`linkinv` against `brms:::inv_link` and `linkfun` against `brms:::link`,
max relative difference over four linear predictors each
(`dev/lk-sanity.R`):

| link | linkinv | linkfun | round trip |
|---|---|---|---|
| probit | 0 | 0 | 5.6e-17 |
| probit_approx | 0 (vs `Phi_approx`) | 0 | not an inverse pair, in either package |
| cauchit | 0 | 4.4e-16 | 4.4e-16 |
| softit | 1.1e-16 | 2.2e-16 | 6.7e-16 |
| softplus | 1.2e-16 | 1.1e-16 | 2.2e-16 |
| squareplus | 0 | 0 | 4.4e-16 |
| sqrt | 0 | 0 | 0 |
| log1p | 0 | 0 | 5.6e-17 |
| 1/mu^2 | 0 | 0 | 0 |

**probit_approx is two different functions inside brms itself.** Its
Stan program uses `Phi_approx()`, which Stan defines as
`inv_logit(0.07056 x^3 + 1.5976 x)`; its R-side `brms:::inv_link()`
answers `pnorm()` for the same link name. Core reproduces the Stan form,
because that is the likelihood a brms fit was actually computed with.
`linkfun` is `qnorm`, brms's choice, so the pair is not an exact inverse
in either package. The test asserts this explicitly rather than working
around it.

## The saturation table

Each cell is the linear predictor at which that path first exceeds 1e-8
relative error against an independent reference (`stats::pnorm(log.p =)`,
`stats::pcauchy(log.p =)`, `stats::plogis(log.p =)`, or an asymptotic
series), or at which it stops being finite. Measured in
`dev/lk-measure2.R` and `dev/lk-measure3.R`.

| link | quantity | plain round trip | robust field | field |
|---|---|---|---|---|
| probit | `log(1 - mu)`, eta > 0 | **6.4** | 38.2 | `logit_eta` |
| probit | `log(mu)`, eta < 0 | 38.2 | 38.2 | same limit: pnorm underflow |
| probit_approx | `log(1 - mu)`, eta > 0 | **5.9** | never | `logit_eta` |
| probit_approx | `log(mu)`, eta < 0 | 21.4 | never | |
| cauchit | `log(1 - mu)`, eta > 0 | 3.2e9 | 1.8e10 | none |
| cauchit | `log(mu)`, eta < 0 | 1.8e10 | 1.8e10 | none |
| softit | `log(1 - mu)`, eta > 0 | 3.2e10 | never | `logit_eta` |
| softit | `log(mu)`, eta < 0 | -745 | -745 | same limit: softplus underflow |
| softplus | `log(mu)` | -745 | -745 | `log_eta`, for the branch |
| squareplus | `log(mu)`, eta < 0 | **-1e5** | never | `log_eta` |
| sqrt | `log(mu)` | \|eta\| = 1e-160 | same | `log_eta`, for the branch |
| 1/mu^2 | `log(mu)` | never | n/a | none |
| cloglog (existing) | `log(1 - mu)`, eta > 0 | 3.2 | 709 | `logit_eta` |

The bold cells are why each field exists. probit's plain round trip is
finished at a linear predictor of 6.4, an ordinary number for a
bernoulli fit; probit_approx's at 5.9; squareplus's at -1e5.

### Where each decision came from

- **probit**: `1 - pnorm(eta)` rounds to zero once `pnorm(eta)` rounds
  to one. `log(pnorm(eta)) - log(pnorm(-eta))` cannot cancel and holds
  until `pnorm` itself underflows near 38.2, which is also where the
  plain `log(mu)` gives out, so the other tail loses nothing. Included,
  with the 38.2 limit stated in the source comment.
- **probit_approx**: the inverse link is a logistic OF the cubic, so the
  log-odds IS the cubic, exact at any finite eta. Included.
- **cauchit**: NO exact branch-free form exists. The upper tail
  `atan(1 / eta) / pi` is exact only for eta > 0; below zero the same
  expression is off by pi, and choosing between them needs a branch. It
  is also unnecessary: the Cauchy tail is polynomial, `1 - mu` is
  `1 / (pi eta)` to leading order, so the plain round trip is still good
  at eta = 3.2e9. Field left absent, as the registry's convention
  allows, with both reasons in the comment.
- **softit**: `mu / (1 - mu)` is the softplus exactly. Included. The
  honest note is that it buys only the far upper tail (eta = 3.2e10),
  because softit's upper tail is polynomial too; both paths end together
  at eta = -745 where the softplus underflows.
- **squareplus**: `(eta + sqrt(eta^2 + 4)) / 2` IS `exp(asinh(eta / 2))`
  exactly, so `log_eta` is `asinh(eta / 2)`. This is the one place an
  identity buys five orders of magnitude: the plain sum catastrophically
  cancels below eta = -1e5 and is exactly zero by -1e9, where asinh
  answers -20.7.
- **1/mu^2**: `-0.5 * log(eta)` is exact, but `1 / sqrt(eta)` can
  neither underflow nor overflow at any representable positive eta, and
  `fam_inverse_gaussian` has no robust branch to select. Left absent.
- **log1p**: brms offers it on one dpar only, the `xi` of
  `gen_extreme_value`, whose support is (-1, Inf). No log or logit form
  applies to a dpar that may be negative. Neither field.

### The softplus finding

`softplus` and `sqrt` were going to have no `log_eta`: neither beats the
plain round trip at computing `log(mu)`. Sweeping `negbinomial()` with
each new positive-mean link at eta = -30 found otherwise. The AD
gradient there was -4.0007 (analytically -4); the central difference of
the same function was -2.33, stable under h = 3e-3, 1e-5 and 1e-6, so
the noise is in the FUNCTION, not the difference (`dev/lk-diag.R`).

The cause is `R/families.R:1075-1083`, which already says it: without a
`log_eta`, `robust_logmu()` returns NULL and negbinomial forms
`dnbinom2(y, mu, mu + mu^2 / shape)`. At mu = 9.4e-14 the excess
`mu^2 / shape = 4.4e-27` is within two significant digits of being lost
against mu, so the variance argument jitters as eta moves. With
`log_eta` present the family takes
`dnbinom_robust(y, lmu, 2 lmu - lsh)`, which never forms that sum.

So `log_eta` has a second job the header did not name: it is the switch
that selects the robust density. `log(softplus(eta))` is exact wherever
the softplus is (to eta = -745) and `log(eta^2)` is exact until eta^2
underflows at 1e-160, so both qualify on exactness, and both now carry
the field. The registry header and both entries say which of the two
reasons applies.

`log(eta^2)` rather than the algebraically equal `2 * log(abs(eta))`:
same number, no kink on the tape at eta = 0.

### The beta shape limit, which is not a link limit

`Beta(link = "probit_approx")` swept at eta = +/-30 gave a non-finite
log-density. The cause is not the link. `mu_pair()`
(`R/families.R:455-462`) has to exponentiate the log-odds back, because
the beta shapes are `mu * phi` and `(1 - mu) * phi`, and
`exp(log_inv_logit(lo))` underflows to exactly zero past `|lo| = 745`.
probit's log-odds reaches only 454 at eta = 30 and is fine;
probit_approx's log-odds is CUBIC in eta and is already at 1953 there.
The wall is at eta = 21.3.

bernoulli and binomial are unaffected at any eta, because
`dbinom_robust()` takes the log-odds directly and never exponentiates
it. The test sweeps beta with probit_approx at +/-20 and asserts the
underflow explicitly, so the limit is recorded rather than hidden.

## AD safety

Every primitive these links need tapes and differentiates under RTMB
(`dev/lk-probe5.R`): `RTMB::pnorm`, `RTMB::qnorm`, `RTMB::logspace_add`,
`atan`, `asinh`, `sqrt`, `expm1`, `log1p`, `tan`. Nothing branches.

`RTMB::pnorm`'s AD derivative is **bit-identical** to `stats::dnorm` at
eta in {-3, -0.5, 0, 1.2, 4}, and agrees with `numDeriv::grad` to 1e-11
relative. probit's `mu_eta` is therefore written as the normal density
in closed form, `exp(-0.5 eta^2) / sqrt(2 pi)`, and the test asserts the
two agree to 1e-15.

`pnorm` and `qnorm` must be written `RTMB::`-qualified: unqualified they
resolve to `stats::`, which has no advector method. Every other function
used here is a base generic RTMB overloads, matching the registry's
existing unqualified idiom.

A test tapes every field of every one of the 17 registry entries
(`linkfun`, `linkinv`, `mu_eta`, `logit_eta`, `log_eta`) and requires a
finite value and a finite gradient. That covers the eight pre-existing
links too.

## mu_eta, predict(se.fit) and conditional_effects

`mu_eta` is what `predict.frmtmb(se.fit = TRUE)` reads
(`R/predict.R:1381`), what `conditional_effects()` inverts for its band
(`R/conditional-effects.R:379`) and what the prior Jacobian uses
(`R/priors.R:1809`).

Against `numDeriv::grad(linkinv)`, max relative difference over each
link's own range: probit 1.6e-11, probit_approx 4.0e-11, cauchit
3.9e-11, softit 2.6e-11, softplus 4.6e-11, squareplus 2.0e-11, sqrt
1.1e-11, log1p 1.1e-11, `1/mu^2` 1.9e-11. That is numDeriv's own
accuracy, not the link's.

End to end on a fitted bernoulli, `predict(type = "response",
se.fit = TRUE)` against `|linkinv'(eta_hat)| * se_eta` with the
derivative taken by numDeriv: probit 3.7e-11, cauchit 2.9e-10,
probit_approx 9.2e-10, softit 4.1e-11. `conditional_effects()` returns
100 finite rows with the whole band strictly inside (0, 1) for all four.

## glm agreement

`stats::glm` and `frm()` both tightened, because at their defaults the
two stopping rules differ by about 2e-6 in the coefficients while
agreeing to 1e-9 in the log-likelihood, which measures the stopping
rules and not the links. frmtmb:
`optCtrl = list(rel.tol = 1e-14, x.tol = 1e-12)`; glm:
`glm.control(epsilon = 1e-14, maxit = 200)`.

| family and link | max abs coefficient difference | logLik difference |
|---|---|---|
| binomial(probit), `y \| trials(nt) ~ x + z` | 3.2e-9 | 2.3e-13 |
| binomial(cauchit), same | 5.0e-7 | 5.1e-11 |
| binomial(cloglog), same (regression check) | 1.0e-8 | 0 |
| binomial(logit), same (regression check) | 1.3e-8 | 2.3e-13 |
| poisson(sqrt), `y ~ x` | 2.6e-8 | 0 |
| inverse.gaussian(1/mu^2), `y ~ x` | 2.0e-8 | |

cauchit is the one that does not improve when frmtmb is tightened: 5.0e-7
with the default control and 5.0e-7 with the tight one. The residual is
glm's IRLS, not frmtmb's optimizer.

Two of these needed a design that keeps the linear predictor away from
zero. IRLS on a `sqrt` or `1/mu^2` link takes a square root of the
linear predictor and diverges the moment a step puts it negative, so the
poisson and inverse Gaussian responses are drawn THROUGH the link with a
bounded predictor. That is a property of `glm`, not of frmtmb, which
fits the unbounded design without complaint.

`probit_approx`, `softit`, `softplus` and `squareplus` have no `stats`
counterpart, so there is no glm row for them; they are covered by the
brms identity rows and by the delta-method and taping tests.

## The stats::make.link seam

`frm(bf(y ~ x), family = binomial(link = "softit"))` fails before frmtmb
sees anything:

    Error in make.link(link) : 'softit' link not recognised

`as_frmtmb_family()` (`R/families.R:357-390`) reads `x$link` off a
`stats::family` object and passes the string to core's own constructor,
which resolves it through `get_link()`. But `stats::binomial()`,
`stats::poisson()` and `stats::Gamma()` validate their link argument
through `stats::make.link()` at CONSTRUCTION, so the four names stats
does not know (`probit_approx`, `softit`, `softplus`, `squareplus`)
never reach that path.

They are reachable today through frmtmb's own family constructors:
`bernoulli()`, `Beta()`, `beta_binomial()`, `negbinomial()`,
`weibull()`, `geometric()`, `exponential()` and the rest of
`family_registry`. Verified fitting for all four. `probit`, `cauchit`,
`cloglog`, `sqrt` and `1/mu^2` are stats link names, so those work
through either route.

Closing this seam means core owning `binomial()`, `poisson()` and
`Gamma()` constructors that return `frmtmb_family` objects instead of
deferring to stats. That is `R/families.R` work, which belongs to a
sibling lane this round. Recorded, tested (the test asserts the error),
and documented in the migration vignette.

## The ordinal seam (recorded, not taken)

`R/families.R` was not touched, so the ordinal families keep their
hard-coded switch. The exact lines:

- `R/families.R:2352-2356`, `fam_cumulative()`:

      Fcdf <- switch(link,
        logit = function(x) 1 / (1 + exp(-x)),
        probit = function(x) RTMB::pnorm(x),
        stop("cumulative() supports links 'logit' and 'probit'", ...)
      )

- `R/families.R:2600-2604`, `ord_link_cdf()`, called at
  `R/families.R:2630` (sratio) and `R/families.R:2677` (cratio), with
  the same two-arm switch and the same error.
- `R/families.R:2721-2723`, `fam_acat()`, narrower still: `logit` only.

Both switches ARE `get_link(link)$linkinv` for the two links they
handle. Routing them through `get_link()` would give `cumulative()`,
`sratio()` and `cratio()` the `cloglog`, `cauchit`, `probit_approx` and
`softit` that brms has for them, and `acat()` all six, at the cost of a
line each. brms's tables say `cumulative` and `hurdle_cumulative` take
{logit, probit, probit_approx, cloglog, cauchit, softit}; `sratio` and
`cratio` take those minus softit; `acat` takes the full six.

One caveat the seam must carry, already noted at
`R/families.R:2357-2359`: cumulative's robust log-space CDF difference
exists for the logistic and not for the normal, because `RTMB::pnorm`
carries no `log.p`. Any link routed through `get_link()` there inherits
the non-robust branch unless the same exception is made for it.

## The link table against brms

`dev/lk-table.R`, computed from brms's own `.family_*` tables (custom
excluded, and `gamma`/`beta` mapped to core's `Gamma`/`Beta`
spellings). "core fam" counts the brms families that offer this link and
that core's `family_registry` also has.

| link | brms families | in core registry | core families | brms families core lacks |
|---|---|---|---|---|
| identity | 40 | yes | 28 | com_poisson dirichlet2 discrete_weibull frechet gen_extreme_value hurdle_negbinomial logistic_normal negbinomial2 wiener xbeta zero_inflated_beta_binomial zero_one_inflated_beta |
| log | 33 | yes | 23 | com_poisson dirichlet2 frechet gen_extreme_value hurdle_negbinomial negbinomial2 wiener xbeta zero_inflated_beta_binomial zero_one_inflated_beta |
| softplus | 25 | **new** | 18 | com_poisson dirichlet2 frechet gen_extreme_value hurdle_negbinomial negbinomial2 wiener |
| squareplus | 25 | **new** | 18 | the same seven |
| logit | 22 | yes | 15 | dirichlet dirichlet_multinomial discrete_weibull hurdle_cumulative xbeta zero_inflated_beta_binomial zero_one_inflated_beta |
| probit | 18 | **new** | 13 | discrete_weibull hurdle_cumulative xbeta zero_inflated_beta_binomial zero_one_inflated_beta |
| probit_approx | 18 | **new** | 13 | the same five |
| cauchit | 18 | **new** | 13 | the same five |
| cloglog | 18 | yes | 13 | the same five |
| softit | 16 | **new** | 11 | the same five |
| inverse | 16 | yes | 14 | frechet gen_extreme_value |
| sqrt | 9 | **new** | 6 | com_poisson hurdle_negbinomial negbinomial2 |
| tan_half | 1 | yes | 1 | |
| 1/mu^2 | 1 | **new** | 1 | |
| logm1 | 1 | yes | 0 | dirichlet2 |

Registry entries brms has no name for in a `.family_*` links vector:
`log1p` (brms has the name, but only on the `xi` dpar of
`gen_extreme_value`) and `power12` (frmtmb's own, for tweedie).

**brms link names still absent from core: none.** The remaining gap is
entirely a FAMILY gap, not a link gap: the 13 brms families core does
not have are `gen_extreme_value`, `xbeta`, `discrete_weibull`,
`hurdle_cumulative`, `com_poisson`, `negbinomial2`, `frechet`,
`wiener`, `dirichlet2`, `dirichlet_multinomial`, `logistic_normal`,
`zero_inflated_beta_binomial` and `zero_one_inflated_beta`. Each one
that arrives will find its links already registered.

### Which dpars the new links actually reach

The brief asked that every new link serve every dpar a (0, 1) parameter
appears in. Measured, it reaches `mu` and only `mu`, and that is a
`R/families.R` constraint rather than a registry one:

- `mu` of `bernoulli()`, `binomial()`, `Beta()`, `beta_binomial()`,
  `zero_inflated_binomial()` and `zero_inflated_beta()` takes the link
  the caller passes. All four new (0, 1) links were fitted through it.
- The zero-inflation and hurdle GATES are pinned. Every such family
  hard-codes `zi = "logit"` or `hu = "logit"` in its `links` list:
  `R/families.R:1391`, `:1430`, `:1478`, `:1746`, `:1787`, `:1835`,
  `:1887`, `:1974`. The gate link is not a user argument, so no link
  choice reaches it. brms does expose it, as `link_zi` / `link_hu`.
- Mixture weights are not a link at all. `theta1 ... theta{K-1}` are a
  multinomial logit against the last component
  (`R/families.R:2799`, `:3400`), so there is no per-dpar link slot to
  register anything into.

None of this is reachable without editing `R/families.R`, which this
lane did not touch. It is the same shape of seam as the ordinal switch
and the auxiliary-dpar note below.

Two dpar-level notes the table cannot show:

- brms lets `softplus` and `squareplus` be chosen for the auxiliary
  dpars `sigma`, `shape`, `phi`, `kappa`, `beta`, `disc`, `bs`, `ndt`
  and `alpha`. Core pins those to `log` inside each family
  (`links = list(mu = lk, shape = "log")`), so the choice is not
  user-reachable. That is an `R/families.R` seam, not a registry one.
- `log1p` is reachable only through `gen_extreme_value`'s `xi`. Core has
  no `gen_extreme_value` family, so the link is registered and correct
  but has no family to serve, and **no brms identity row is possible for
  it**. Recorded by name rather than tested.

## Tests added

`tests/testthat/test-numerical-robustness.R:365` onward. The file
held 17 blocks before this lane and holds 37 after, so the lane adds 20.
The first 15:

1. the new links reproduce brms's own definitions (gated on brms)
2. every registry field tapes with a finite value and gradient (all 17)
3. mu_eta is the derivative predict(se.fit) thinks it is (numDeriv), and
   RTMB::pnorm's AD derivative equals probit's mu_eta to 1e-15
4. probit's log-odds holds where the round trip has saturated
5. probit_approx's log-odds is the cubic and never saturates
6. softit's log-odds is the softplus's log
7. cauchit needs no robust field and says so by having none
8. squareplus's log mean is asinh, where the sum has cancelled
9. softplus and sqrt log_eta pick the robust density branch
10. the links with no exact robust form leave the field absent
11. the new (0, 1) links survive a separated predictor (bernoulli,
    binomial, beta at eta = +/-30, with the beta/probit_approx limit
    asserted)
12. the new positive-mean links survive an extreme predictor
    (negbinomial at eta = +/-30)
13. the new links fit a GLM and match stats::glm; predict(se.fit) is the
    delta method through every new link; the links stats::make.link
    rejects still reach frmtmb

## brms identity rows (row 22)

`tests/testthat/test-brms-likelihood.R`, in the tier's own numbering.
Row 21 is the family roster, so row 22 is the LINK roster: the plainest
possible model through each link name, `y ~ x`, checks A and B only, so
a divergence can only be the inverse link or its log-odds. The matrix in
`dev/brms-likelihood-tests.md` gains entry 22 and the results table
gains the rows below.

Residuals measured with `options(frmtmb.brms_lp_report = TRUE)`
(`dev/lk-brms-rows.R`). "constant" is `log_prob(brms) - logLik(frmtmb)`;
the tier admits only zero here unless brms carries a known flat
Dirichlet. "max abs gradient" is check B, tolerance 1e-3.

| row | result | constant | max abs gradient |
|---|---|---|---|
| 22a bernoulli(probit) | pass | 0 | 5.89e-06 |
| 22b bernoulli(probit_approx) | pass | 0 | 1.17e-04 |
| 22c bernoulli(softit) | brms cannot compile its own program (below) | | |
| 22d binomial(cauchit), trials(nt) | pass | 1.31e-12 | 6.83e-05 |
| 22e Beta(probit) | pass | -9.71e-12 | 1.55e-04 |
| 22f poisson(sqrt) | pass | -4.55e-13 | 2.30e-06 |
| 22g negbinomial(softplus) | pass | 5.68e-13 | 6.58e-04 |
| 22h negbinomial(squareplus) | pass | 3.41e-13 | 1.31e-04 |
| 22i inverse.gaussian(1/mu^2) | pass | -1.62e-13 | 2.06e-05 |

Every runnable row is an identity: the constant is zero to 1e-11 on all
eight, and the largest gradient is 6.58e-04 against the tier's 1e-3.
No new exemption is claimed for a numeric mismatch, because there is
none.

22b is the row that pins down the probit_approx reading. brms's Stan
program is `Phi_approx()` and its R `inv_link()` is `pnorm()`; frmtmb
computes the Stan one, and this row would fail against the other.

### 22c: brms 2.23.0 cannot compile its own softit program

Not a numeric divergence and not a skip. brms emits both softit helpers
with a vector-by-vector `/`, which Stan rejects:

    vector softit(vector p) {
      return log(expm1(-p / (p - 1)));
    }
    vector inv_softit(vector y) {
      return log1p_exp(y) / (1 + log1p_exp(y));
    }

    Semantic error in 'string', line 19, column 22 to column 34:
    Ill-typed arguments supplied to infix operator /.
    Instead supplied arguments of incompatible type: vector, vector.

Stan needs `./` for element-wise division of two vectors. `stanc` stops
on the FORWARD helper first, at the `softit(vector)` body on line 19;
`inv_softit` is wrong the same way one function later, so both are
asserted. NO brms model with this link compiles in 2.23.0, so there is
no program to check an identity against, for frmtmb or for anyone. frmtmb's softit is checked
against brms's R-side `inv_link()`/`link()` instead, where it agrees to
1.1e-16, and the row asserts the defect textually so that it fails the
day brms fixes the emitted code and a real identity can be written.

This follows the tier header's rule: an exemption is recorded by
asserting the structural difference rather than by skipping the row.

### Links with no possible identity row

- `log1p`: brms offers it only on `gen_extreme_value`'s `xi`, and core
  has no `gen_extreme_value` family. Untestable by name.
- `softplus` and `squareplus` on an auxiliary dpar (`sigma`, `shape`,
  `phi`, ...): brms allows the choice there, core pins those dpars to
  `log` inside each family. The rows above put both links on `mu` of
  `negbinomial()`, which is where core can reach them.

## Punch round (review `dev/reviews/2026-09-06-links.md`)

Verdict PUNCH, four items. All four addressed; nothing left standing.

### H1 (blocker): the conditional_effects() band on a decreasing link

Reproduced exactly as the review wrote it. `1/mu^2`: **99 of 100 rows
inverted, 1 row NaN**. `Gamma(inverse)`: **100 of 100 inverted** with
`R/links.R` irrelevant, which confirms the inversion is pre-existing and
that only the NaN is new with this link. `poisson(log)` control: 0.

Two defects meeting in one place, as the review said:

1. The band is built on the scale it is symmetric on and pushed through
   `linkinv` without reordering. `1/mu^2`'s `linkinv` is `1 / sqrt(eta)`,
   which DECREASES, so the ends arrive swapped.
2. `1 / sqrt(eta)` has no value at `eta <= 0`, so once the band reaches
   zero the transform gives NaN, not merely a misordered number.

Three sites transform band endpoints. Two already sorted with
`pmin`/`pmax` and one did not, and the unsorted one is the branch this
family takes:

| site | before | after |
|---|---|---|
| `R/conditional-effects.R:1922` (expected-response Wald band) | sorted with `pmin`/`pmax` | `ce_band_ends()` |
| `R/conditional-effects.R:1934` (link-scale Wald band) | **not ordered at all** | `ce_band_ends()` |
| `R/conditional-effects.R:1975` (profile band) | sorted with `pmin`/`pmax` | `ce_band_ends()` |

The `pmin`/`pmax` the first and third sites used is not sufficient
either: with one end NA it fills the missing bound from the other and
reports a band the data does not support.

**The fix** is one helper, `ce_band_ends()` at
`R/conditional-effects.R:414`, and it does three things:

- **Orders by direction, not by sorting.** A monotone decreasing link
  has the ends SWAPPED wholesale (`dm < 0` throughout), so a bound that
  is missing stays on the side it belongs to. For `1/mu^2` the
  unreachable end is the UPPER one, because the mean runs away as eta
  falls to zero. The measured row now reads `lower = 5.658,
  upper = NA`, where before it was `lower = NaN, upper = 5.658`.
- **Refuses an end outside the domain**, returning NA rather than NaN: a
  finite linear predictor whose response is not finite.
- **Refuses a pole BETWEEN two finite ends.** This is `inverse` with the
  band straddling `eta = 0`: both ends are finite, neither bounds the
  response, and a naive sort would report a plausible-looking interval
  that is wrong. Detected without naming any link, by checking that the
  midpoint lies between the ends.

One warning, accumulated across panels (`R/conditional-effects.R:1816`
for the counter, `:1991` for the warning), naming the link and the
count, in the idiom the existing `pfail` warning already uses:

    The '1/mu^2' link does not reach the band at 1 of 100 grid
    point(s): the interval runs past the link's domain, so that bound
    is NA rather than a number

Measured after the fix: `1/mu^2` **0 of 100 inverted**, one NA upper
bound, one warning. `Gamma(inverse)` **0 of 100 inverted**, 0 non-finite.
`poisson(log)` unchanged at 0.

`predict(se.fit = TRUE)` does NOT go through this path and has no
endpoints to order: it returns `|mu_eta| * se_eta`, non-negative by
construction. Asserted rather than assumed, in
"predict(se.fit) reports a standard error, not an interval".

The hunk is confined to `R/conditional-effects.R` and adds one function;
it does not touch `ce_pred_dpar()`, which main gained from the sample-ce
lane after this branch started.

Five new test blocks cover it: a band assertion for **all nine** new
links on a family that accepts each (the review correctly noted only the
four (0, 1) links had one), `Gamma(inverse)` and
`inverse.gaussian("1/mu^2")` ordered-band regressions, the helper's own
four cases, and the `predict(se.fit)` claim. `log1p` has no brms family
core carries, so its band is asserted on `negbinomial(link = "log1p")`,
where `expm1(eta)` is a valid positive mean.

### H1 follow-up: the swap had to be per row, not per call

Raised on the re-check, low severity. `ce_band_ends()` decided the swap
once for the whole vector with `all(dm[fin] < 0)`. That is right for
`inverse` and `1/mu^2`, which decrease everywhere, and wrong for `sqrt`,
whose `mu_eta` is `2 * eta` and changes sign at zero.

Measured with the old rule, `lo = c(-3, 1)`, `hi = c(-1, 3)`, so
`dm = c(-4, 4)`: `all(dm < 0)` is FALSE, nothing swaps, and the band
comes back `lower = c(9, 1)`, `upper = c(1, 9)` - **1 of 2 rows
inverted**. Per row it is `lower = c(1, 1)`, `upper = c(9, 9)`, 0
inverted.

No fitted model reaches this today: a `sqrt` fit's likelihood depends on
`eta^2` only, so it settles on the positive branch and `dm` has one sign.
It is a latent inversion in a helper, not a live defect, and it is fixed
rather than documented.

`R/conditional-effects.R:432` now reads

    dec <- rep_len(is.finite(dm) & dm < 0, length(a))

and swaps `a[dec]` against `b[dec]`. The `rep_len()` is the second half
of the same bug: a length-1 derivative against a longer band would have
indexed only the first row and left the rest unswapped. Both are
asserted, at `tests/testthat/test-numerical-robustness.R:835` (mixed
sign) and `:846` (length-1 derivative against a two-row band).

### M1: NEWS said four robust fields, six exist

`NEWS.md:40`. Six carry one (`probit`, `probit_approx`, `softit` with
`logit_eta`; `softplus`, `squareplus`, `sqrt` with `log_eta`); four of
the six move a saturation boundary. Rewritten to say exactly that.

### L1: counts

- new blocks in `test-numerical-robustness.R`: the report said 13, the
  review measured 15, and after the punch round's five it is **20**
  (17 before the lane, 37 after).
- gated brms tier: **394**, not 393. Corrected.
- `dev/brms-likelihood-tests.md` said "zero to 1e-11" in prose and
  "0 (to 1e-9)" in the table. Reconciled to the measurement: the largest
  constant over the eight identity rows is **9.71e-12** (row 22e), so
  the nine added rows now read "0 (to 1e-11)" and the prose names the
  9.71e-12.

### L2: the softit assertion was over-broad

`expect_false(grepl("./", code, fixed = TRUE))` forbade the two
characters anywhere in the program, including in an unrelated
elementwise division or a comment, and would have blamed softit. Now
scoped to the two helper bodies:
`log(expm1(-p ./ (p - 1)))` and `log1p_exp(y) ./ (1 + log1p_exp(y))`.
Row 22c goes from 4 assertions to 5; the tier total is unchanged at 394.

### L3: prose named the wrong helper

`stanc` stops on the FORWARD `softit(vector)` at line 19, not on
`inv_softit`. Both are the same defect and the test asserted both, so
only the prose was wrong. Corrected here and in the row 22c comment.

### Recorded, not actioned

M2 (`summary()` prints no link name for any family) and M4 (no
`link_sigma` / `link_shape`, so a non-`mu` dpar link is unreachable) are
follow-ups the coordinator assigned elsewhere. M3 (`probit_approx` is
not an inverse pair) is brms's own behavior, already documented at
`R/links.R:80-85` and in the vignette; no change. I agree with the
review on all three.

One number I read differently: the review counts row 22's positive-mean
test at 9 assertions; I measure 8, which is 4 `brms_lp_check()` calls at
2 assertions each. The tier total of 394 agrees either way.

## Verification

Counts audited by name, each in a fresh process against the build
installed in `lk-lib`.

| what | result |
|---|---|
| `test-numerical-robustness.R` | 673 pass, 0 fail, 0 error, 0 skip (37 blocks) |
| `test-brms-likelihood.R` (gated, one process) | 394 pass, 0 fail, 0 error, 0 skip |
| `test-message-uniqueness.R` | 6 pass, 0 fail, 0 error, 0 skip |
| `test-ce-bands.R` | 167 pass, 0 fail, 0 error, 0 skip |
| `test-ce-facets.R` | 92 pass, 0 fail, 0 error, 0 skip |
| `test-custom-family.R` | 41 pass, 0 fail, 0 error, 0 skip |
| `test-input-validation.R` | 43 pass, 0 fail, 0 error, 0 skip |
| `test-families.R` | 32 pass, 0 fail, 0 error, 0 skip |
| `test-famgaps.R` | 99 pass, 0 fail, 0 error, 0 skip |
| `test-ce-bands.R` | 167 pass, 0 fail, 0 error, 0 skip |
| `test-effects.R` | 54 pass, 0 fail, 0 error, 0 skip |
| `devtools::document()` | idempotent; no `man/` or `NAMESPACE` change |
| `R CMD check --as-cran` | **0 errors, 0 warnings, 1 NOTE** (vignettes built with pandoc 3.8.3, suite 16m). The NOTE is "Skipping checking math rendering: package 'V8' unavailable", an absent optional package in this environment and not attributable to the diff. |

### Two harness traps worth recording for the sibling lanes

1. **Running the tier file needs the package namespace.** Invoked as
   `Rscript -e 'library(testthat); library(frmtmb); test_file(...)'`
   the tier reports 9 errors that are not failures at all:
   `could not find function "ord_tau_from_raw"`,
   `object 'covstruct_registry' not found`,
   `could not find function "autocor_natural"`, `"car_rho"`, `"coef_b"`.
   `helper-brms.R` calls those internals unqualified, and a bare
   `test_file()` evaluates in the global environment where they are not
   visible. Under `R CMD check` testthat runs inside the namespace, so
   the file is green there. The invocation that reproduces it outside
   the check is

       test_file("test-brms-likelihood.R",
                 env = new.env(parent = asNamespace("frmtmb")))

   which gives 393 pass, 0 error. `test-numerical-robustness.R` does not
   have this problem because it spells every internal `frmtmb:::`.

2. **`R CMD build --no-build-vignettes` costs three CRAN warnings.**
   `checking files in 'vignettes'`, `checking package vignettes` and
   `Directory 'inst/doc' does not exist` are all artifacts of the flag,
   not of any change. A real `--as-cran` run needs the vignettes built,
   which needs pandoc: it is not on PATH in this environment, and
   RStudio's bundled copy is at

       C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools

   (pandoc 3.8.3), exported as `RSTUDIO_PANDOC`. Without it every
   vignette fails identically with "Pandoc is required to build R
   Markdown vignettes but not available".
