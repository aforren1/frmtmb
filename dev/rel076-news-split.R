# The sentences over 25 words in the 0.76.0 NEWS sections
# (dev/rel076-news-long.R, dev/rel076-log/news-long-1.txt: core 17,
# frmtmb.sample 14, frmtmb.spatial 2), split into shorter ones, as
# ASD-STE100 asks, and the edits of fact the merge needs.
# dev/rel076-news-merge.R applies each pair (lane text, new text) to the
# top section of the file named, before it rewraps the section; a run of
# whitespace in the lane text matches any whitespace, so a pair does not
# depend on where the lane wrapped its lines. The numbers and the facts
# are the lanes'; only the sentence boundaries and the order of clauses
# change, except where a comment says otherwise. Adapted from
# dev/rel075-news-split.R.
list(
  "NEWS.md" = list(
    c(r"---[Since `mo()` was added in 0.15.0, frmtmb coded the integers
      themselves, so the simplex had one weight for each value between 0
      and the minimum that no row has.]---",
      r"---[Since `mo()` was added in 0.15.0, frmtmb coded the integers
      themselves. So the simplex had one weight for each value between 0
      and the minimum that no row has.]---"),
    c(r"---[Now the fit on 1..4 is the fit on 0..3, with the same
      log-likelihood, coefficients and finite standard errors, and
      `frm_sample()` agrees with brms within 2 Monte Carlo standard errors
      on every parameter.]---",
      # with an edit of fact, a release fix (the release review, finding
      # 1): a two-valued integer is now refused, which no lane's bullet
      # said (dev/review-rel076-mo2.R)
      r"---[Now the fit on 1..4 is the fit on 0..3, with the same
      log-likelihood, coefficients and finite standard errors.
      `frm_sample()` agrees with brms within 2 Monte Carlo standard errors
      on every parameter. An integer `mo()` variable with only two values,
      such as 1 and 2, is now refused ("mo() needs at least 3 ordered
      categories"). 0.75.0 fitted it with an empty category 0, so its
      coefficients had no standard errors. brms fits it as a linear term
      in the code.]---"),
    c(r"---[A value that is not a whole number, or that is below the fitted
      minimum or above the fitted maximum, is refused, as brms's
      `validate_newdata()` refuses it.]---",
      r"---[A value that is not a whole number is refused, as brms's
      `validate_newdata()` refuses it. So is a value below the fitted
      minimum or above the fitted maximum.]---"),
    c(r"---[A missing value, a factor or a character vector in new data is
      refused with its own reason, and `mo(mi(x))` is refused at the fit,
      as brms refuses it.]---",
      r"---[A missing value, a factor or a character vector in new data is
      refused with its own reason. `mo(mi(x))` is refused at the fit, as
      brms refuses it.]---"),
    c(r"---[It held the variable at the mean and the mean +/- sd, which the
      rounding above then moved to another category, so a display was off
      by a whole category and its labels named values the model was not
      evaluated at.]---",
      r"---[It held the variable at the mean and the mean +/- sd, which the
      rounding above then moved to another category. So a display was off
      by a whole category, and its labels named values the model was not
      evaluated at.]---"),
    c(r"---[On 200 rows with the values 1, 2 and 5, the standard errors at
      1, 2 and 5 were 0.121, 0.0987 and 0.121, and are now 0.122, 0.114
      and 0.121; brms's posterior sds are 0.124, 0.116 and 0.123.]---",
      r"---[On 200 rows with the values 1, 2 and 5, the standard errors at
      1, 2 and 5 were 0.121, 0.0987 and 0.121. They are now 0.122, 0.114
      and 0.121, and brms's posterior sds are 0.124, 0.116 and
      0.123.]---"),
    c(r"---[now has a `NaN` standard error, as the fit's warning says;
      through 0.75.0 it was finite (0.0859 and 0.0950).]---",
      r"---[now has a `NaN` standard error. The fit's warning says so.
      Through 0.75.0 it was finite (0.0859 and 0.0950).]---"),
    c(r"---[and `get_vcov()`: on 300 rows of `y ~ mo(k) + z`,
      `predictions()` at k = 2 and 3 gave 0.0846 and 0.1035 where the exact
      standard errors are 0.1220 and 0.1142, and it gives those now.]---",
      r"---[and `get_vcov()`. On 300 rows of `y ~ mo(k) + z`,
      `predictions()` at k = 2 and 3 gave 0.0846 and 0.1035. The exact
      standard errors are 0.1220 and 0.1142, and it gives those now.]---"),
    c(r"---[A fit saved by frmtmb 0.75.0 or earlier with an integer `mo()`
      variable that does not start at 0 is refused, with a message to fit
      the model again, by every method that reports an estimate, a degree
      of freedom, a standard error or a comparison (`fixef()`, `coef()`,
      `logLik()`, `AIC()`, `BIC()`, `vcov()`, `confint()`, `summary()`,
      `anova()`, `drop1()`) and by `fitted()`, `predict()`, `simulate()`,
      `residuals()`, `conditional_effects()` and `frm_sample()`.]---",
      r"---[A fit saved by frmtmb 0.75.0 or earlier is refused when its
      integer `mo()` variable does not start at 0. The message says to fit
      the model again. Every method that reports an estimate, a degree of
      freedom, a standard error or a comparison refuses it: `fixef()`,
      `coef()`, `logLik()`, `AIC()`, `BIC()`, `vcov()`, `confint()`,
      `summary()`, `anova()` and `drop1()`. So do `fitted()`, `predict()`,
      `simulate()`, `residuals()`, `conditional_effects()` and
      `frm_sample()`.]---"),
    c(r"---[This is the rule that `sar()`'s `lagsar` got in 0.75.0, and it
      is now one rule for core's classes and for each bounded parameter of
      a term registered with `frmtmb_register_autocor()`.]---",
      r"---[This is the rule that `sar()`'s `lagsar` got in 0.75.0. It is
      now one rule for core's classes and for each bounded parameter of a
      term registered with `frmtmb_register_autocor()`.]---"),
    c(r"---[For `cosy` it lists the window `frm()` fits it on, `-1/(d - 1)`
      to 1 for `d` time points, and on `route = "sample"` brms's 0 to 1,
      where frmtmb.sample now samples it.]---",
      r"---[For `cosy` it lists the window `frm()` fits it on, `-1/(d - 1)`
      to 1 for `d` time points. On `route = "sample"` it lists brms's 0 to
      1, where frmtmb.sample now samples it.]---"),
    c(r"---[`sigma()` and `VarCorr()` report it; brms's `sigma` is the
      innovation sd,]---",
      r"---[`sigma()` and `VarCorr()` report it. brms's `sigma` is the
      innovation sd,]---"),
    c(r"---[was a density on the marginal sd; now it is brms's density at
      `sigma * c`,]---",
      r"---[was a density on the marginal sd. Now it is brms's density at
      `sigma * c`,]---"),
    c(r"---[A `lb` above 0 or a `ub` on `sigma` there is refused, because a
      box on the innovation sd moves with `ar` and is no box on the
      marginal sd frmtmb holds; brms 2.23.0 takes such a bound on its
      innovation sd (`real<lower=0.2> sigma`).]---",
      r"---[A `lb` above 0 or a `ub` on `sigma` there is refused. A box on
      the innovation sd moves with `ar`, so it is no box on the marginal
      sd frmtmb holds. brms 2.23.0 takes such a bound on its innovation sd
      (`real<lower=0.2> sigma`).]---"),
    c(r"---[the pointwise log-density of a response with a residual
      covariance term (the covariance form of `ar()`, `ma()` and `arma()`,
      `cosy()`, `unstr()`, `fcor()` and a registered term).]---",
      r"---[the pointwise log-density of a response with a residual
      covariance term. These terms are the covariance form of `ar()`,
      `ma()` and `arma()`, `cosy()`, `unstr()`, `fcor()` and a registered
      term.]---"),
    c(r"---[One draw costs one small inverse per pattern of time points:
      0.64 ms against 17.8 ms for brms's per-group dense inverses on 400
      groups of 10 rows, and 2.6 ms against 97.5 ms for `fcor()` on 600
      rows.]---",
      r"---[One draw costs one small inverse per pattern of time points. On
      400 groups of 10 rows a draw takes 0.64 ms, against 17.8 ms for
      brms's per-group dense inverses. For `fcor()` on 600 rows it takes
      2.6 ms, against 97.5 ms.]---"),
    c(r"---[A given `jacobian` is checked against the map when the frame is
      built, at the start and at a generic point (sigma not 1 and varying
      by row, every parameter moved off the start), and a wrong one is
      refused by the term's name.]---",
      r"---[A given `jacobian` is checked against the map when the frame is
      built. The check runs at the start and at a generic point: sigma not
      1 and varying by row, and every parameter moved off the start. A
      wrong one is refused by the term's name.]---")
  ),
  "extensions/frmtmb.sample/NEWS.md" = list(
    # an edit of fact: the lane could not write the release's number
    c(r"---[With the frmtmb that codes an integer `mo()` variable from its
      minimum, as brms 2.23.0 does, `frm_sample()` samples brms's model on
      such data.]---",
      r"---[With frmtmb 0.76.0, which codes an integer `mo()` variable from
      its minimum as brms 2.23.0 does, `frm_sample()` samples brms's model
      on such data.]---"),
    c(r"---[where they had 4, and every posterior mean agrees with brms's
      within 2 Monte Carlo standard errors,]---",
      r"---[where they had 4. Every posterior mean agrees with brms's
      within 2 Monte Carlo standard errors,]---"),
    # an edit of fact: the release's number
    c(r"---[its tests need that frmtmb.]---",
      r"---[its tests need frmtmb 0.76.0.]---"),
    c(r"---[internal parameter: for a first-order `cov = TRUE` `ar` or
      `ma` that is `(1 - rho^2)^(-3/2)` on `rho`, with its mass at -1 and 1,
      and for `cosy` a density flat on the logit of its whole window.]---",
      r"---[internal parameter. For a first-order `cov = TRUE` `ar` or
      `ma` that is `(1 - rho^2)^(-3/2)` on `rho`, with its mass at -1 and
      1. For `cosy` it is a density flat on the logit of its whole
      window.]---"),
    c(r"---[Against brms 2.23.0's own draws (4 chains of 4000 on each
      side, seed 3), the posterior means were 42.3 Monte Carlo standard
      errors apart for `ar[1]` (AR(1) of 0.9, 8 groups of 6 time points)
      and 19.6 for `sigma`; an MA(1) of 0.6 had every `ma[1]` draw at
      1.]---",
      r"---[The reference is brms 2.23.0's own draws, 4 chains of 4000 on
      each side (seed 3). Against them, the means of `ar[1]` were 42.3
      Monte Carlo standard errors apart, on an AR(1) of 0.9 in 8 groups of
      6 time points. For `sigma` they were 19.6 apart. An MA(1) of 0.6 had
      every `ma[1]` draw at 1.]---"),
    c(r"---[48 percent of the draws were below 0, where brms has none, and
      the mean was 0.011 against brms's 0.076; now there are none, and the
      means are]---",
      r"---[48 percent of the draws were below 0, where brms has none. The
      mean was 0.011 against brms's 0.076. Now there are none, and the
      means are]---"),
    c(r"---[is flat on its stationary or invertible region, and its `sigma`
      default stays on the marginal sd.]---",
      r"---[is flat on its stationary or invertible region. Its `sigma`
      default stays on the marginal sd.]---"),
    # a release fix (acprior's review, fix round 2 check, consolidation
    # note 2): the headline says "every prior on class `sigma`", since a
    # prior on the sigma intercept under `sigma ~ x` is not read at the
    # innovation sd (?set_prior, Residual correlation)
    c(r"---[is the marginal sd, and every prior on it is read at the
      innovation sd, as brms's is.**]---",
      r"---[is the marginal sd. Every prior on class `sigma` is read at the
      innovation sd, as brms's is.**]---"),
    c(r"---[report frmtmb's marginal sd (on the AR(1) above, a posterior
      mean of 1.16), where brms's `sigma` is the innovation sd (0.43); the
      model is the same.]---",
      r"---[report frmtmb's marginal sd: on the AR(1) above, a posterior
      mean of 1.16. brms's `sigma` is the innovation sd (0.43); the model
      is the same.]---"),
    c(r"---[is read where the default is, at the innovation sd, so a
      `default_prior(route = "sample")` table]---",
      r"---[is read where the default is, at the innovation sd. So a
      `default_prior(route = "sample")` table]---"),
    c(r"---[Before, a user prior there was read on the marginal sd, a
      density up to 0.72 nats from brms's]---",
      r"---[Before, a user prior there was read on the marginal sd. That
      density was up to 0.72 nats from brms's]---"),
    c(r"---[is flat on [0.5, 1] in brms, and was flat on frmtmb's internal
      parameter here (mean 0.946 against brms's 0.916, z 43.8); now 0.916
      (z 0.5).]---",
      r"---[is flat on [0.5, 1] in brms. Here it was flat on frmtmb's
      internal parameter (mean 0.946 against brms's 0.916, z 43.8); now
      the mean is 0.916 (z 0.5).]---"),
    c(r"---[the other rows of its group, normal or Student-t, and equals
      brms's `log_lik()`]---",
      r"---[the other rows of its group, normal or Student-t. It equals
      brms's `log_lik()`]---"),
    c(r"---[On a group with a missing time point the conditional is the
      fitted model's, whose `ar()`, `ma()` and `arma()` lags count time
      levels where brms's count positions.]---",
      r"---[On a group with a missing time point, the conditional is the
      fitted model's. Its `ar()`, `ma()` and `arma()` lags count time
      levels, where brms's count positions.]---"),
    c(r"---[gave 5 to 9 divergent transitions per 4 chains of 2000 on a
      `cosy` near its bound at 0 (seeds 3 to 5, mode and random starts),
      where brms has none; now there are none.]---",
      r"---[gave 5 to 9 divergent transitions per 4 chains of 2000. That
      was on a `cosy` near its bound at 0 (seeds 3 to 5, mode and random
      starts), where brms has none. Now there are none.]---"),
    c(r"---[the bound is the default, brms's declaration, and names]---",
      r"---[the bound is the default, brms's declaration. It names]---")
  ),
  "extensions/frmtmb.spatial/NEWS.md" = list(
    c(r"---[as a sparse matrix, so frmtmb.sample's `log_lik()` returns brms
      2.23.0's leave-one-out conditional density for both types, gaussian
      and student, equal to brms's `log_lik()`]---",
      r"---[as a sparse matrix. So frmtmb.sample's `log_lik()` returns brms
      2.23.0's leave-one-out conditional density for both types, gaussian
      and student. It equals brms's `log_lik()`]---"),
    c(r"---[as `lb` and `ub`, as brms's does, and `frm()` refused them as
      outside the interval.]---",
      r"---[as `lb` and `ub`, as brms's does. `frm()` refused them as
      outside the interval.]---")
  )
)
