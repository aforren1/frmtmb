# The sentences over 25 words in the 0.73.0 NEWS sections
# (dev/rel073-news-long.R, dev/rel073-log/news-long-1.txt: core 21,
# frmtmb.sample 7), split into shorter ones, as ASD-STE100 asks.
# dev/rel073-news-merge.R applies each pair (lane text, new text) to the
# top section of the file named, before it rewraps the section; a run
# of whitespace in the lane text matches any whitespace, so a pair does
# not depend on where the lane wrapped its lines. The numbers and the
# facts are the lanes'; only the sentence boundaries and the order of
# clauses change.
list(
  "NEWS.md" = list(
    c(r"---[200 poisson draws at `lambda = 4` recorded as "6 or more" and
      right censored at 6 move by 30.2 in log-likelihood and 5.5 percent
      in `lambda`; right censored at 5 they give the old fit exactly.]---",
      r"---[200 poisson draws at `lambda = 4` recorded as "6 or more" and
      right censored at 6 move by 30.2 in log-likelihood and 5.5 percent
      in `lambda`. Right censored at 5, they give the old fit exactly.]---"),
    c(r"---[names the family** and points at `frm_compat()` for the term
      that was used, the list of the families that have a CDF, instead of
      listing them.]---",
      r"---[names the family.** It points at `frm_compat()` for the term
      that was used, which lists the families that have a CDF. It no
      longer lists them itself.]---"),
    c(r"---[and `(1 | g)` in both `mu` and `sigma` used to fit, but the
      Gauss-Kronrod rule integrated only some of the effects and left the
      rest to the Laplace approximation: over `(1 | ga/gb)` the
      log-likelihood]---",
      r"---[and `(1 | g)` in both `mu` and `sigma` used to fit. But the
      Gauss-Kronrod rule integrated only some of the effects and left the
      rest to the Laplace approximation. Over `(1 | ga/gb)` the
      log-likelihood]---"),
    c(r"---[When every row with a positive weight has a random effect of its
      own in `mu` and the weights sum to more than the number of rows, the
      likelihood of a family whose dispersion formula has an intercept
      grows without bound as the dispersion goes to its limit (`sigma` to
      0, `shape` to infinity), like the peak height of the density to the
      power `sum(w) - n`.]---",
      r"---[The refusal applies when two conditions hold. First, every row
      with a positive weight has a random effect of its own in `mu`.
      Second, the weights sum to more than the number of rows. If the
      dispersion formula has an intercept, the likelihood grows without
      bound as the dispersion goes to its limit (`sigma` to 0, `shape` to
      infinity). It grows like the peak height of the density to the
      power `sum(w) - n`.]---"),
    c(r"---[under `profile = TRUE` and `REML = TRUE`, or report a
      log-likelihood of 1.4e11 under ML with only convergence
      warnings.]---",
      r"---[under `profile = TRUE` and `REML = TRUE`. Under ML it reported
      a log-likelihood of 1.4e11 with only convergence warnings.]---"),
    c(r"---[A prior exempts the model only when it provably bounds the
      dispersion, as `gamma()` on `shape` does, or a `normal()` prior on
      every coefficient of an intercept-free cell-means formula (`sigma ~
      0 + f`); a `student_t()` or `gamma(2, 1)` on `sigma`, or a prior on
      the slopes of a dispersion formula with an intercept, does not.]---",
      r"---[A prior exempts the model only when it provably bounds the
      dispersion. `gamma()` on `shape` does, and so does a `normal()`
      prior on every coefficient of an intercept-free cell-means formula
      (`sigma ~ 0 + f`). A `student_t()` or `gamma(2, 1)` on `sigma` does
      not, and neither does a prior on the slopes of a dispersion formula
      with an intercept.]---"),
    c(r"---[so the data give such an effect no curvature, and the
      approximate likelihood grew without bound (log-likelihood 449 with
      `sigma` at exp(-7.1) on 80 unweighted rows, where the exact integral
      gives -56.7).]---",
      r"---[so the data give such an effect no curvature. The approximate
      likelihood grew without bound: log-likelihood 449 with `sigma` at
      exp(-7.1) on 80 unweighted rows, where the exact integral gives
      -56.7.]---"),
    c(r"---[`loo_subsample()`**, with a refusal on a fit, as `log_lik()`,
      `loo()` and `bayes_R2()` have.]---",
      r"---[`loo_subsample()`**, with a refusal on a fit. `log_lik()`,
      `loo()` and `bayes_R2()` refuse a fit the same way.]---"),
    c(r"---[attached at "no applicable method"; `loo_moment_match(fit)`
      with loo attached]---",
      r"---[attached at "no applicable method". `loo_moment_match(fit)`
      with loo attached]---"),
    c(r"---[`geometric()`, `brmsfamily("com_poisson")`, the zero-inflated
      poisson,]---",
      r"---[`geometric()` and `brmsfamily("com_poisson")`. The others are
      the zero-inflated poisson,]---"),
    c(r"---[keeps its digits and its gradient; three are the complement of
      the other tail in part of the range, `com_poisson`'s survivor above
      `1e-3`, a hurdle family's CDF where the survivor is below one half,
      and the von Mises tails on the near side of the mean.]---",
      r"---[keeps its digits and its gradient. Three are the complement of
      the other tail in part of the range: `com_poisson`'s survivor above
      `1e-3`, a hurdle family's CDF where the survivor is below 1/2, and
      von Mises tails on the near side of the mean.]---"),
    c(r"---[Against `stats`, the closed forms, 300-bit `Rmpfr` references
      and adaptive quadrature, over grids that include both tails and the
      edges of the support, the relative error of the log of either tail
      is 1e-12 or better on 554 of 582 comparisons and 3.7e-9 or better
      on the rest, except two negative-binomial far-tail points where
      `stats::pnbinom()` is the one in error; the incomplete gamma's is
      8.04e-14]---",
      r"---[The references are `stats`, the closed forms, 300-bit `Rmpfr`
      arithmetic and adaptive quadrature, over grids that include both
      tails and the edges of the support. Against them, the relative error
      of the log of either tail is 1e-12 or better on 554 of 582
      comparisons. It is 3.7e-9 or better on the rest, except two
      negative-binomial far-tail points where `stats::pnbinom()` is the
      one in error. The incomplete gamma's is 8.04e-14]---"),
    c(r"---[A censored row is brms's likelihood: at frmtmb's estimates the
      log density of brms's Stan program with flat priors equals
      `logLik()` for every family with all four censoring codes, given the
      same data frame, and for truncated fits including a zero-inflated
      count under a constant bound and a right-censored count inside
      `trunc(lb = )`; the gradient is not asserted]---",
      r"---[A censored row is brms's likelihood. At frmtmb's estimates, the
      log density of brms's Stan program with flat priors equals
      `logLik()`, given the same data frame. This holds for every family
      with all four censoring codes, and for truncated fits. Among those
      are a zero-inflated count under a constant bound and a
      right-censored count inside `trunc(lb = )`. The gradient is not
      asserted]---"),
    c(r"---[Two exceptions are brms's: its inverse-gaussian survivor adds a
      term it should subtract (-0.466 where the survivor is -0.988 at `y =
      2`, `mu = 2`, `shape = 4`), so frmtmb's right-censored
      inverse.gaussian fit differs from brms's; and Stan's gradient of the
      gamma CDF in the shape, measured at frmtmb's estimates of one
      censored hurdle_gamma model, is -0.18 where finite differences of
      its own log density give 0.0002.]---",
      r"---[Two exceptions are brms's. Its inverse-gaussian survivor adds a
      term it should subtract: -0.466 where the survivor is -0.988 at
      `y = 2`, `mu = 2`, `shape = 4`. So frmtmb's right-censored
      inverse.gaussian fit differs from brms's. Also, Stan's gradient of
      the gamma CDF in the shape is -0.18 where finite differences of its
      own log density give 0.0002. That was measured at frmtmb's
      estimates of one censored hurdle_gamma model.]---"),
    c(r"---[and as `F(y2) - F(y)` below, so an interval in an upper tail
      keeps its digits,]---",
      r"---[and as `F(y2) - F(y)` below. So an interval in an upper tail
      keeps its digits,]---"),
    c(r"---[Fit to the same data as brms, the maximum-likelihood estimates
      lie within 0.05 posterior sd of brms's posterior means for the
      coefficients and `phi`, with and without a group term.]---",
      r"---[Fit to the same data as brms, with and without a group term,
      the maximum-likelihood estimates lie within 0.05 posterior sd of
      brms's posterior means. This holds for the coefficients and
      `phi`.]---"),
    c(r"---[Its log density is within 6e-12 relative of an exact reference
      over 384 points with `sigma` from 1e-8 to 1e6, where brms 2.23.0's
      `neg_binomial_2_log_lpmf()` is within 1.4e-10.]---",
      r"---[Over 384 points with `sigma` from 1e-8 to 1e6, its log density
      is within 6e-12 relative of an exact reference. brms 2.23.0's
      `neg_binomial_2_log_lpmf()` is within 1.4e-10.]---"),
    c(r"---[when a row's fitted distribution kept only a few percent of its
      mass inside the window:]---",
      r"---[when a row kept only a few percent of its fitted mass inside the
      window:]---"),
    c(r"---[had an objective of `-Inf` and a weibull `trunc(lb = 0.5)` fit
      with a `diag()` term stopped at a NaN gradient; that weibull model
      now fits at a log-likelihood of -145.97.]---",
      r"---[had an objective of `-Inf`. A weibull `trunc(lb = 0.5)` fit
      with a `diag()` term stopped at a NaN gradient, and now fits at a
      log-likelihood of -145.97.]---"),
    c(r"---[The objective was then rounding noise, and its gradient 21 and
      45 percent off for `skew_normal()` and `von_mises()` and eight times
      its size for `hurdle_gamma()`, at points away from the
      optimum.]---",
      r"---[The objective was then rounding noise. At points away from the
      optimum its gradient was 21 and 45 percent off for `skew_normal()`
      and `von_mises()`, and eight times its size for
      `hurdle_gamma()`.]---"),
    c(r"---[now scores -2773.28, and a Poisson fit right censored at
      75000,]---",
      r"---[now scores -2773.28. A Poisson fit right censored at
      75000,]---")),
  "extensions/frmtmb.sample/NEWS.md" = list(
    c(r"---[take the relative efficiency on one chain**, as brms 2.23.0 does
      and as `pp_check()`'s `loo_*` types already did.]---",
      r"---[take the relative efficiency on one chain.** brms 2.23.0 does
      the same, and so did the `loo_*` types of `pp_check()`.]---"),
    c(r"---[On 60-row gaussian and random-intercept models with 4000 draws,
      `loo(ndraws =)` moved by at most 0.0065 in `elpd_loo` (standard
      error 5.0) at 500, 1000 and 3000 draws.]---",
      r"---[On 60-row gaussian and random-intercept models with 4000 draws,
      `loo(ndraws =)` at 500, 1000 and 3000 draws moved by at most 0.0065
      in `elpd_loo`. Its standard error was 5.0.]---"),
    c(r"---[When every row with a positive weight has a random effect of its
      own in `mu`, the likelihood grows like `sigma^-(sum(w) - n)` as
      `sigma` shrinks to 0, and a prior that is finite at 0, as brms's
      default `student_t(3, 0, s)` on `sigma` is, leaves the posterior
      improper once `sum(w) - n` reaches 1.]---",
      r"---[Suppose that every row with a positive weight has a random
      effect of its own in `mu`. The likelihood grows like
      `sigma^-(sum(w) - n)` as `sigma` shrinks to 0. A prior that is
      finite at 0, as brms's default `student_t(3, 0, s)` on `sigma` is,
      then leaves the posterior improper once `sum(w) - n` reaches
      1.]---"),
    c(r"---[Each computes one `log_lik()` and one `loo::psis()` per call
      (per response for `loo_R2()`, per model for `loo_model_weights()`),
      with the weights of `psis()` and of `loo(save_psis = TRUE)`, on the
      same draws as the prediction.]---",
      r"---[Each computes one `log_lik()` and one `loo::psis()` per call,
      on the same draws as the prediction (per response for `loo_R2()`,
      per model for `loo_model_weights()`). The weights are those of
      `psis()` and of `loo(save_psis = TRUE)`.]---"),
    c(r"---[every `loo_epred()` and `loo_linpred()` value and every
      `loo_model_weights()` weight (a stacking pair of 0.155 and 0.845) is
      identical to brms's under R's reference BLAS, and the `loo_R2()`
      bootstrap draws]---",
      r"---[every `loo_epred()` and `loo_linpred()` value is identical to
      brms's under R's reference BLAS. So is every `loo_model_weights()`
      weight (a stacking pair of 0.155 and 0.845). The `loo_R2()`
      bootstrap draws]---"),
    c(r"---[which `log_lik()` and `posterior_epred()` take; `loo_R2()` takes
      brms's spelling only,]---",
      r"---[which `log_lik()` and `posterior_epred()` take. `loo_R2()` takes
      brms's spelling only,]---"),
    c(r"---[`draws x observations x (K - 1)` array, `posterior_epred()` and
      `posterior_predict()` return `draws x observations x K` arrays that
      read `trials()` from `newdata`, and `predictive_error()` and]---",
      r"---[`draws x observations x (K - 1)` array. `posterior_epred()` and
      `posterior_predict()` return `draws x observations x K` arrays that
      read `trials()` from `newdata`. `predictive_error()` and]---"))
)
