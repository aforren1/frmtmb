# The sentences over 25 words in the 0.77.0 NEWS sections
# (dev/rel077-news-long.R, dev/rel077-log/news-long-1.txt: core 18,
# frmtmb.sample 20), split into shorter ones, as ASD-STE100 asks, and the
# edits of fact the merge needs. dev/rel077-news-merge.R applies each
# pair (lane text, new text) to the top section of the file named,
# before it rewraps the section; a run of whitespace in the lane text
# matches any whitespace, so a pair does not depend on where the lane
# wrapped its lines. The numbers and the facts are the lanes'; only the
# sentence boundaries and the order of clauses change, except where a
# comment says otherwise. Adapted from dev/rel076-news-split.R.
list(
  "NEWS.md" = list(
    # an edit of fact, release fix S1 (sigmainn's review, short check):
    # get_sigma() keeps insight's default, and under sigma ~ x sigma()
    # returns NA, it does not stop
    c(r"---[`insight::get_sigma()` reads `sigma()` for a frmtmb fit.
      Before, it fell back on `sqrt(deviance / df)` where `sigma()`
      stopped.]---",
      r"---[For a current fit, `insight::get_sigma()` gives insight's
      default, as before: `sigma()` where it answers, and
      `sqrt(deviance / df)` under `sigma ~ x`, where `sigma()` returns
      `NA`. Only a saved fit with one of these terms is refused.]---"),
    c(r"---[code: on 200 rows, `y ~ mo(x) + z` with `x` in 1 and 2 has the
      log-likelihood of `lm(y ~ I(x - 1) + z)` to 12 digits, and
      `summary()` reports the weight `mox1[1]` as 1 with an error of 0, as
      brms does.]---",
      r"---[code. On 200 rows, `y ~ mo(x) + z` with `x` in 1 and 2 has the
      log-likelihood of `lm(y ~ I(x - 1) + z)` to 12 digits. `summary()`
      reports the weight `mox1[1]` as 1 with an error of 0, as brms
      does.]---"),
    c(r"---[Its coefficient is not determined by the data, so its standard
      error is `NaN`, with the usual warning; it is not counted in the
      degrees of freedom of `logLik()`, `AIC()` and `anova()`, as `lm()`
      does not count a column it aliases, and `predict()` holds it at its
      estimate.]---",
      r"---[Its coefficient is not determined by the data, so its standard
      error is `NaN`, with the usual warning. It is not counted in the
      degrees of freedom of `logLik()`, `AIC()` and `anova()`, as `lm()`
      does not count a column it aliases. `predict()` holds it at its
      estimate.]---"),
    c(r"---[Now a matrix over `get_coef()` (matched by its names), a matrix
      over the whole outer parameter vector (refused for a REML fit, whose
      vector has no fixed effects), a function, `FALSE`, `"HC0"`, `"HC1"`
      (and `"stata"`, which marginaleffects maps to HC1) and a one-sided
      cluster formula work.]---",
      r"---[Now a matrix over `get_coef()`, matched by its names, works. So
      do a matrix over the whole outer parameter vector, a function,
      `FALSE`, `"HC0"`, `"HC1"` and a one-sided cluster formula.
      marginaleffects maps `"stata"` to HC1. The whole-vector matrix is
      refused for a REML fit, whose vector has no fixed effects.]---"),
    # an edit of fact, release fix S2 (mointerop's review, short check)
    c(r"---[`"HC2"` to `"HC5"` and `"robust"` need leverages,]---",
      r"---[`"HC2"` to `"HC5"` and `"robust"` need leverages, which
      frmtmb does not compute,]---"),
    c(r"---[On a fit whose `mo()` term has a category that no row is in,
      `get_coef()` carries that simplex as its cumulative weights at the
      observed categories (`zeta1_cum1`, ...), not as its chart
      coordinates.]---",
      r"---[Take a fit whose `mo()` term has a category that no row is in.
      There `get_coef()` carries that simplex as its cumulative weights at
      the observed categories (`zeta1_cum1`, ...), not as its chart
      coordinates.]---"),
    c(r"---[On 200 rows with values 1, 2 and 5, `predictions()` at 1, 2 and
      5 now gives 0.1222210, 0.1135415 and 0.1212624, the standard errors
      of `lm()` on `factor(x)` with the maximum-likelihood sigma, and with
      `vcov = "HC0"` 0.1224166, 0.1103619 and 0.1248376, `lm()`'s
      HC0.]---",
      r"---[On 200 rows with values 1, 2 and 5, `predictions()` at 1, 2 and
      5 now gives 0.1222210, 0.1135415 and 0.1212624. These are the
      standard errors of `lm()` on `factor(x)` with the maximum-likelihood
      sigma. With `vcov = "HC0"` it gives 0.1224166, 0.1103619 and
      0.1248376, `lm()`'s HC0.]---"),
    c(r"---[A standard deviation is class `sd` (a `gp()` term's `sdgp`, and
      `lscale` for its length scales), a correlation that an LKJ density
      fits is class `cor`, and a residual correlation is class
      `rescor`.]---",
      r"---[A standard deviation is class `sd` (a `gp()` term's `sdgp`, and
      `lscale` for its length scales). A correlation that an LKJ density
      fits is class `cor`, and a residual correlation is class
      `rescor`.]---"),
    c(r"---[On `(1 + x + z | g)`, `lkj(2)` on class `cor` and then
      `normal(0, 1)` on `theta_4` kept the normal on one correlation
      parameter and no density on the other two, and a MAP fit took those
      two to 559 and 1741.]---",
      r"---[On `(1 + x + z | g)`, `lkj(2)` on class `cor` and then
      `normal(0, 1)` on `theta_4` kept the normal on one correlation
      parameter and no density on the other two. A MAP fit took those two
      to 559 and 1741.]---"),
    c(r"---[Class `theta` keeps the parameters that no brms class names, such
      as a `toep()` correlation or a parameter of a structure that an
      extension registers (a `car()` dependence); with no `coef` it covers
      only those.]---",
      r"---[Class `theta` keeps the parameters that no brms class names.
      Examples are a `toep()` correlation and a parameter of a structure
      that an extension registers (a `car()` dependence). With no `coef` it covers
      only those.]---"),
    c(r"---[Class `Intercept` with `nlpar` is stored as brms's slot of the
      same parameter, class `b` with `coef = "Intercept"`, so the two
      spellings together are refused as a duplicate.]---",
      r"---[Class `Intercept` with `nlpar` is stored as brms's slot of the
      same parameter, class `b` with `coef = "Intercept"`. So the two
      spellings together are refused as a duplicate.]---"),
    c(r"---[a `(flat)` row with `lb` and `ub` is a box with no density, as
      every bounds-only specification is there, and a box at the ends of
      the interval bounds nothing.]---",
      r"---[a `(flat)` row with `lb` and `ub` is a box with no density, as
      every bounds-only specification is there. A box at the ends of the
      interval bounds nothing.]---"),
    c(r"---[The prior object that `prior_summary()` returns for a
      `frm_sample()` run carries the densities, and passed to `frm()`
      still gives 0.7249, the mode of the sampled density.]---",
      r"---[The prior object that `prior_summary()` returns for a
      `frm_sample()` run carries the densities. Passed to `frm()`, it
      still gives 0.7249, the mode of the sampled density.]---"),
    c(r"---[The named list of prior objects,
      `prior = list(theta_1 = prior_normal(0, 1))`, which put densities on
      internal parameters by name, is refused by `frm()`, `update()`,
      `validate_prior()` and `frm_simulate()`, with how to write the same
      density with `set_prior()`.]---",
      r"---[The named list of prior objects,
      `prior = list(theta_1 = prior_normal(0, 1))`, put densities on
      internal parameters by name. `frm()`, `update()`, `validate_prior()`
      and `frm_simulate()` refuse it, and say how to write the same
      density with `set_prior()`.]---"),
    # an edit of fact, release fix R1 (reslot's review, short check of
    # fix round 2): the coerced table is marked too
    c(r"---[A priorlist read from a `default_prior(route = "fit")` table
      keeps that route through `+`, `c()` and `as.brmsprior()`, so
      frmtmb.sample can refuse it; `frm()` drops the mark from the prior
      it keeps.]---",
      # and the release review's finding 1: a row or column selection
      # and subset() keep the route ([.frmtmb_prior_rows)
      r"---[A `default_prior(route = "fit")` table keeps that route through
      a row or column selection and `subset()`. A priorlist read from it
      keeps it through `+`, `c()` and `as.brmsprior()`. So frmtmb.sample
      can refuse it. `as.data.frame()` drops the route. So
      `as.brmsprior()` marks a data frame that holds frmtmb's printed
      `(flat)`. `frm()` drops either mark from the prior it keeps.]---"),
    c(r"---[`c()` of a prior specification with something that is not one (a
      density object, a named list of them) is refused with the
      `set_prior()` spelling to write, where it built a prior that `frm()`
      then stopped on with an internal error.]---",
      r"---[`c()` of a prior specification with something that is not one,
      such as a density object or a named list of them, is refused. The
      message gives the `set_prior()` spelling to write. Before, it built a prior that
      `frm()` then stopped on with an internal error. `c()` now reads a
      `default_prior()` table or a brms prior as `as.brmsprior()` reads
      it, with its route. That holds when a specification or a frmtmb
      table comes first.
      brms's `c()` combines the rows the same way. With a brms table
      first, brms's own `c()` refuses the frmtmb specification.]---"),
    # an edit of fact, the release review's finding 2: the sandwich
    # spellings under a residual correlation term
    c(r"---[`"kenward-roger"` denominator degrees of freedom. See
      `?frmtmb-marginaleffects`.]---",
      r"---[`"kenward-roger"` denominator degrees of freedom. Under a
      residual correlation term, such as `ar(cov = TRUE)`, `"HC0"`,
      `"HC1"`, `"stata"` and a cluster formula are refused too. Its density
      is joint over each group. `vcov = TRUE` works there. See
      `?frmtmb-marginaleffects`.]---"),
    c(r"---[A prior density can be brms's `(flat)`, which the prior resolver
      reads as no density on the parameter a class names, with the
      Jacobian of the map to the internal parameter, and an LKJ density can
      be flat on the Cholesky factor. frmtmb.sample's `prior = "flat"` uses
      both; the frmtmb.sample of the same release needs this frmtmb.]---",
      r"---[A prior density can be brms's `(flat)`. The prior resolver
      reads it as no density on the parameter a class names, with the
      Jacobian of the map to the internal parameter. An LKJ density can be
      flat on the Cholesky factor. frmtmb.sample's `prior = "flat"` uses
      both, so frmtmb.sample 0.24.0 needs frmtmb 0.77.0.]---"),
    c(r"---[With `prior_only = TRUE` the objective leaves out every response
      density and keeps the densities of the group-level coefficients, of
      the `me()` latent values and their measurements, and of a `mi(sdy)`
      response's measurements, which are the terms brms 2.23.0 keeps
      outside `if (!prior_only)` under `brm(sample_prior = "only")`.]---",
      r"---[With `prior_only = TRUE` the objective leaves out every response
      density. It keeps the densities of the group-level coefficients, of
      the `me()` latent values and their measurements, and of a `mi(sdy)`
      response's measurements. These are the terms brms 2.23.0 keeps
      outside `if (!prior_only)` under `brm(sample_prior = "only")`.]---"),
    c(r"---[marginaleffects 1.0.0 still warns that these arguments are "not
      known to be supported" for this class, because it keeps a fixed list
      of known arguments per class; they have their effect.]---",
      r"---[marginaleffects 1.0.0 still warns that these arguments are "not
      known to be supported" for this class. It keeps a fixed list of
      known arguments per class. The arguments have their effect.]---"),
    c(r"---[`vcov_cluster()` on a fit that lost a parameter (a `mo()`
      category no row is in, a variance at its boundary) gave `NaN` in
      every entry, because the bread held `NaN` in the lost rows.]---",
      r"---[`vcov_cluster()` gave `NaN` in every entry on a fit that lost a
      parameter (a `mo()` category no row is in, a variance at its
      boundary). The bread held `NaN` in the lost rows.]---")
  ),
  "extensions/frmtmb.sample/NEWS.md" = list(
    c(r"---[A standard deviation, `sigma` and the other positive
      distributional parameters are now flat on their own scale, a
      correlation matrix (`cor`, `rescor`, `cortime`, `lncor`) flat on its
      Cholesky factor, and the ordered thresholds and `delta` flat on
      themselves; the defaults that are already flat on brms's parameter
      (`ar`, `ma`, `cosy`, `lagsar`, and `dirichlet(1)` on a simplex, a
      mixture's weights included) stay.]---",
      r"---[A standard deviation, `sigma` and the other positive
      distributional parameters are now flat on their own scale. A
      correlation matrix (`cor`, `rescor`, `cortime`, `lncor`) is flat on
      its Cholesky factor, and the ordered thresholds and `delta` are flat
      on themselves. The defaults that are already flat on brms's
      parameter stay: `ar`, `ma`, `cosy`, `lagsar`, and `dirichlet(1)` on
      a simplex, a mixture's weights included.]---"),
    c(r"---[every slot, 4 chains of 3000 draws each: on `y ~ 1 + (1 | g)`
      with 6 groups the posterior mean of `sd_g__Intercept` is 1.458 in
      brms, 1.455 now and was 1.197; on `y ~ x + (1 + x | g)` with 8 groups
      the correlation's mean is -0.283 in brms, -0.265 now, and was -1.000
      (sd below 0.0001, 529 divergent transitions); on `ar(cov = TRUE)`
      with 6 series of 5 the mean of `ar[1]` is 0.481 in brms, 0.480 now
      and was 0.541.]---",
      r"---[every slot, 4 chains of 3000 draws each. On `y ~ 1 + (1 | g)`
      with 6 groups the posterior mean of `sd_g__Intercept` is 1.458 in
      brms, 1.455 now and was 1.197. On `y ~ x + (1 + x | g)` with 8
      groups the correlation's mean is -0.283 in brms, -0.265 now, and was
      -1.000 (sd below 0.0001, 529 divergent transitions). On
      `ar(cov = TRUE)` with 6 series of 5 the mean of `ar[1]` is 0.481 in
      brms, 0.480 now and was 0.541.]---"),
    c(r"---[The `prior_summary()` of a `prior = "flat"` run lists the
      intercepts as `(flat)` too, so passed back as `prior =`, as the
      object or as its table through `as.brmsprior()`, it gives the same
      density (log density and gradient to 1e-12 relative at 40
      points).]---",
      r"---[The `prior_summary()` of a `prior = "flat"` run lists the
      intercepts as `(flat)` too. Passed back as `prior =`, as the object
      or as its table through `as.brmsprior()`, it gives the same density:
      log density and gradient to 1e-12 relative at 40 points.]---"),
    c(r"---[It took over the internal parameters it named after the class
      priors, and a name inside a correlation block removed the default
      `lkj(1)` from the whole block:]---",
      r"---[It took over the internal parameters it named after the class
      priors. A name inside a correlation block removed the default
      `lkj(1)` from the whole block:]---"),
    c(r"---[brms declares it on (0, 1) and leaves it flat, and here it carries
      no density, which is flat on its logit and not integrable at 0.]---",
      r"---[brms declares it on (0, 1) and leaves it flat. Here it carries
      no density, which is flat on its logit and not integrable at 0.]---"),
    c(r"---[`prior =` refuses a `default_prior()` or `validate_prior()` table
      of `route = "fit"`, the default route, edited or not, and a priorlist
      made from one with `as.brmsprior()`, `+` or `c()`.]---",
      r"---[`prior =` refuses a `default_prior()` or `validate_prior()` table
      of `route = "fit"`, the default route, edited, cut down to some rows
      or columns, or not. It refuses a
      priorlist made from one with `as.brmsprior()`, `+` or `c()`
      too.]---"),
    c(r"---[Its `(flat)` rows are `frm()`'s defaults; read as brms's `(flat)`
      they replaced brms's defaults on the intercept, sigma, the standard
      deviations and the correlations without a word.]---",
      r"---[Its `(flat)` rows are `frm()`'s defaults. Read as brms's
      `(flat)`, they replaced brms's defaults on the intercept, sigma, the
      standard deviations and the correlations without a word.]---"),
    # an edit of fact, release fix R1 (reslot's review, short check of
    # fix round 2): the table through as.data.frame()
    c(r"---[A table built by hand, or brms's own, is read as written.]---",
      r"---[`as.data.frame()` drops a table's route. So `as.brmsprior()`
      marks a data frame that holds frmtmb's printed `(flat)`. The
      priorlist is refused where its
      `(flat)` rows would replace one of brms's proper defaults. A
      `route = "sample"` table never does that. A table built by hand in
      brms's spelling, with `""` for `(flat)`, is read as written, and so
      is brms's own.]---"),
    c(r"---[`prior_summary()` and the default-prior message print brms's
      `(flat)` with `lb` and `ub` for a first-order `cov = TRUE` `ar` or
      `ma`, for `cosy` and for `lagsar`, where they printed
      `uniform(-1, 1)`.]---",
      r"---[`prior_summary()` and the default-prior message print brms's
      `(flat)` with `lb` and `ub` for `cosy`, for `lagsar` and for a
      first-order `cov = TRUE` `ar` or `ma`. They printed
      `uniform(-1, 1)`.]---"),
    c(r"---[`hypothesis()` divides the kernel density of the posterior draws
      at the point by the kernel density of the prior draws there (brms's
      `evidence_ratio()`), and returns the prior draws of each hypothesis
      as `prior_samples`.]---",
      r"---[`hypothesis()` divides the kernel density of the posterior draws
      at the point by the kernel density of the prior draws there, as
      brms's `evidence_ratio()` does. It returns the prior draws of each
      hypothesis as `prior_samples`.]---"),
    c(r"---[Before, the denominator was the prior density read off the prior
      specification, so draws made without prior draws gave a number where
      brms gives `NA`; and `sigma`, a standard deviation, a correlation, a
      nonlinear function and the intercept were refused with a warning,
      where brms gives a ratio for all but the intercept.]---",
      r"---[Before, the denominator was the prior density read off the prior
      specification. So draws made without prior draws gave a number where
      brms gives `NA`. And `sigma`, a standard deviation, a correlation, a
      nonlinear function and the intercept were refused with a warning.
      brms gives a ratio for all but the intercept.]---"),
    c(r"---[Against brms 2.23.0 on `y ~ x + (1 | g)` (4 chains of 6000 a
      side, seed 11), the ratios of `x = 0`, `x = 0.2`, `sigma = 1` and
      `sd_g__Intercept = 0.5` are within 1.7 Monte Carlo standard errors
      of brms's.]---",
      r"---[We compared with brms 2.23.0 on `y ~ x + (1 | g)` (4 chains of
      6000 a side, seed 11). The ratios of `x = 0`, `x = 0.2`,
      `sigma = 1` and `sd_g__Intercept = 0.5` are within 1.7 Monte Carlo
      standard errors of brms's.]---"),
    c(r"---[`frm_sample()` refuses a `mo()` term on one value, which is 0 on
      every row, so its coefficient has an improper posterior under the
      default flat prior.]---",
      r"---[`frm_sample()` refuses a `mo()` term on one value. That term is
      0 on every row, so its coefficient has an improper posterior under
      the default flat prior.]---"),
    c(r"---[as `brm(sample_prior = "yes")` does: one column per prior
      statement of brms's Stan program for the model (`prior_b` for a
      class-wide prior, `prior_b_x` for each coefficient of a vector with a
      coefficient-specific prior, `prior_Intercept`, `prior_sigma`,
      `prior_sd_g`, `prior_cor_g`, `prior_simo_mox1[k]`, ...), none for a
      flat prior, a bounded parameter drawn truncated to its bounds.]---",
      r"---[as `brm(sample_prior = "yes")` does. That is one column per
      prior statement of brms's Stan program for the model: `prior_b` for
      a class-wide prior, `prior_b_x` for each coefficient of a vector with
      a coefficient-specific prior, `prior_Intercept`, `prior_sigma`,
      `prior_sd_g`, `prior_cor_g`, `prior_simo_mox1[k]`, and so on. A flat
      prior has no column. A bounded parameter is drawn truncated to its
      bounds.]---"),
    c(r"---[On 19 models (population-level and monotonic coefficients,
      intercepts of the location and of a distributional parameter,
      nonlinear parameters, `sigma`, `nu`, `phi`, `hu`, smooth standard
      deviations in two predictors, group-level standard deviations,
      correlations, `rescor`, `mo()` simplexes, mixture weights, ordinal
      thresholds, `ar` of two responses, multivariate and bounded priors),
      the names are brms 2.23.0's, in brms's order.]---",
      r"---[We checked 19 models. They cover population-level and monotonic
      coefficients, intercepts of the location and of a distributional
      parameter, nonlinear parameters, `sigma`, `nu`, `phi` and `hu`.
      They also cover smooth standard deviations in two predictors,
      group-level standard deviations, correlations, `rescor` and `mo()`
      simplexes. And they cover mixture weights, ordinal thresholds, `ar`
      of two responses, and multivariate and bounded priors. On each, the names are brms 2.23.0's, in brms's
      order.]---"),
    c(r"---[generated quantities: on 39 columns of 9 models the one-sample
      Kolmogorov-Smirnov p-value against the exact (truncated) prior is
      between 0.125 and 0.986, and against brms's own draws of two models
      (14 columns, 4000 draws a side) between 0.033 and 0.962.]---",
      r"---[generated quantities. On 39 columns of 9 models the one-sample
      Kolmogorov-Smirnov p-value against the exact (truncated) prior is
      between 0.125 and 0.986. Against brms's own draws of two models (14
      columns, 4000 draws a side), it is between 0.033 and 0.962.]---"),
    c(r"---[`"only"` samples the prior of the whole model with the likelihood
      left out, as brms does, so `posterior_predict()`, `pp_check()` and
      `conditional_effects()` give the prior predictive distribution.]---",
      r"---[`"only"` samples the prior of the whole model with the
      likelihood left out, as brms does. So `posterior_predict()`,
      `pp_check()` and `conditional_effects()` give the prior predictive
      distribution.]---"),
    c(r"---[Unlike brms, a vector with some coefficients flat is refused too,
      and so is a model with missing values of an `mi()` term, which have
      no density without the likelihood.]---",
      r"---[Unlike brms, a vector with some coefficients flat is refused
      too. So is a model with missing values of an `mi()` term, which have
      no density without the likelihood.]---"),
    c(r"---[On `y ~ x + (1 | g)` with a `normal(0, 1)` prior on `b` (4 chains
      of 4000), every mean and 5, 50 and 95 percent quantile of
      `b_Intercept`, `b_x`, `sigma`, `sd_g__Intercept`, `r_g[1,Intercept]`
      and of a prior predictive draw is within 3.0 Monte Carlo standard
      errors of brms's, and at two more seeds every quantile of `b_x`,
      `sigma` and `sd_g__Intercept` is within 1.3 of the exact prior's
      (brms's: 2.4).]---",
      r"---[We checked `y ~ x + (1 | g)` with a `normal(0, 1)` prior on
      `b` (4 chains of 4000). We compared each mean and each 5, 50 and 95
      percent quantile of `b_Intercept`, `b_x`, `sigma`,
      `sd_g__Intercept`, `r_g[1,Intercept]` and of a prior predictive
      draw. Each is within 3.0 Monte Carlo standard errors of brms's. At two more seeds every
      quantile of `b_x`, `sigma` and `sd_g__Intercept` is within 1.3 of
      the exact prior's (brms's: 2.4).]---"),
    c(r"---[A name gets the prior column whose name is the longest match at
      its start, in random order (`b_x` reads `prior_b`); the overall
      intercept `b_Intercept` is left out with brms's warning.]---",
      r"---[A name gets the prior column whose name is the longest match at
      its start, in random order (`b_x` reads `prior_b`). The overall
      intercept `b_Intercept` is left out with brms's warning.]---"),
    c(r"---[The warmup iterations rstan saves (`save_warmup = TRUE`, the
      default of rstan and of brms) stay in the draws' stanfit, as in a
      brmsfit, and `inc_warmup = TRUE` on]---",
      r"---[The warmup iterations rstan saves (`save_warmup = TRUE`, the
      default of rstan and of brms) stay in the draws' stanfit, as in a
      brmsfit. `inc_warmup = TRUE` on]---"),
    c(r"---[`rhat()` and `neff_ratio()` converts them to the draws' variables
      and puts them before each chain's draws, as brms does.]---",
      r"---[`rhat()` and `neff_ratio()` converts them to the draws'
      variables. It puts them before each chain's draws, as brms
      does.]---")
  )
)
