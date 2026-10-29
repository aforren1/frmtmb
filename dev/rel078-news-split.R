# The sentences over 25 words in the 0.78.0 NEWS sections
# (dev/rel078-news-long.R, dev/rel078-log/news-long-1.txt: core 6,
# frmtmb.sample 11), split into shorter ones, as ASD-STE100 asks.
# dev/rel078-news-merge.R applies each pair (lane text, new text) to the
# top section of the file named, before it rewraps the section; a run of
# whitespace in the lane text matches any whitespace, so a pair does not
# depend on where the lane wrapped its lines. The numbers and the facts
# are the lanes'; only the sentence boundaries and the order of clauses
# change. Adapted from dev/rel077-news-split.R.
list(
  "NEWS.md" = list(
    c(r"---[Now that row applies, so a brms `get_prior()` table passed to
      `frm()` as `prior =` gives a different MAP fit on a model with
      correlated `me()` terms.]---",
      r"---[Now that row applies. So a brms `get_prior()` table passed to
      `frm()` as `prior =` gives another MAP fit on a model with
      correlated `me()` terms.]---"),
    c(r"---[`set_prior(class = "meanme")` is a density on the latent mean of
      a `me()` term, and `class = "sdme"` one on the latent standard
      deviation itself, with the log-Jacobian of the log map, as for class
      `"sd"`.]---",
      r"---[`set_prior(class = "meanme")` is a density on the latent mean of
      a `me()` term. `class = "sdme"` is a density on the latent standard
      deviation itself, with the log-Jacobian of the log map, as for class
      `"sd"`.]---"),
    c(r"---[As in brms, the three classes take no `resp`, `dpar` or `nlpar`,
      `meanme` and `sdme` take no `group`, `corme` takes no `coef` and no
      bound, and a model without `me()` terms refuses them.]---",
      r"---[As in brms, the three classes take no `resp`, `dpar` or `nlpar`.
      `meanme` and `sdme` take no `group`, and `corme` takes no `coef` and
      no bound. A model without `me()` terms refuses them.]---"),
    c(r"---[list brms's rows for these classes: a class row and one row per
      `me()` variable for `meanme` and for `sdme`, `sdme` bounded below by
      0, and one `corme` row, plus one per `gr =` variable of a correlated
      set.]---",
      r"---[list brms's rows for these classes. `meanme` and `sdme` have a
      class row and one row per `me()` variable, and `sdme` is bounded
      below by 0. `corme` has one row, plus one per `gr =` variable of a
      correlated set.]---"),
    c(r"---[equals brms 2.23.0's Stan `log_prob` up to the Jacobian of the map
      between the two parameterizations and brms's `log(2)` normalization
      of each half-Student-t.]---",
      r"---[equals brms 2.23.0's Stan `log_prob` up to two terms. These are
      the Jacobian of the map between the two parameterizations and brms's
      `log(2)` normalization of each half-Student-t.]---"),
    c(r"---[On one `me()` term, two terms with `gr =` and a `corme` prior, and
      a term beside `(1 | h)`, at four points each, the constant agrees to
      6.4e-12 relative and the gradient to 4.3e-10 relative.]---",
      r"---[We compared three models at four points each: one `me()` term,
      two terms with `gr =` and a `corme` prior, and a term beside
      `(1 | h)`. The constant agrees to 6.4e-12 relative and the gradient
      to 4.3e-10 relative.]---")
  ),
  "extensions/frmtmb.sample/NEWS.md" = list(
    c(r"---[on two `me(gr = g)` terms, a class `corme` row edited to `lkj(4)`
      sampled `lkj(1)` (a Kolmogorov-Smirnov p-value of 1.3e-32 against
      the `lkj(4)` correlation), and the same held for `cor` and
      `sd`.]---",
      r"---[on two `me(gr = g)` terms, a class `corme` row edited to `lkj(4)`
      sampled `lkj(1)`. Its Kolmogorov-Smirnov p-value against the
      `lkj(4)` correlation was 1.3e-32. The same held for `cor` and
      `sd`.]---"),
    c(r"---[A bound on one coefficient, `coef` with `lb` or `ub`, is refused
      with brms's message, "Prior argument 'coef' may not be specified
      when using boundaries.", for the brms classes.]---",
      r"---[For the brms classes, a bound on one coefficient (`coef` with
      `lb` or `ub`) is refused with brms's message: "Prior argument 'coef'
      may not be specified when using boundaries."]---"),
    c(r"---[On two `me(gr = g)` terms with 10 levels and one row-level `me()`
      term (4 chains of 3000 draws), the posterior means of `sdme_mexg`
      and `sdme_mezg` were 12.6 and 13.6 Monte Carlo standard errors from
      brms's, and the slopes up to 8.0.]---",
      r"---[We sampled two `me(gr = g)` terms with 10 levels and one
      row-level `me()` term (4 chains of 3000 draws). The posterior means
      of `sdme_mexg` and `sdme_mezg` were 12.6 and 13.6 Monte Carlo
      standard errors from brms's, and the slopes up to 8.0.]---"),
    c(r"---[without `resp` the method returned the first response only, and
      draws of two responses taken one call each correlated at -0.05 to
      0.08 on a fit whose residual correlation was 0.82.]---",
      r"---[without `resp` the method returned the first response only.
      Draws of two responses taken one call each correlated at -0.05 to
      0.08, on a fit whose residual correlation was 0.82.]---"),
    c(r"---[are now 0.750, 0.554 and 0.315, within 1.3 Monte Carlo standard
      errors of what the draws imply and within 0.7 of brms's own draw
      from the same parameter draws.]---",
      r"---[are now 0.750, 0.554 and 0.315. They are within 1.3 Monte Carlo
      standard errors of what the draws imply, and within 0.7 of brms's
      own draw from the same parameter draws.]---"),
    c(r"---[`resp` takes several responses, in the order you name them, in
      `posterior_predict()`, `posterior_epred()`, `posterior_linpred()`,
      `predictive_error()`, `fitted()`, `predict()`, `residuals()`,
      `log_lik()`, `loo()`, `waic()` and `loo_predict()`, `loo_epred()`,
      `loo_linpred()` and `loo_predictive_interval()`.]---",
      r"---[`resp` takes several responses, in the order you name them, in
      these methods: `posterior_predict()`, `posterior_epred()`,
      `posterior_linpred()`, `predictive_error()`, `fitted()`,
      `predict()`, `residuals()`, `log_lik()`, `loo()`, `waic()`,
      `loo_predict()`, `loo_epred()`, `loo_linpred()` and
      `loo_predictive_interval()`.]---"),
    c(r"---[Each takes brms's spelling of a name (`ya` for `y_a`) as well as
      frmtmb's, refuses a name given twice, as brms does, and reads an
      empty `resp` as every response, as brms does.]---",
      r"---[Each takes brms's spelling of a name (`ya` for `y_a`) as well as
      frmtmb's. As in brms, each refuses a name given twice and reads an
      empty `resp` as every response.]---"),
    c(r"---[with the correlations in the order named; brms 2.23.0 stops there,
      because it looks up the correlations by the names in the order
      given.]---",
      r"---[with the correlations in the order named. brms 2.23.0 stops
      there, because it looks up the correlations by the names in the
      order given.]---"),
    c(r"---[Where a response predicts more than one value per row (the
      category probabilities of an ordinal or categorical response), the
      stacked layout has no place for it.]---",
      r"---[A response can predict more than one value per row: the category
      probabilities of an ordinal or categorical response. The stacked
      layout has no place for it.]---"),
    c(r"---[as a t with the shared `nu`, which is the margin of the
      multivariate t.]---",
      r"---[as a t with the shared `nu`. That t is the margin of the
      multivariate t.]---"),
    c(r"---[numbered by the group of terms that share a `gr =`, and, where a
      term has a prior of its own, one column per term,
      `prior_meanme_1_meanme_mex`, as brms's renaming spells it.]---",
      r"---[numbered by the group of terms that share a `gr =`. Where a term
      has a prior of its own, there is one column per term,
      `prior_meanme_1_meanme_mex`, as brms's renaming spells it.]---")
  )
)
