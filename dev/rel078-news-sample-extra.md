## Breaking changes

* Where `a ~ 1`, `set_prior(class = "Intercept", nlpar = "a")` and
  brms's class row `set_prior(class = "b", nlpar = "a")` are one prior
  slot. With both, the `Intercept` spelling's density stands, in either
  order, as `frm()` and brms resolve them. A bound on either row is the
  class bound, the `Intercept` spelling's where both set one. The fit's
  own prior and the call's `prior =` combine in the same way.
  `prior_summary()` of the draws lists both rows as they were written,
  as brms does.

## Bug fixes

* `bayes_R2()` and `loo_R2()` refuse a response named twice in `resp`,
  as the other methods do. Before, `bayes_R2()` gave two identical rows,
  and `loo_R2()` gave two rows with different values.
