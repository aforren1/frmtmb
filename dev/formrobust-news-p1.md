# frmtmb (development version)

## Breaking changes

* **`update()` with a complete formula keeps the stored parameter
  formulas**, as brms's `update.brmsfit()` keeps them through
  `update.brmsformula()`. `update(fit, y ~ x2)` on `bf(y ~ x, sigma ~ z)`
  fitted `y ~ x2` with a constant `sigma`; it now fits
  `bf(y ~ x2, sigma ~ z)`. A formula of the update replaces a stored one
  of the same name, with brms's message. On a nonlinear fit the update
  replaces the body and keeps the nonlinear parameters' formulas, with
  brms's message "Argument 'formula.' will completely replace the
  original formula in non-linear models", so
  `update(fit, bf(count ~ a + b, nl = TRUE))` now refits where it was
  refused (ledger row `brmsfit-methods:955`). The updated fit's call
  holds the pooled formula as a `frmtmb::bf()` call. A complete formula
  on a multivariate fit is refused, as brms refuses it ("Updating
  formulas of multivariate models is not yet possible"); it replaced
  every response's formulas before.
* **`bf(nl = TRUE)` without a parameter formula is no longer refused at
  `bf()`**; `frm()` refuses it, as brms refuses it at `brm()`. The
  formulas may arrive afterwards, `bf(y ~ a * x, nl = TRUE) + lf(a ~ 1)`,
  or from the fit `update()` keeps.
* **A `bernoulli()` response can hold any two values.** They are coded
  0 and 1 by level order, as brms codes them (`-2`, `-1` become 0, 1;
  so do `1`, `2`, `0.1`, `0.9`, `"no"`, `"yes"` and `FALSE`, `TRUE`).
  The coding is kept with the fit and applied to the response of
  newdata in `residuals()`, `pp_check()` and frmtmb.sample's
  `predictive_error()`. Before, anything but 0/1 was refused (ledger
  row `standata:75`). Two values that both lie strictly between 0 and 1
  are coded too, with a warning that they look like a proportion, which
  `binomial()` with `trials()` or `Beta()` fits (lme4#682); effect
  coding such as -0.5 and 0.5 does not warn. `refit()` takes the 0/1
  codes that `simulate()` returns, and refuses any other value by name.
  One divergence: brms codes an all-`TRUE` logical response 0, and
  frmtmb codes it 1.
* **An addition term's expression that gives `NA` is refused by name**,
  as brms refuses it at `standata()`: `trunc(lb = ifelse(x > 1, NA,
  -5))` and `trials(ifelse(x > 1, NA, n))`. At 0.66.0 the expression
  was a column of the model frame, and `na.omit` dropped those rows.

## Formula grammar

* **An addition term takes an expression**, as brms evaluates it on
  the data rows: `weights(wt * 2)`, `rate(time * 2)`, `se(s / 2)`,
  `trunc(lb = lb - 1)`, `subset(x > 0)`, `index(id * 3)`. Before, the
  expression went into the model frame as a formula term, where
  `wt * 2` is an interaction, and the fit stopped with "invalid model
  formula in ExtractVars". A single value is used for every row, so
  `trunc(lb = min(y) - 1)` works too, and a value of any other length is
  refused by name. Fit, `predict()`, `fitted()`, `simulate()` and
  `conditional_effects()` read the expression the same way, in a
  multivariate model too. A constant of the formula environment is
  read as the constant it is, bare or in an expression: `trials(k)`,
  `trials(k + 0)` and `weights(wt * k)` with `k <- 10` all work (brms
  refuses all three, since it reads variables from the data alone). A
  name that R finds only as a function, as `t` in `weights(t * 2)`
  without a column `t`, is refused by name. The interval bound of
  `cens(x, y2)` takes one value per row, as brms requires.
* **`weights(w, scale = TRUE)`** scales the weights to a mean of one on
  the rows the response reads, as brms does.
* **brms's deprecated `0 + intercept`** is accepted with brms's warning
  ("Reserved variable name 'intercept' is deprecated") and read as brms
  reads it, a data column of ones whose coefficient is `intercept`.
  Newdata needs no such column (ledger rows `standata:970`, `:974`).
* **`ar()`, `ma()`, `arma()`, `cosy()` and `unstr()` exist as functions**
  that return the term, as in brms. Added to a formula,
  `bf(y ~ 1) + arma(x)`, the object is refused with brms's message
  (ledger row `brm:112`). brms's two ways to write the terms apart from
  the formula work: `bf(y ~ x, autocor = ~ ar(t, g))` and
  `bf(y ~ x) + acformula(~ ar(t, g))`.

## New arguments and answers

* **`frm(drop_unused_levels = )`**, brms's argument. With `FALSE` an
  unused level of a factor predictor gives a column of zeros, which
  frmtmb drops with its rank-deficiency message (ledger row
  `standata:1133`: brms keeps the column).
* **`autocor()` on a fit**, brms's deprecated accessor, returns `NULL`
  with brms's deprecation warning, as brms does for every fit since
  2.11 (ledger rows `brmsfit-methods:112`, `:113`).
* **`predict(newdata = )` under `ar()`/`arma()` with `cov = FALSE`
  answers when the response is `NA`** in some rows, or absent. As in
  brms's `.predictor_arma()`, the recursion runs in each group's time
  order and a missing response is filled with a draw from the family
  at its shifted mean, so the rows after it read that draw's residual.
  Before, the call was refused by name (ledger row
  `brmsfit-methods:747`). **`fitted()` diverges from brms here, on
  purpose:** brms fills with a draw per posterior draw, so its
  `fitted()` carries the fill's spread, and a maximum likelihood fit has
  no draws to carry the fill through. `fitted()` fills with the
  expected value, which is the mean of brms's fill at fixed parameters,
  and its `Est.Error` is the parameter uncertainty alone (see
  `?fitted.frmtmb_fit`). frmtmb.sample's `posterior_epred()` fills with
  draws, as brms does.

## Bug fixes

* **`conditional_effects()` and `emmeans()` on a model with an
  `offset()`** stopped with "non-numeric argument to mathematical
  function" and "undefined columns selected": the model frame held the
  column `offset(log(time))` and not `time`, so the grids had no
  `time`. The frame now holds the offset's variables, and both grids
  hold them at their means, as brms's do. `emmeans()` includes a
  predictor's offset at the grid's value, `log(mean(time))` or the
  value `at =` sets, which is emmeans's own `.offset.` column and what
  brms's `emmeans()` gives. An offset inside a nonlinear parameter's
  formula is left out of the body's emmean, as brms leaves it out, and
  `epred = TRUE` includes every offset once. An offset's variable is
  not a default display of `conditional_effects()`, as in brms, and a
  factor read inside an offset, as in `offset(log(as.numeric(f)))`, no
  longer lets model.frame()'s "is not a factor" warning escape.

## Extension API

* `arma_cond_fill_dpars()` gives a response's parameters on newdata
  for a predictive draw, filling a missing `cov = FALSE` response as
  brms does, `arma_cond_fill_epred()` gives one draw's expected
  response with that fill, and `response_codes_newdata()` codes a
  response read from newdata as the fit coded its own. All three are on
  `?frmtmb-sampling-api`, for frmtmb.sample.
