# frmtmb.survival 0.2.0

Needs frmtmb 0.75.0, whose `frmtmb_register_aterm(option =)` this
package calls when it loads.

## New features

* **brms's `bhaz()` addition term.** `time | cens(c) + bhaz(df = 6) ~ x`
  with `family = cox()` is the model of `cox(df = 6)`: the two fits have
  `identical()` log likelihoods, estimates and covariance matrices.
  `bhaz()` takes brms 2.23.0's arguments `df`, `degree`, `intercept`,
  `knots` and `Boundary.knots`, and sets the argument of the same name
  of `cox()`. Its basis equals the `Zbhaz` and `Zcbhaz` of brms's
  `standata()` to 6.7e-16 of their size or less (the largest difference
  over the largest entry). At each of 1000 posterior draws of three brms
  fits, the log likelihood equals the sum of brms's `log_lik()` to
  1.5e-16 of its size.

* **`cox(knots =, Boundary.knots =)`** set the internal and the boundary
  knots of the baseline basis, which only `bhaz()` could set in brms.
  With `knots`, `df` is `length(knots) + degree + intercept` and need
  not be given.

### Where `bhaz()` differs from brms

* An argument that `bhaz()` and `cox()` both set must have one value.
  `bhaz(df = 6)` with `family = cox(df = 5)` stops and names the two.
  `df` and `knots` that disagree stop too, in either place: brms ignores
  `df` when `knots` is given. The basis is checked after the two are
  merged, so `cox(df = 3)` with `bhaz(intercept = FALSE)` fits.

* `bhaz(gr = g)`, brms's stratified Cox model with one baseline per
  level of `g`, stops. The message says that `cox()` has one baseline.
  It tells you to fit each stratum as a model of its own, or to add `g`
  to the formula. An unnamed first argument, `bhaz(g)`, is `gr` as in
  brms; `bhaz(6)`, which brms reads as one stratum named 6, stops and
  suggests `bhaz(df = 6)`.

* A shortened argument name that brms reads is read the same way:
  `bhaz(d = 6)` is `df = 6`, `kn =` is `knots`, `deg =` is `degree`. A
  name that brms drops without a message stops and names the argument
  that was probably meant: `bhaz(int = FALSE)` and `bhaz(inter = FALSE)`
  (brms keeps the intercept), `bhaz(Bound = c(0, 30))` (brms keeps the
  default boundary knots). Other names brms does not document, such as
  `bhaz(foo = 1)`, stop too. `warn.outside` is taken and has no effect,
  as in brms.

* Boundary knots that do not cover every observed time stop. brms
  extrapolates the basis past them. With `Boundary.knots = c(0.5, 10)`
  on times from 0.0044 to 14.4, 62 rows lie below the lower knot and 1
  above. At equal weights, brms's cumulative baseline hazard reaches -1.08
  there, so the survivor function exceeds one. The hazard is negative
  for some weights.

## Other changes

* `cox(bhaz = list(df = 6))`, the spelling of brms before 2.22.0, stops
  with the spellings to use instead. Before, it stopped with "unused
  argument". brms 2.23.0's `cox()` has no `bhaz` argument either.

* A `df` below `degree + intercept` is still refused when the model is
  fitted, and the message now gives the right minimum without an
  intercept: "df must be at least degree (3) without an intercept".

# frmtmb.survival 0.1.0

The first release. It holds the survival families that were in two
other packages, so that each has one home. Needs frmtmb 0.74.0, whose
`frmtmb_family()` takes `simplex_pars =` and `aterm_data =`.

## Breaking changes

* **`cox()` and `cox_baseline()` moved here from frmtmb.** Load this
  package with `library(frmtmb.survival)`, which attaches frmtmb as
  well. `family = "cox"` and `brmsfamily("cox")` reach this package's
  `cox()` once it is loaded, and without it stop with a message that
  names it. A `cox()` fit is `identical()` to one from frmtmb 0.73.0, in
  its log likelihood, estimates and standard errors.

* **`royston_parmar()` and `rp_floored()` moved here from
  frmtmb.spline**, which is retired, with the article "How to fit a
  flexible parametric survival model". A fit is `identical()` to one
  from frmtmb.spline 0.10.2. The curve functions that read a fit,
  `frm_curve()`, `frm_curve_deriv()` and `frm_curve_feature()`, are in
  frmtmb now. They still refuse a fit whose cumulative hazard is not
  monotone, through the family's `post$curve_check` hook.

* The refusals of these functions carry the class
  `frmtmb_survival_error`. They carried `frmtmb_error` (`cox()`) and
  `frmtmb_spline_error` (`royston_parmar()`) before. All of them still
  carry `frmtmb_error`.

* frmtmb.spline's NEWS, the history of `royston_parmar()` before the
  move, is the article "History of frmtmb.spline" on this package's
  site.
