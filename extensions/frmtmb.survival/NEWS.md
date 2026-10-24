# frmtmb.survival (development version)

The first release. It holds the survival families that were in two
other packages, so that each has one home.

## Breaking changes

* **`cox()` and `cox_baseline()` moved here from frmtmb.** Load this
  package with `library(frmtmb.survival)`, which attaches frmtmb as
  well. `family = "cox"` and `brmsfamily("cox")` reach this package's
  `cox()` once it is loaded, and without it stop with a message that
  names it. A `cox()` fit is `identical()` to one from frmtmb 0.73.0,
  in its log likelihood, estimates and standard errors.

* **`royston_parmar()` and `rp_floored()` moved here from
  frmtmb.spline**, which is retired, with the article "How to fit a
  flexible parametric survival model". A fit is `identical()` to one
  from frmtmb.spline 0.10.2. The curve functions that read a fit,
  `frm_curve()`, `frm_curve_deriv()` and `frm_curve_feature()`, are in
  frmtmb now, and still refuse a fit whose cumulative hazard is not
  monotone, through the family's `post$curve_check` hook.

* The refusals of these functions carry the class
  `frmtmb_survival_error`. They carried `frmtmb_error` (`cox()`) and
  `frmtmb_spline_error` (`royston_parmar()`) before. All of them still
  carry `frmtmb_error`.


* frmtmb.spline's NEWS, the history of `royston_parmar()` before the
  move, is the article "History of frmtmb.spline" on this package's
  site.
