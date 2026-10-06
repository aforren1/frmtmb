## New features

* Through frmtmb 0.68.0, draws of a `cs()` model on `cumulative()`, as
  brms samples it: the per-threshold coefficients are
  `bcs_x[k]`, and a draw's `posterior_epred()` equals brms's own
  R-side density at that draw to 1.9e-16 (`dev/rel068-cs-probe.R`).
  A row whose offsets cross two thresholds has `NaN` probabilities and
  an `NA` draw, as in core.
