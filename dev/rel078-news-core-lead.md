The main changes:

- BREAKING: brms's prior classes `meanme`, `sdme` and `corme` of a
  `me()` term are honored, as in brms 2.23.0. A brms prior table passed
  to `frm()` now applies its `corme` row.
- BREAKING: in a `default_prior(route = "sample")` table, a group row
  prints brms's `(flat)`, so an edited class row reaches the groups.
- New sampling-API export `predict_rescor_draw()`, the joint predictive
  draw of a `set_rescor(TRUE)` model, and a `resp` argument to
  `rescor_row_loglik()`.
- Update the extensions together with frmtmb. frmtmb.sample 0.25.0
  needs frmtmb 0.78.0.
