The main changes:

- BREAKING: `mo()` codes an integer variable from its minimum, as brms
  2.23.0 does. Every fit whose integer `mo()` variable does not start
  at 0 changes. A fit saved by 0.75.0 or earlier with such a variable
  is refused.
- BREAKING: each residual correlation parameter has one prior slot, and
  the prior tables bound the residual correlation classes as brms does.
  Under a first-order `cov = TRUE` `ar()`, `ma()` or `arma()`, a prior
  on class `sigma` is read at brms's innovation sd.
- The standard error of a prediction includes the uncertainty of the
  `mo()` simplex.
- New in the sampling API: `autocor_flat_priors()`,
  `autocor_bound_chart()`, `autocor_innov_sd()` and
  `autocor_row_loglik()`. New seam slot:
  `frmtmb_register_autocor(jacobian =)`.
- Update the extensions together with frmtmb. frmtmb.sample 0.23.0 and
  frmtmb.spatial 0.3.0 need frmtmb 0.76.0.
