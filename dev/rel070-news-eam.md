# frmtmb.eam 0.11.3

Needs frmtmb 0.70.0. No change to the package's code.

* `test-defects.R` requires the prediction warning on the residuals of
  a mixture whose Wiener component's weight has run to 6e-11. Its
  parameters are flat there, and frmtmb 0.70.0 reads the prediction
  test in the unit-diagonal metric of the Hessian, where every row
  moves along that direction (a cosine of up to 8.3e-5); in the
  optimizer's units it was 3.8e-12, below the test's tolerance, and
  0.69.1 gave the residuals no warning. The estimates are still
  finite.

* rstan is a suggested package: `test-sampling.R` calls `rstan::`, and
  `R CMD check` warned that the import was not declared.
