Needs frmtmb 0.77.0. That release reads brms's `(flat)` as a prior
density and gives `build_objective()` its `prior_only` field. Its
`sigma` is the innovation sd under the covariance form of `ar()`,
`ma()` and `arma()`.
