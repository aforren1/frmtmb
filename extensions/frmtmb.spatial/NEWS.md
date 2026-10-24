# frmtmb.spatial 0.1.0

The first release. Load it with `library(frmtmb.spatial)`, which
attaches frmtmb as well. Needs frmtmb 0.74.0, which exports
`frmtmb_register_covstruct()`.

## New features

* New package. Two fields moved here from frmtmb.
  `car(M, gr = g, type = )` is brms's conditional autoregressive field
  in its four types (`"escar"`, `"esicar"`, `"icar"`, `"bym2"`).
  `spde(fem, gr = node)` is a Matern field over a finite-element mesh.
  Both are registered through `frmtmb::frmtmb_register_covstruct()`, and
  the package reads no frmtmb internal. Fits are unchanged. On the moved
  tests' models, these are identical to frmtmb 0.73.0's: the
  log-likelihood, every estimate and standard error, `fitted()` with and
  without an unseen location, `ranef()` conditional standard deviations,
  `confint_varcorr()` and `simulate()`. `frmtmb.sample::frm_sample()`
  gives the same draws from the same seed. `frm_compat()` gives `car`
  and `spde` the statuses frmtmb gave them.
