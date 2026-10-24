# frmtmb.spatial (development version)

* New package. `car(M, gr = g, type = )`, brms's conditional
  autoregressive field in its four types (`"escar"`, `"esicar"`,
  `"icar"`, `"bym2"`), and `spde(fem, gr = node)`, a Matern field over a
  finite-element mesh, moved here from frmtmb. Both are registered
  through `frmtmb::frmtmb_register_covstruct()`, and the package reads
  no frmtmb internal.
  Fits are unchanged: on the moved tests' models the log-likelihood,
  every estimate and standard error, `fitted()` with and without an
  unseen location, `ranef()` conditional standard deviations,
  `confint_varcorr()` and `simulate()` are identical to frmtmb 0.73.0's.
  `frmtmb.sample::frm_sample()` gives the same draws from the same seed.
  `frm_compat()` gives `car` and `spde` the statuses frmtmb gave them.
