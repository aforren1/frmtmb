# frmtmb.ode 0.1.0

First release, extracted from frmtmb 0.46.0.

* `frm_ode()` writes ordinary differential equation dynamics inside a
  nonlinear frmtmb formula, solved by RTMBode with adjoint
  differentiation through the dynamics parameters and initial states.
  `frm_ode_failures()` reads the solver failure log of a fit.
* The package registers its formula check with the core at load
  through `frmtmb::frmtmb_register_frame_check()`, so a misuse of
  `frm_ode()` is refused with the same message it always was.
* RTMBode is not on CRAN; this package carries the
  `Additional_repositories` pointer to kaskr's r-universe, so core
  frmtmb's dependencies stay CRAN-only.
