# Changelog

## frmtmb.ode 0.1.0

First release, extracted from frmtmb 0.46.0.

- [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  writes ordinary differential equation dynamics inside a nonlinear
  frmtmb formula, solved by RTMBode with adjoint differentiation through
  the dynamics parameters and initial states.
  [`frm_ode_failures()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode_failures.md)
  reads the solver failure log of a fit.
- The package registers its formula check with the core at load through
  [`frmtmb::frmtmb_register_frame_check()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_frame_check.html),
  so a misuse of
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  is refused with the same message it always was.
- RTMBode is not on CRAN; this package carries the
  `Additional_repositories` pointer to kaskr’s r-universe, so core
  frmtmb’s dependencies stay CRAN-only.
