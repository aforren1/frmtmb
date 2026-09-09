# frmtmb.ode (development version)

* **The package registers itself with the compatibility matrix.**
  `frm_ode()` is a feature of kind `special` in
  [frmtmb::frm_compat()]'s vocabulary, and the pair
  `frm_ode() x frm_sample` is registered as `refused`.
  `frmtmb.sample::frm_sample()` and `frmtmb.sample::as_tmbstan()` both
  read that row and stop before taping, so sampling an ODE fit now
  fails in the time a registry read takes and says why, instead of
  aborting a chain at warmup iteration 1 and leaving rstan's autodiff
  arena unusable for the rest of the session. `as_tmbstan()` was the
  worse of the two doors: it returned an empty `stanfit` with no error
  at all.
  The defect is upstream, is reproduced without frmtmb in
  `dev/upstream/`, and still reproduces against RTMBode at commit
  5242257: `RTMBode::ode()` calls `deSolve::ode()` with no guard, so a
  failed solve escapes as an R error where a sampler needs `NaN`.
  `?frm_ode` gains a "Sampling an ODE fit" section with the evidence
  and the two ways forward. The row goes away when a fixed RTMBode is
  released.

# frmtmb.ode 0.1.1

Requires frmtmb 0.55.0, for the hazard-container lint that now
runs in this package's own check.

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
