# frmtmb.ode 0.3.0

* **`frm_lincmt()`, the analytic one- to three-compartment model.** A
  linear compartment model has a closed form, and this one writes it:
  the state is a superposition of the impulse responses of the doses
  that precede each observation, with a steady-state record summed as
  the geometric series it is rather than approached by a run-in. There
  is no integrator, so the fit needs neither RTMBode nor deSolve.
  The Phase 0 design of `dev/extension-gaps-plan.md` (100 subjects x 8
  samples, a depot and a central compartment, twice-daily dosing for
  seven days through `ii`/`addl` and one `ss` row) falls from 3717 s to
  **55 s** at the default `n_ss = Inf` and 119 s with the run-in
  written out. Taking each arm's slowest measurement over two passes
  against the other's fastest, that is **at least 67x** and at least
  31x (3716.6 / 55.0 = 67.6 and 3716.6 / 119.1 = 31.2). At a size
  where the same fit is cheap enough to repeat and
  `nlminb` takes identical iteration counts in every arm the factors
  are 95x and 47x, so the true figure is near 100x; the Phase 0
  design's arms run one after another rather than interleaved, so the
  conservative bound is what its own design proves. Every estimate
  matches to seven significant digits what `dev/scale-findings.md`
  records for that row. It also converges: `nlminb` returns code 0 with
  a maximum gradient of 4.7e-04, where the solver arm returned its
  false-convergence code 1 at 6.1e-03. The load-independent version of
  the same claim is the solve count: 33 `ADjoint` atomic solves per
  subject baked into the tape and replayed on every evaluation,
  against 0.

* The schedule grammar is `frm_ode()`'s and is not re-implemented:
  `time`, `value`, `state`, `method`, `duration`, `ii`, `addl`, `ss`,
  `group` and `event_scale` mean what `?frm_ode` says they mean,
  through the same validation. What superposition cannot carry is
  refused by name: `"replace"` and `"multiply"` rows, a non-zero
  `"reset"`, a dose into or an output from a peripheral compartment, an
  infusion into the depot, `tv`, and an infusion that spans a restart.

* **The numerics are the substance.** The textbook difference-of-
  exponentials form divides by differences of rate constants, and
  `ka == ke` is both the ordinary flip-flop of an oral model and the
  first evaluation of any fit that starts two log rates at the same
  value. Nothing in `frm_lincmt()` divides by such a difference.
  Against a 240-bit `Rmpfr` reference that never forms an eigenvalue,
  the worst relative error over 180 draws from a wide parameter box is
  7.1e-11, and over the sweeps through every coalescence (`ka` to `ke`,
  a two-compartment double root, a three-compartment double and triple
  root, each down to exact equality) it is 7.2e-15. Against
  `frm_ode()` at `atol = rtol = 1e-12` over 118 schedules the worst
  disagreement is 2.6e-12 of the trajectory's own scale, which is the
  solver's tolerance and not the closed form's.

* `frmtmb.sample::frm_sample()` samples a `frm_lincmt()` model, and
  `frm_lincmt()` is registered in the compatibility vocabulary with NO
  refusal row. That is a claim rather than an omission: the refusal on
  `frm_ode()` is about RTMBode calling deSolve unguarded, and the
  closed form calls neither. Measured: one chain of 400 iterations on a
  six-subject oral model returns in 7.7 s with the session intact.

# frmtmb.ode 0.2.0

* This package registers `frm_ode()` as refused for `frm_sample()`, so
  frmtmb.sample stops before taping instead of handing an empty
  `stanfit` back. The tmbstan and RTMBode boundary is still broken and
  it fails on a FIXED-EFFECT-ONLY model, which `dev/feature-gaps.md`
  had recorded as an open question.

* MEASURED, and it decides where this package goes next: one
  `frm(se = TRUE)` on 100 subjects by 8 samples with twice-daily
  dosing for 7 days takes 3948 seconds, with a 40 second tape build
  and a 14 second gradient. The truth is recovered, so it is the
  segmented sensitivity solve rather than the conditioning. An
  analytic one- to three-compartment path is now a prerequisite for
  population work rather than an option.

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
