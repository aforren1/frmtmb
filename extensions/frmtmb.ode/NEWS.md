# frmtmb.ode 0.4.0

* **`frm_ode()`'s steady-state run-in sums the tail it used to drop.**
  An `ss` row is reached by repeating the dosing cycle `n_ss` times, and
  truncating there leaves the state short of the limit by
  `exp(-n_ss * lambda_z * ii)`, where `lambda_z` is the slowest
  disposition eigenvalue. That shortfall grows with the terminal
  half-life measured in dosing intervals, no check could run during a
  fit, and it moved estimates: the gradient with respect to `log(k21)`
  on a 107 hour two-compartment oral model had the WRONG SIGN, +10.37
  where the exact steady state gives -1.30, which is how a fit at
  `n_ss = 20` moved `k21` by a factor of 3.2. At population scale, on a
  30-subject dataset simulated from the exact steady state, truncating
  the run-in moves the objective by 60.66 units and turns
  `d/dlog(k21)` from +2.915 into -294.63.

  The run-in already computes the last cycle-start states, and their
  successive differences give the per-cycle contraction, so the sum
  that is left is `d * r / (1 - r)`. Adding it is arithmetic on tape
  variables, so unlike a convergence test it runs DURING a fit, at
  whatever parameters the fit has reached, and it costs no extra solve:
  the count is 21 solves per group at `n_ss = 20` either way. Measured
  worst over one dosing interval against the exact limit, at the
  shipped default and `atol = rtol = 1e-8`, on two-compartment oral
  models: 3.9e-03 to 9.9e-09 at a 23 hour half-life dosed every 8
  hours, 1.2e-02 to 2.6e-09 at 107 hours dosed daily, and 2.1e-01 to
  1.4e-08 at 670 hours dosed daily. On 22 ordinary population
  pharmacokinetic schedules the error is now 6.6e-10 to 6.1e-08
  everywhere, which is the integrator's own tolerance rather than the
  run-in's.

  **This changes the numbers a model with `ss` rows returns.**
  `ss_extrapolate = FALSE` restores the truncated run-in, bit for bit
  against 0.3.0 on 26 schedules, `identical()` on 26 of 26.

  **The default can lose, and where it does is bounded.** The ratio the
  run-in reads is a weighted mean of the cycle map's modes, and when
  two of the weights have opposite signs it overshoots. Opposite signs
  are the normal arrangement in an oral model, so the hazard is `ka`
  near `lambda_z`: **flip-flop kinetics**, which extended-release and
  depot formulations are written to produce. Over 338 two-compartment
  oral cycle maps the correction loses on 8, worst by a factor of 2.27,
  and every one of those has `lambda_z * ii` = 0.05. The bound is what
  makes the default defensible: the smallest error truncation leaves on
  a losing case is 0.40 and the largest it leaves on a winning one is
  0.96, so the correction never loses where truncation was usable, and
  the warning fires on every losing case in BOTH arms. Through
  `frm_ode()` on the worst such model, a 333 hour terminal half-life
  dosed daily: 5.2e-01 truncated against 8.0e-01 extrapolated. On a
  Michaelis-Menten system deep in its saturated regime it improved a
  1.1e-01 shortfall to 7.2e-03 and no further.

  The correction is applied in full while the measured ratio is between
  0.05 and 0.9875, gated to nothing below the first and stood down to
  nothing as it reaches 1. The gate is the degree-7 smootherstep, whose
  first three derivatives vanish at both ends; the stand-down is that
  same shape divided by its argument, degree 6, whose slope at 1 is -1,
  the slope that joins `r / (1 - r)` below the cap. Measured rather
  than asserted, **the objective has no corner in the parameters** at
  any of the four junctions, which matters because a fit whose `k21`
  runs to zero walks through one of them. It costs about 260 nodes on
  the tape, a fixed cost independent of `n_ss`, which is +124 percent
  of the outer tape at `n_ss` = 20 on a three-state model; the
  truncated arm's tape is 0.3.0's node for node.

* **The run-in warning now reports a distance to the limit.** It used
  to report the movement between the last two cycles, which understates
  the distance by about `1 / (lambda_z * ii)` and so understated it
  most exactly where the error was largest: measured 0.79x, 3.1x, 6.6x
  and 10.8x as the half-life grows. What it prints now is within 0.65x
  to 0.84x of the true error on the same four rows. It also says so by
  name when a state is not contracting between cycles at all, which
  means no `n_ss` reaches a steady state for it, and that is reported
  on its own evidence rather than through the distance. On 28
  schedules classified against the exact limit there is **no false
  alarm and no miss in either arm**, including three constructions
  where each of the correction's own guards fires, two where `n_ss` is
  too small for a tail to be read at all, and three where it is in the
  hundreds and inside the stand-down band.

  **Treat the number as a detector and not as a measurement.** It is
  built out of the same geometric model the correction is, so where
  that model is poor it is poor with it, and it is a lower bound
  rather than an estimate: measured across 64 runs the ratio of what
  it prints to the true error spans 0.0104 to 9.17e+09. It is good for
  deciding whether to look, not for how much to trust the third digit.

  A very long run-in has a limit that raising `n_ss` cannot pass. The
  cycles are chained solves, so the integrator's own error accumulates
  over all `n_ss + 1` of them: at `atol = rtol = 1e-8` and
  `n_ss` = 1000 it contributes about 2e-05, larger than the run-in
  shortfall the extra cycles were bought to remove. Past a few hundred
  cycles, tighten the tolerances too.

* `ss_tol` is documented as a distance to the steady state rather than
  as a relative change between cycles, which is what it now compares
  against.

* `frm_ode()` gains `ss_extrapolate`, default `TRUE`. It sits between
  `ss_tol` and `method`, so a call that matched `method`, `atol`,
  `rtol`, `on_error` or `penalty` POSITIONALLY now matches one argument
  earlier. Name them.

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
