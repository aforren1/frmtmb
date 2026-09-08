# frmtmb.coupling 0.2.0

* The coherence complement now comes from frmtmb's public accessor,
  `dpar_log1m()`, rather than from this package's own copy of the same
  arithmetic. The copy existed because the accessor was
  internal to frmtmb: a family defined outside that package could see
  the reserved `.eta_<dpar>` entry and had no sanctioned way to read
  it. The two agree to the last bit on the tape, which is the path a
  fit runs on: measured on the full `cross_wishart()` log density over
  961 values of the coherence linear predictor from 10 to 700, the
  maximum difference is 0. Off the tape, where both fall back to the
  plain arithmetic, `1 / (1 - C)` is now `exp(-log(1 - C))` and the
  log density moves by at most 2.3e-13 on a value of -1693.4, one part
  in 7.4e15, at a coherence of 0.999. Deviance residuals are the only
  quantity that reads it.

  Nothing about the density changes: it still never subtracts from 1,
  and `log(1 - C)` is still exact to `eta = 709` against 36.74 for the
  naive form. This is a change of source, not of arithmetic.

* **This release needs a newer frmtmb.** The `Depends:` floor has to
  rise to the frmtmb release that exports `dpar_log1m()`.
  Building against an older frmtmb fails at install time, where the
  import cannot be resolved, rather than at run time.

* `RTMB` leaves `Imports:`. Its only use was the `logspace_add()` call
  the accessor replaces; the density is taped by frmtmb and calls no
  RTMB function of its own.

# frmtmb.coupling 0.1.0

First release. A complex Wishart family for the cross-spectrum of a
signal pair, with coherence and phase as distributional parameters.

* `cross_wishart()`, the likelihood of a two-channel cross-spectral
  matrix. Its four dpars are the two channel powers, the magnitude
  squared coherence and the phase, so a random effect on coherence or a
  smooth in coherence over frequency is an ordinary formula. The logit
  link on coherence IS the positive definiteness constraint: no linear
  predictor can produce an invalid spectral matrix.

* `frm_cross_spectrum()`, which turns two signals into the rows the
  family reads. It drops frequency zero and the Nyquist frequency,
  whose coefficients are real rather than complex, and it counts the
  degrees of freedom from segments, sine tapers and frequency smoothing
  together. Overlapping segments are refused, because overlap raises
  the nominal count without raising the independent one.

* `frm_coherence()` and `frm_phase()`, which read the fit through
  `predict(se.fit = TRUE)` and return intervals. The coherence interval
  is formed on the logit scale, so it cannot leave the unit interval.

* `frm_cross_simulate()`, which draws whole Hermitian matrices.
  `simulate()` refuses on this family for a stated reason: a draw is a
  matrix and the response slot holds one of its four numbers.

* At two channels the complex Wishart density is real arithmetic in
  closed form: `det S = S11 S22 (1 - C)` and the trace term expands to
  four real products. RTMB 1.9's complex automatic differentiation was
  used to derive and check the general form, and it tapes correctly at
  four channels including under a Laplace approximation, but at two it
  measured 5.8 times slower than the closed form for the same answer to
  six decimals. The only thing this package takes from RTMB is
  `logspace_add()`, for the next item.

* The coherence complement is computed on the LOG scale, from the linear
  predictor rather than by subtracting from 1. `plogis(eta)` is exactly
  1 in double precision above about `eta = 36.74`, so `1 - C` there is
  exactly 0 and a density written over it returns `NaN` for both value
  and gradient; below that it is merely wrong, by 1e-3 relative at
  `eta = 30`. That region is reachable: a random effect on `coh` with
  one group near a coherence of 1 drives the estimate past `eta = 34`,
  and two signals differing by 1e-5 of noise put an intercept-only fit
  at `eta = 23`. Core stores each dpar's linear predictor beside it, so
  `log(1 - C)` is `-logspace_add(0, eta)`, exact to `eta = 709`.
  Measured: the log density is finite and correct at every `eta` from 10
  to 700 where it was `NaN` above 36.74.

* `frm_coherence()` and `frm_phase()` refuse a standard error that is
  zero or not finite instead of returning it. A 95 percent interval
  whose lower, estimate and upper are the same number is a wrong answer
  rather than a narrow one, and that is what the two-signals-differing-
  by-1e-5 case used to produce.

* There is NO automatic warning for a power model too simple for the
  data, and the reason is worth stating. A draft carried a `fit_check`
  hook on the spread of the power residuals; that statistic is not
  monotone in the coherence displacement it was meant to catch. It went
  silent at a 2.7-fold power spread where the reported coherence was
  already nearly double the truth, and it fired at a ratio of 20.5 on a
  correctly specified model whose coherences were right. It was
  withdrawn rather than retuned. Compare AIC with and without a
  frequency term on the powers instead: measured across power spreads
  from 1.6-fold to 385-fold, AIC rejects the flat model by 80 or more in
  every case, including the two where the hook was silent.

* More than two channels is refused as a SCOPE decision of this package.
  Core is not the obstacle: `frm()` carries a matrix-valued response and
  a custom family reading `y[, 1]` and `y[, 2]` fits today. An earlier
  draft of this package said otherwise and was wrong.

* The measurements behind all of this are in `dev/xspec-findings.md`:
  the naive coherence bias table, the effective degrees of freedom of
  each route to `n`, and the interval coverage that argues for putting
  a random effect on every dpar rather than only on the coherence.
