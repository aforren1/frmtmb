# Changelog

## frmtmb.coupling 0.4.0

- Requires frmtmb 0.60.0, for `frm_stop()` and the other condition
  helpers.

- Requires the frmtmb release that exports `frm_stop()`; the
  `frmtmb (>= 0.59.0)` floor must move to it.

- **BREAKING:** every error, warning and message that frmtmb.coupling
  raises is classed. An error has the class
  `c("frmtmb_coupling_error", "frmtmb_error", "error", "condition")`,
  and warnings and messages follow the same pattern, so
  `tryCatch(frmtmb_error = )` catches any refusal. The class vector no
  longer contains `simpleError`, `simpleWarning` or `simpleMessage`. See
  `?frmtmb::frmtmb-conditions`.

- **BREAKING:** a value that matches none of the choices of
  `frm_cross_spectrum(window =)` is refused with a
  `frmtmb_coupling_error` that names the argument, the value and the
  choices. The old text was `'arg' should be one of ...`. A partial
  value still matches.

- A refusal that frmtmb raises about
  [`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md),
  such as a missing addition term, is a `frmtmb_coupling_error`.

## frmtmb.coupling 0.3.2

- Requires frmtmb 0.59.0.
- Tests follow frmtmb 0.59.0’s `VarCorr()` structure, through
  `varcorr_matrices()`.

## frmtmb.coupling 0.3.1

- **A within-subject contrast needs `(1 | id:cond)`, and without it the
  interval is 2.7 times too narrow rather than in the wrong place.**
  Measured over 148 replicates of 40 subjects by 2 conditions by 60
  frequencies, with a true condition contrast of 0.5 on the logit scale:
  `coh ~ cond + s(freq, by = cond) + (1 | id) + (1 | id:cond)` covers
  141 of 148, 0.953 (0.906, 0.977), and recovers both components, while
  the same model without `(1 | id:cond)` covers 78 of 148, 0.527 (0.447,
  0.606). Paired on the same data its estimate sits 0.0036 from the
  correct model’s, so what fails is the width: on 63 of the 148
  replicates it missed where the correct model covered, and it never
  covered where the correct model missed.
  [`?cross_wishart`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md)
  and
  [`vignette("coherence")`](https://aforren1.github.io/frmtmb/frmtmb.coupling/articles/coherence.md)
  carry the table, and `dev/coh-findings.md` the construction, the seeds
  and a null arm in which the two models agree to a median width ratio
  of 1.000.

  Dropping `(1 | id)` instead, and keeping `(1 | id:cond)`, goes the
  other way: 1.9 times too wide, coverage 148 of 148. So the rule is
  neither more random effects nor fewer. A within-subject contrast needs
  both terms, the subject one for the intercept and the crossed one for
  the contrast.

- **A model with no random effects estimates the MARGINAL contrast, and
  estimates it correctly.** `coh ~ cond` and
  `coh ~ cond + s(freq, by = cond)` sit -0.0197 and -0.0133 from the
  subject-level contrast over the same 148 replicates. That is not a
  bias to be widened away: a logit contrast is not collapsible, and the
  factor `1 / sqrt(1 + 0.346 V)` with `V` the omitted variance predicts
  -0.0196 and -0.0135 from the simulator’s constants with nothing
  fitted.
  [`?cross_wishart`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md)
  and the vignette now say which of the two failures a reader is looking
  at, since the remedies differ. `dev/scale-findings.md` recorded the
  attenuation, at 0.016 across the ladder, before this study ran.

- The gated scale tier records `coh_cond_se`, the standard error of the
  condition contrast, beside the interval it produced. That is the
  column the five rungs of the ladder differ in.

## frmtmb.coupling 0.3.0

- [`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md)
  splits at runs of `NA` instead of refusing them, so an
  artifact-rejected record is one call rather than one per clean span. A
  record with two rejected spans gives the sum of the three clean
  pieces.

- `window = "hann"` beside the sine tapers, measured to cost no degrees
  of freedom (8.021 against 8.008 untapered at a nominal 8, 3000
  replicates) and identical to
  [`stats::spec.pgram()`](https://rdrr.io/r/stats/spec.pgram.html) to
  4.8e-16. Hann with `smooth` is refused, because it delivers 4.67 to
  5.66 where it claims 8.

- The help now states the cross-row correlation the ordinates carry,
  which was understated: 0.397 at the default and 0.753 at `tapers = 4`,
  which is one independent frequency in four. The degrees of freedom `n`
  reports were re-measured and are right.

- **[`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md)
  splits at NA runs instead of refusing them.** An artifact-rejected
  recording arrives with gaps in it, and it used to refuse the whole
  record for one missing sample. `NA` now marks a sample the pair does
  not have: the record is cut at it and segments are laid inside the
  clean spans, so no transform ever crosses a gap and an
  artifact-rejected recording is one call. A sample is usable only where
  BOTH signals have it. `NaN` counts as `NA`; an infinity is refused by
  name, because it is a value the arithmetic cannot use rather than a
  value the record is missing.

  The property that makes one call correct is an identity: a record with
  two rejected spans returns exactly the sum of the three clean pieces
  transformed separately, in `w11`, `w22`, `w12r`, `w12i` and in `n`,
  and the test file asserts it. The segment length comes from the usable
  sample count, each span supplies as many whole segments as fit in it
  in record order, and `n` is the count the record actually supplied
  rather than the count `segments` asked for. That count falls only when
  the record breaks into SEVERAL spans, each dropping its own remainder:
  a record with one clean span always yields the full `segments` however
  much was rejected from its ends, because `usable %/% segments` rounds
  the segment length down. Nothing warns about the shortfall, because
  where it happens it is the normal case; `n` is per-row data the
  density reads and every standard error downstream is formed from it.

- **`window = "hann"`, a data window beside the sine tapers**, for a
  spectrum steep enough that the untapered transform’s high frequencies
  are leakage from its low ones.

  **The degrees of freedom it costs are none, and that is measured
  rather than assumed.** Segments are disjoint, so a window inside each
  one leaves the number of independent complex draws alone. At a true
  coherence of zero, 3000 replicates on white noise at 512 samples, 8
  segments deliver an effective 8.021 (se 0.127) with the window against
  8.008 (se 0.130) without it. The scaling is `sqrt(N / sum(w^2))`,
  exactly `sqrt(8/3)` for the periodic raised cosine, so `w11 / n`
  estimates the same spectral density either way: 1.0013 (se 0.0022)
  untapered and 0.9963 (se 0.0032) windowed, on white noise of unit
  variance over 400 replicates. The arithmetic is pinned against
  [`stats::spec.pgram()`](https://rdrr.io/r/stats/spec.pgram.html) on
  the windowed segments, which agree to 4.8e-16 of the largest ordinate.

  **It cannot be combined with `smooth`, and that is refused rather than
  counted wrong**, the way `tapers` and `smooth` already are. A Hann
  window makes each ordinate a weighted sum of itself and its two
  neighbors, so adjacent bins are correlated and averaging them adds
  less than it claims: at nominal `n = 8` throughout, 4 segments over 2
  bins delivers 5.66, 2 segments over 4 bins delivers 4.98 and 1 segment
  over 8 bins delivers 4.68, against 8.19, 8.02 and 8.05 for the same
  three unwindowed. Declaring the shortfall instead was considered and
  rejected: there is no one number to declare, since the same window
  loses a different amount at each width, and frmtmb’s `whittle()`
  measured the same thing from the other side, where an honestly
  declared equivalent degrees of freedom does not rescue a smoothed
  spectrum from its raw-periodogram check at kernel widths of 7 and
  above. It cannot be combined with `tapers` either, because both are a
  taper on the same segment.

  **What `n` does not say** is now a section of the help page. `n`
  counts the draws behind one row and a window leaves that alone; what a
  window changes is the relation BETWEEN rows. See the next bullet,
  which carries the numbers for the window and for `tapers` together.

- **`tapers` above 1 correlates adjacent frequencies, and nothing said
  so.** This is not a change of behavior and `n` is not wrong. That was
  checked before this bullet was written, because a wrong `n` would be a
  correctness bug in a released package. Re-measured at a true coherence
  of zero over every retained ordinate, 600 replicates at 512 samples:

  | configuration                   | nominal | effective n       |
  |---------------------------------|---------|-------------------|
  | `segments = 8`                  | 8       | 7.968 (se 0.051)  |
  | `segments = 8, window = "hann"` | 8       | 7.964 (se 0.051)  |
  | `segments = 2, tapers = 4`      | 8       | 8.033 (se 0.026)  |
  | `segments = 1, tapers = 8`      | 8       | 8.006 (se 0.018)  |
  | `segments = 1, tapers = 16`     | 16      | 16.144 (se 0.039) |

  So the per-row degrees of freedom are right at every taper count this
  package documents, and a window costs none of them either. What was
  missing is the OTHER quantity, the relation between rows, which a user
  fitting `s(freq)` on a tapered frame has been assuming since 0.1.0.

  Measured on
  [`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md)’s
  own output, 400 replicates of white noise at 1024 samples, as the
  lag-1 correlation of `log(w11)` and as the variance inflation factor
  of a mean over the rows:

  | configuration                   | lag 1  | inflation |
  |---------------------------------|--------|-----------|
  | `segments = 8`                  | -0.010 | 1.00      |
  | `segments = 8, window = "hann"` | 0.397  | 1.79      |
  | `segments = 8, tapers = 2`      | 0.525  | 2.24      |
  | `segments = 4, tapers = 4`      | 0.745  | 4.07      |
  | `segments = 1, tapers = 8`      | 0.876  | 8.02      |

  So a `tapers = 4` frame carries about one independent frequency in
  four, against one in 1.8 for a Hann window, and a smooth in frequency
  fitted on either has fewer effective points than it has rows. The
  inflation factor is summed over twelve lags: at `tapers = 8` the
  correlation is still positive past lag 4, so stopping there reports
  6.45 for a row that is really 8.02, and every other row is the same
  either way. The help page now carries the table, and an earlier draft
  of this release quoted 0.301 for the window, which is the value at
  `segments = 1`, a call the function refuses.

## frmtmb.coupling 0.2.0

- The coherence complement now comes from frmtmb’s public accessor,
  `dpar_log1m()`, rather than from this package’s own copy of the same
  arithmetic. The copy existed because the accessor was internal to
  frmtmb: a family defined outside that package could see the reserved
  `.eta_<dpar>` entry and had no sanctioned way to read it. The two
  agree to the last bit on the tape, which is the path a fit runs on:
  measured on the full
  [`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md)
  log density over 961 values of the coherence linear predictor from 10
  to 700, the maximum difference is 0. Off the tape, where both fall
  back to the plain arithmetic, `1 / (1 - C)` is now `exp(-log(1 - C))`
  and the log density moves by at most 2.3e-13 on a value of -1693.4,
  one part in 7.4e15, at a coherence of 0.999. Deviance residuals are
  the only quantity that reads it.

  Nothing about the density changes: it still never subtracts from 1,
  and `log(1 - C)` is still exact to `eta = 709` against 36.74 for the
  naive form. This is a change of source, not of arithmetic.

- **This release needs a newer frmtmb.** The `Depends:` floor has to
  rise to the frmtmb release that exports `dpar_log1m()`. Building
  against an older frmtmb fails at install time, where the import cannot
  be resolved, rather than at run time.

- `RTMB` leaves `Imports:`. Its only use was the `logspace_add()` call
  the accessor replaces; the density is taped by frmtmb and calls no
  RTMB function of its own.

## frmtmb.coupling 0.1.0

First release. A complex Wishart family for the cross-spectrum of a
signal pair, with coherence and phase as distributional parameters.

- [`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md),
  the likelihood of a two-channel cross-spectral matrix. Its four dpars
  are the two channel powers, the magnitude squared coherence and the
  phase, so a random effect on coherence or a smooth in coherence over
  frequency is an ordinary formula. The logit link on coherence IS the
  positive definiteness constraint: no linear predictor can produce an
  invalid spectral matrix.

- [`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md),
  which turns two signals into the rows the family reads. It drops
  frequency zero and the Nyquist frequency, whose coefficients are real
  rather than complex, and it counts the degrees of freedom from
  segments, sine tapers and frequency smoothing together. Overlapping
  segments are refused, because overlap raises the nominal count without
  raising the independent one.

- [`frm_coherence()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_coherence.md)
  and
  [`frm_phase()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_coherence.md),
  which read the fit through `predict(se.fit = TRUE)` and return
  intervals. The coherence interval is formed on the logit scale, so it
  cannot leave the unit interval.

- [`frm_cross_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_simulate.md),
  which draws whole Hermitian matrices.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) refuses on this
  family for a stated reason: a draw is a matrix and the response slot
  holds one of its four numbers.

- At two channels the complex Wishart density is real arithmetic in
  closed form: `det S = S11 S22 (1 - C)` and the trace term expands to
  four real products. RTMB 1.9’s complex automatic differentiation was
  used to derive and check the general form, and it tapes correctly at
  four channels including under a Laplace approximation, but at two it
  measured 5.8 times slower than the closed form for the same answer to
  six decimals. The only thing this package takes from RTMB is
  `logspace_add()`, for the next item.

- The coherence complement is computed on the LOG scale, from the linear
  predictor rather than by subtracting from 1. `plogis(eta)` is exactly
  1 in double precision above about `eta = 36.74`, so `1 - C` there is
  exactly 0 and a density written over it returns `NaN` for both value
  and gradient; below that it is merely wrong, by 1e-3 relative at
  `eta = 30`. That region is reachable: a random effect on `coh` with
  one group near a coherence of 1 drives the estimate past `eta = 34`,
  and two signals differing by 1e-5 of noise put an intercept-only fit
  at `eta = 23`. Core stores each dpar’s linear predictor beside it, so
  `log(1 - C)` is `-logspace_add(0, eta)`, exact to `eta = 709`.
  Measured: the log density is finite and correct at every `eta` from 10
  to 700 where it was `NaN` above 36.74.

- [`frm_coherence()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_coherence.md)
  and
  [`frm_phase()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_coherence.md)
  refuse a standard error that is zero or not finite instead of
  returning it. A 95 percent interval whose lower, estimate and upper
  are the same number is a wrong answer rather than a narrow one, and
  that is what the two-signals-differing- by-1e-5 case used to produce.

- There is NO automatic warning for a power model too simple for the
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

- More than two channels is refused as a SCOPE decision of this package.
  Core is not the obstacle: `frm()` carries a matrix-valued response and
  a custom family reading `y[, 1]` and `y[, 2]` fits today. An earlier
  draft of this package said otherwise and was wrong.

- The measurements behind all of this are in `dev/xspec-findings.md`:
  the naive coherence bias table, the effective degrees of freedom of
  each route to `n`, and the interval coverage that argues for putting a
  random effect on every dpar rather than only on the coherence.
