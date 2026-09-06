# Frequency-domain models in frmtmb: what works, what is missing

Date: 2026-09-06. frmtmb 0.53.0, RTMB 1.9, TMB 1.9.25, R 4.6.1,
Windows 11. **No package code was changed to produce any number below.**
Everything measured here runs through the released grammar.

## Verdict

Most of a frequency-domain offering already exists, because of one
identity: the Whittle likelihood is an exponential GLM on the
periodogram with a log link. `exponential()` is mean-parameterized
(verified: `lpdf(y = 2, mu = 3)` equals `dexp(2, rate = 1/3, log =
TRUE)`), and a periodogram ordinate of a stationary Gaussian series is
exponential about the spectral density. So a spectral model is an
ordinary `frm()` fit whose response is the periodogram and whose linear
predictor is the log spectrum.

What is missing is not the likelihood. It is the data preparation, the
post-processing vocabulary, and the multivariate case.

## What was measured

| model | data | result |
| --- | --- | --- |
| AR(1) by Whittle, nonlinear log spectrum | n = 2048, 1023 ordinates | phi 0.702534 against `arima()` 0.701628, truth 0.7 |
| the same, standard error on phi | | 0.01577 against `arima()` 0.01572 |
| nonparametric spectrum, `I ~ s(w, k = 20)` | same series | max relative error 0.113 against the truth, against 0.0115 for the correct parametric model |
| hierarchical AR(1), `z ~ 1 + (1 \| id)` | 12 subjects, 3060 ordinates | 1.05 s; between-subject SD 0.317 (empirical 0.336, truth 0.35); per-subject correlation 0.985 |
| aperiodic plus peak, three random effects | 40 subjects, 7080 rows | 2.2 s; group effect on the exponent 0.262 se 0.083 (truth 0.35) |
| the same across a montage | 30 subjects, 20 channels, 53,400 rows | 7.95 s; subject exponent SD 0.241, channel 0.177 (truth 0.25, 0.15) |
| individual alpha frequency | 40 subjects | model r = 0.81, worst error 1.95 Hz; band argmax r = 0.52, worst error 2.99 Hz |

The last row is the argument for doing this at all. The peak frequency
of each subject, estimated as a random effect with shrinkage, beats the
argmax of that subject's periodogram inside a fixed band, which is what
analysis pipelines do today.

`R/covstruct.R:1608` already evaluates a squared-exponential spectral
density at multi-index frequencies for the Hilbert-space Gaussian
process, so "the parameters enter through a spectral density" is
machinery this package already carries and tests.

## What a real offering needs

Ordered by how much each one costs.

1. **A periodogram constructor.** `frm_periodogram(y, ...)`: remove the
   mean, optionally taper, return the Fourier frequencies with their
   ordinates, and **drop frequency zero and Nyquist**, which are
   chi-square with one degree of freedom rather than two and so break
   the exponential likelihood. Every one of these is a correctness
   detail a user gets wrong once and never notices. Should carry the
   sampling rate so frequencies come back in Hz.

2. **A `whittle()` family wrapper.** Thin: `exponential(link = "log")`
   for a raw periodogram, `Gamma(link = "log")` with the shape fixed at
   the taper count for a multitaper or smoothed one. Fixing the shape
   needs `map`, which the fit already supports. The wrapper exists to
   name the assumption and to refuse a response that cannot be a
   periodogram.

3. **Honest post-processing.** `predict()`, `residuals()` and
   `simulate()` on such a fit all speak about the periodogram, not the
   series. `simulate()` returning spectral ordinates where a user
   expects a time series is the trap. Either document it plainly in the
   family, or refuse `simulate()` by name and provide an inverse
   transform as a separate verb.

4. **The small-sample bias.** Whittle is biased when the series is
   short or close to a unit root, and an EEG epoch is often one or two
   seconds. The debiased and tapered variants are the standard remedy
   and both change the likelihood, not just the data preparation. A
   vignette that recommends Whittle without measuring this on realistic
   epoch lengths would be misleading.

5. **Complex cross-spectra.** The one piece that is genuinely new. The
   cross-spectrum between two series is complex, and the multivariate
   periodogram at each frequency is complex Wishart, so coherence and
   phase need a new family rather than a new formula. **RTMB 1.9 already
   has complex automatic differentiation**: an `adcomplex` class with
   `fft`, `Re`, `Im`, `Mod`, `Arg`, `Conj` and `solve_complex`. Nothing
   in frmtmb uses any of it. That also unlocks a spectral density with
   no closed form, computed by transforming a parametric autocovariance
   on the tape, and circulant embedding for stationary fields.

## Where it should live

Items 1 to 4 are a vignette plus two helper functions, not a package.
They belong in core, beside the reinforcement-learning vignette, because
they add no families and no estimation machinery.

Item 5 is an extension. It needs a complex Wishart likelihood, a
multivariate response whose rows are Hermitian matrices, and its own
post-processing vocabulary for coherence, phase and directed measures.
Name it for what it does with signals rather than for the transform.

## Applications this serves

The demand is in psychology and neuroscience, where the current practice
is a two-stage pipeline: fit each spectrum alone, extract numbers, then
test the numbers. Replacing that with one hierarchical fit is the
selling point, and the aperiodic exponent and the individual alpha
frequency are the two quantities people most want.

Items 1 to 4 serve: the aperiodic exponent as a dependent variable in
aging, attention deficit disorder, schizophrenia and anesthesia depth;
the theta to beta ratio, whose group differences are substantially
aperiodic rather than periodic; individual alpha frequency for
individualizing bands and stimulation frequencies; developmental
trajectories through a smooth in a nonlinear parameter; single-trial
spectra with trial-level covariates; and, outside the brain, heart rate
variability, pupillometry, circadian actigraphy, postural sway and
tremor frequency.

Item 5 serves everything that involves two signals: corticomuscular
coherence, cortico-cortical connectivity, neural tracking of the speech
envelope, hyperscanning, and directed measures built from a parametric
cross-spectrum. Coherence is biased upward when segments are few, and
group comparisons of it are awkward precisely because the sampling
distribution is not modeled.

**Out of scope, and worth saying so in the vignette.** Phase-amplitude
coupling is a higher-order property that needs the bispectrum, not this
framework. Time-frequency analysis of a non-stationary epoch needs
segmentation or a different formulation.

## Blockers found while measuring

Four usability defects surfaced while writing the models above. They are
not specific to spectra, they are what a nonlinear hierarchical model
runs into, and they are recorded in the round's API lane rather than
here.

1. A nonlinear parameter may not be called `mu` or `b`, and neither
   error message says the name is reserved.
2. `fixef()` names use a dot separator and `vcov()` names an underscore.
3. `ranef()` returns one block per nonlinear parameter, all named after
   the grouping factor, distinguishable only by an attribute.
4. A Gaussian peak term with no starting value is flat over the whole
   frequency range, so the fit returns undefined standard errors under a
   warning that says the model is probably overparameterized.

Item 4 matters most for this work: every peak model needs a starting
value or a centering offset, and the vignette must say so.
