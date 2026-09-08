# Whittle likelihood for a periodogram response

Names the assumption behind a spectral fit. A periodogram ordinate of a
stationary Gaussian series is exponential about the spectral density, so
with a log link the linear predictor IS the log spectrum and nothing
else in the fit changes. `whittle()` is that exponential family, with
the periodogram-specific refusals attached.

## Usage

``` r
whittle(tapers = 1, link = "log")
```

## Arguments

- tapers:

  Number of independent periodograms averaged into each ordinate. `1`
  (the default) is the raw periodogram and gives the exponential
  likelihood; `k > 1` gives Gamma with the shape fixed at `k`.

- link:

  Link for the mean. The default `"log"` is what makes the linear
  predictor the log spectrum, and is very nearly always what is wanted.

## Value

A `frmtmb_family`.

## Details

`tapers = k` for `k > 1` is the Gamma likelihood with the shape held at
`k`, which is the distribution of an average of `k` independent
ordinates: a Welch or Bartlett estimate from
`frm_periodogram(segments = k)`, or a `k`-taper multitaper estimate. The
shape is data preparation, not a parameter, so it is fixed rather than
estimated - through the same mechanism as `bf(y ~ x, shape = k)`, and an
explicit `shape` in
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) still wins.

## Post-processing speaks about the periodogram

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
fitted spectral density at each row's frequency,
[`predict()`](https://rdrr.io/r/stats/predict.html) the same for new
frequencies, and [`residuals()`](https://rdrr.io/r/stats/residuals.html)
compares an ordinate with that density. None of them is about the
series, and none of them can be: the phases are not in a periodogram.

[`simulate()`](https://rdrr.io/r/stats/simulate.html) is **refused by
name** rather than documented, and so are the four entry points built on
it:
[`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md),
[`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md),
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md),
and `conditional_effects(method = "predict")`. Draws from this family
are ordinates, a vector of them under the response's name reads as a
simulated series, and that is a mistake a user makes once and never
notices. Each refusal names the two things to do instead: draw ordinates
in the open, with `rexp(nobs(fit), 1 / fitted(fit))` (or
`rgamma(nobs(fit), shape = k,` `scale = fitted(fit) / k)` at
`tapers = k`), which is also the one-line replacement for a DHARMa-style
check; or draw a series whose spectrum is the fitted one with
[`frm_series_draw()`](https://aforren1.github.io/frmtmb/reference/frm_series_draw.md).

## What the family refuses, and what it cannot see

Two things are checked, and both look only at the response.

A response that is not positive is refused outright: a periodogram
ordinate is a squared modulus, so a negative one means log power or
power in dB, which this likelihood is not about.

A response too SMOOTH to have the declared shape is also refused.
`var(diff(log(y), lag = 3))` is `2 trigamma(tapers)` for ordinates of
that shape whatever the spectrum is, and only the spectrum's own
variation across the step adds to it, so a value far below that is
evidence the ordinates were averaged more than `tapers` says, or that
they are spectral leakage rather than signal. The threshold is a
calibrated fraction of the model value that shrinks as the number of
ordinates falls.

The step is three ordinates rather than one because of the taper. A Hann
window's transform is three bins wide, so it correlates NEIGHBORING
ordinates (about 0.31 on the log scale) and leaves ordinates three apart
uncorrelated (-0.003). At lag one a legitimate Hann-tapered response was
refused about 2% of the time; at lag three the measured rate is 0 in
2000 replicates at each of nine cells, and the same holds for Hamming
and Blackman windows, which this package does not apply and which were
refused 1.0% and 18.1% of the time at lag one. The wider step costs
leakage detection at the shortest usable length, 1.3 points at exponent
3 and 2.3 at exponent 2.5 with 127 ordinates, and under one point at 255
ordinates and above (2000 replicates a cell).

**What it cannot see.** It needs the rows in frequency order, which is
what
[`frm_periodogram()`](https://aforren1.github.io/frmtmb/reference/frm_periodogram.md)
returns: a periodogram sorted by its own power is refused every time,
wrongly. It does not run at all below 24 ordinates. It separates
`tapers = 2` from a raw periodogram only above about 200 ordinates,
because Gamma(2) and Gamma(1) are close. It reads an estimate smoothed
ACROSS FREQUENCY (a Daniell window, a multitaper) poorly at both ends.
Declared raw, such an estimate is caught every time above about 100
ordinates and often missed below it: a three-bin Daniell smooth is
caught 18% of the time at 32 ordinates and 58% at 64. Declared honestly,
by giving `tapers` the equivalent degrees of freedom, it passes up to
about five bins of smoothing and is then refused anyway, 37% to 89% of
the time at seven bins and 92% or more at eleven, because the check
assumes ordinates independent across frequency and a wide smooth is not.
Fit a frequency-smoothed estimate with care, or fit the raw or
segment-averaged periodogram it came from. And it is a backstop, not a
test to rely on: at 127 ordinates (a one-second epoch at 256 Hz) an
untapered exponent-3 power law is caught about 92% of the time and an
exponent-2 one much less often. Choose the taper from the shape of the
spectrum, not from whether an error appeared.
[`vignette("spectral")`](https://aforren1.github.io/frmtmb/articles/spectral.md)
has the measured rates.

## What it is not

The likelihood treats ordinates as independent, which is exact only in
the limit. It is biased for a short series or a spectrum near a unit
root;
[`vignette("spectral")`](https://aforren1.github.io/frmtmb/articles/spectral.md)
measures how much and says where it is safe. Tapering reduces that bias
and correlates the rows.

## See also

[`frm_periodogram()`](https://aforren1.github.io/frmtmb/reference/frm_periodogram.md),
[`frm_series_draw()`](https://aforren1.github.io/frmtmb/reference/frm_series_draw.md),
[`vignette("spectral")`](https://aforren1.github.io/frmtmb/articles/spectral.md).

## Examples

``` r
set.seed(1)
y <- as.numeric(arima.sim(list(ar = 0.6), 512))
pg <- frm_periodogram(y)

# a nonparametric log spectrum: the smoothing parameter is estimated
# by the same Laplace marginal likelihood as any other smooth
fit <- frm(bf(pgram ~ s(freq, k = 8)), family = whittle(), data = pg)
fixef(fit)
#> $mu
#> (Intercept) s(freq).fx1 
#>  0.02596863 -0.82825927 
#> 

# an averaged periodogram needs the shape it was averaged with
pg4 <- frm_periodogram(y, segments = 4)
frm(bf(pgram ~ log(freq)), family = whittle(tapers = 4), data = pg4)
#> frmtmb fit: pgram ~ log(freq) 
#> Family: whittle   Method: ML 
#>  Links: mu = log; shape = log
#> 
#> logLik: -43.121  AIC: 90.242  nobs: 63 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)   log(freq) 
#>     -1.5726     -0.9801 
#>  shape:
#> (Intercept) 
#>       1.386 
```
