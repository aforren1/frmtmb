# Periodogram for a Whittle-likelihood fit

Turns one or more time series into the data frame a spectral model is
fitted to: one row per retained Fourier frequency, with the ordinate
that
[`whittle()`](https://aforren1.github.io/frmtmb/reference/whittle.md)
treats as the response.

## Usage

``` r
frm_periodogram(
  x,
  fs = NULL,
  taper = c("none", "hann", "split_cosine"),
  p = 0.1,
  detrend = c("mean", "none", "linear"),
  segments = 1L,
  group = NULL
)
```

## Arguments

- x:

  A numeric vector (one series), a `ts`, or a matrix or data frame whose
  **columns** are series of equal length.

- fs:

  Sampling rate in Hz. Defaults to 1 (cycles per sample), or to
  `stats::frequency(x)` when `x` is a `ts`.

- taper:

  `"none"`, `"hann"`, or `"split_cosine"` (the split cosine bell of
  [`stats::spec.taper()`](https://rdrr.io/r/stats/spec.taper.html),
  applied to a proportion `p` of each end).

- p:

  Proportion tapered at each end for `"split_cosine"`.

- detrend:

  `"mean"` (the default), `"none"`, or `"linear"`. `"mean"` and `"none"`
  return identical ordinates, because the retained frequencies are
  exactly those where a constant contributes nothing; `"linear"` does
  change them.

- segments:

  Number of non-overlapping blocks to average. See Segment averaging.

- group:

  Optional grouping vector as long as `x` (a vector only), splitting it
  into series that need not have equal lengths.

## Value

A data frame with columns `freq` (Hz), `pgram` (the ordinate), and, when
the input names more than one series, `series` (a factor). Rows are
ordered by series and then by frequency.

## Details

The zero frequency and the Nyquist frequency are **dropped**. Their
ordinates are chi-square with one degree of freedom rather than two, so
they are not exponential and do not belong in the likelihood; a fit that
keeps them is wrong in a way nothing downstream reports. A series of
length `m` therefore returns `floor((m - 1) / 2)` rows.

## Scaling

The ordinate is
`|sum_t h_t (x_t - trend) exp(-2i pi f t / fs)|^2 / (fs sum_t h_t^2)`
and `freq` is in cycles per unit time (Hz when `fs` is in Hz), which is
the scaling and the frequency axis of
[`stats::spec.pgram()`](https://rdrr.io/r/stats/spec.pgram.html). So the
ordinates of white noise average to `var(x) / fs`, and a spectrum
written for the model must be a density per Hz. With `fs = 1` the
frequency axis runs to just under `0.5`, in cycles per sample.

## What a taper does to the likelihood

A taper trades bias for dependence. Multiplying the series by a window
before the transform suppresses the leakage that makes a steep
spectrum's high frequencies too large, which is the main small-sample
bias of the Whittle likelihood; in exchange, neighboring ordinates are
no longer close to independent, so the likelihood is a working one and
standard errors from it are optimistic. The exponential marginal
survives: dividing by `sum(h^2)` keeps `E I(f)` at the spectral density
smoothed by the window, so the mean model is unchanged and only the
correlation between rows is not.

## Segment averaging

`segments = k` splits the series into `k` non-overlapping blocks of
`floor(n / k)` points each (a remainder at the end is dropped),
periodograms each block, and averages them. The average of `k`
independent exponential ordinates is Gamma with shape `k` and the same
mean, so the matching family is `whittle(tapers = k)`: the two arguments
are one decision written twice, and this is the seam to get right. The
price is resolution, since the grid belongs to the block length
`floor(n / k)` rather than to `n`.

## See also

[`whittle()`](https://aforren1.github.io/frmtmb/reference/whittle.md)
for the family, and
[`vignette("spectral")`](https://aforren1.github.io/frmtmb/articles/spectral.md).

## Examples

``` r
set.seed(1)
y <- as.numeric(arima.sim(list(ar = 0.7), 512))
pg <- frm_periodogram(y, fs = 256)
str(pg)
#> 'data.frame':    255 obs. of  2 variables:
#>  $ freq : num  0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 ...
#>  $ pgram: num  0.00444 0.0383 0.00306 0.03211 0.11917 ...

# the ordinates average to the variance, in per-Hz units
c(mean(pg$pgram), var(y) / 256)
#> [1] 0.007353254 0.007339650

# one call for many series: columns are series
pg2 <- frm_periodogram(cbind(a = y, b = rev(y)), fs = 256)
table(pg2$series)
#> 
#>   a   b 
#> 255 255 
```
