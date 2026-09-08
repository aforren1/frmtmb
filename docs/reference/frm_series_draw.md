# Draw a time series from a fitted spectrum

The inverse of the transform a spectral fit is built on, and the verb to
reach for when [`simulate()`](https://rdrr.io/r/stats/simulate.html)
gave ordinates where a series was wanted. It draws a stationary Gaussian
series whose spectral density is the one the fit estimated, by giving
each Fourier frequency an independent complex Gaussian amplitude with
the fitted variance and transforming back.

## Usage

``` r
frm_series_draw(
  object,
  nsim = 1,
  seed = NULL,
  newdata = NULL,
  freq = "freq",
  ...
)
```

## Arguments

- object:

  A `frmtmb_fit` with a
  [`whittle()`](https://aforren1.github.io/frmtmb/reference/whittle.md)
  family.

- nsim:

  Number of series to draw.

- seed:

  Optional integer seed.

- newdata:

  Rows to take the spectrum from, for a fit that pools several series.
  Defaults to the fitted rows.

- freq:

  The frequencies, either as a column name in the data or as the numbers
  themselves. A model written in terms of a transformed frequency keeps
  no `freq` column, so passing `pg$freq` is the general route.

- ...:

  Passed to [`predict()`](https://rdrr.io/r/stats/predict.html), for
  example `re.form`.

## Value

A `ts` matrix with `2 * nf + 1` rows and `nsim` columns, at the sampling
rate the grid implies.

## Details

A simulated PERIODOGRAM cannot be inverted: the phases are not in it.
What is well defined is a draw from the fitted MODEL, which is what this
is, and it is a different object - two calls with the same fitted
spectrum give two unrelated series.

The length is fixed by the grid at `2 * nf + 1`, an odd number, so that
the drawn series has no Nyquist ordinate to invent and exactly zero
mean. The sampling rate follows from the grid spacing, and the result
carries it as a `ts` frequency, so
`frm_periodogram(frm_series_draw(fit))` lands back on the same grid.

## See also

[`whittle()`](https://aforren1.github.io/frmtmb/reference/whittle.md),
[`frm_periodogram()`](https://aforren1.github.io/frmtmb/reference/frm_periodogram.md).

## Examples

``` r
set.seed(1)
y <- as.numeric(arima.sim(list(ar = 0.6), 512))
pg <- frm_periodogram(y, fs = 128)
fit <- frm(bf(pgram ~ s(freq, k = 8)), family = whittle(), data = pg)

# a new series with the same estimated spectrum, not the same series
z <- frm_series_draw(fit, nsim = 2, seed = 1)
dim(z)
#> [1] 511   2
stats::frequency(z)
#> [1] 127.75
```
