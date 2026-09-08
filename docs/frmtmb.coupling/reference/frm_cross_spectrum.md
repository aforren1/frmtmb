# The cross-spectrum of a pair of signals, as rows a model can read

Turns two time series into one row per frequency holding the Hermitian
cross-spectral matrix and the number of independent complex draws that
went into it. Those rows are the response
[`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md)
consumes.

## Usage

``` r
frm_cross_spectrum(
  x,
  y,
  sfreq = 1,
  segments = 8L,
  tapers = 1L,
  smooth = 1L,
  frange = NULL
)
```

## Arguments

- x, y:

  The two signals. Either two numeric vectors of equal length, or two
  matrices with the same dimensions whose COLUMNS are UNITS (subjects,
  trials, sessions) and never channels; a matrix pair adds an `id`
  column naming the column. A three-column pair is three units of the
  same two signals, not a three-channel recording:
  [`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md)
  models a channel PAIR, and a third channel has no route in. To relate
  three channels, fit the three pairs separately and say so.

- sfreq:

  Sampling rate in Hz. `freq` comes back in Hz; the default of 1 returns
  cycles per sample.

- segments:

  Number of disjoint, non-overlapping blocks the record is cut into.
  Each contributes one complex draw per frequency.

- tapers:

  Number of orthogonal sine tapers applied within each segment. Each
  contributes one further draw. `1`, the default, means no taper at all
  rather than the first sine taper.

- smooth:

  Number of adjacent Fourier bins averaged together. The returned `freq`
  is the middle of each group. `tapers` and `smooth` may not both be
  above 1; see "How many degrees of freedom".

- frange:

  Optional `c(low, high)` in the same units as `freq`, applied after
  everything else.

## Value

A data frame with one row per retained frequency and columns `freq`,
`w11`, `w22`, `w12r`, `w12i` and `n`, plus `id` for the matrix form.
`w11` and `w22` are the two auto-spectra summed over draws, `w12r` and
`w12i` the real and imaginary parts of the summed cross-spectrum, and
`n` the degrees of freedom. The matrix is
`[[w11, w12r + 1i w12i], [w12r - 1i w12i, w22]]`, and it is the SUM over
draws rather than the average, which is the scaling the complex Wishart
density is written for. Each draw carries the `1/N` periodogram
convention, so `w11 / n` estimates a spectral density rather than a
squared transform. Any constant rescaling of the pair is absorbed by the
two fitted powers and leaves coherence and phase untouched, so the
convention matters for reading `mu` and `pow2` and for nothing else.

## What is dropped, and why

Frequency zero and the Nyquist frequency are removed. Their Fourier
coefficients are real rather than complex, so the pair's periodogram
there is real Wishart with half the degrees of freedom and a different
density. Keeping them is a mistake a user makes once and never notices,
because the fit still converges.

The mean is removed from every segment before its transform, which is
what makes frequency zero uninformative in the first place.

## How many degrees of freedom, and where to get them

The density needs `n >= 2`, because a cross-periodogram from a single
draw is rank one: its coherence is exactly 1 by arithmetic, whatever the
signals did. Below 4 the estimate is nearly worthless even though it
exists.

Segments, tapers and smoothing all buy real degrees of freedom. At a
true coherence of zero the naive estimate has mean exactly `1/n`, so
`1 / mean(coherence)` reads the effective count straight off a
simulation. Every accepted configuration was measured that way, 4000
replicates on white noise, nominal `n = 8` throughout:

|                            |            |
|----------------------------|------------|
| configuration              | measured n |
| `segments = 8`             | 7.92       |
| `segments = 4, smooth = 2` | 8.17       |
| `segments = 2, smooth = 4` | 7.99       |
| `segments = 1, smooth = 8` | 7.98       |
| `segments = 4, tapers = 2` | 8.03       |
| `segments = 2, tapers = 4` | 7.76       |
| `segments = 1, tapers = 8` | 7.98       |

Two combinations are refused rather than counted wrong.

**Tapers and smoothing together.** Both widen the same spectral window,
so their product is not the degrees of freedom: 4 tapers smoothed over 4
bins measures 5.80 against a nominal 16, and 2 tapers over 2 bins
measures 2.70 against 4. Raise `segments` instead.

**Overlapping segments.** Overlap raises the nominal count without
raising the independent one, so `n` would be a lie and every standard
error downstream would be too small.

Frequency smoothing assumes the spectrum is flat across the bins it
averages, so it is the route with a cost that a white-noise measurement
cannot see. Measured on an AR(1) spectrum with `phi = 0.9`, series
length 4096, every retained bin of 250 replicates pooled, **that cost is
not detectable at these widths**: nominal 2, 4, 8 and 16 deliver 1.994,
3.993, 8.008 and 16.069, against 2.002, 3.992, 8.011 and 15.965 on white
noise. It appears only at `smooth = 32`, which delivers 31.5.

An earlier draft of this package reported a 5 percent loss at
`smooth = 16` on that spectrum. That number came from reading one
frequency bin per replicate instead of all of them, and it did not
survive replication. Smoothing is safe at the widths anyone uses; what
is not safe is smoothing across a band where the COHERENCE varies, which
is a modelling error rather than a degrees-of-freedom one and no count
will catch it.

## Examples

``` r
set.seed(1)
src <- rnorm(2048)
a <- src + rnorm(2048)
b <- 0.8 * src + rnorm(2048)
xs <- frm_cross_spectrum(a, b, sfreq = 256, segments = 8)
head(xs)
#>   freq       w11      w22       w12r        w12i n
#> 1    1 15.587842 11.56755 10.1503351  1.49247664 8
#> 2    2  9.024724 11.10212  5.4034368  0.57131947 8
#> 3    3  9.366658 11.07669  0.1945243 -0.85311662 8
#> 4    4 17.254757 15.84840  7.4288153  0.07170445 8
#> 5    5 22.136941 17.15566  7.6355380 -9.01984914 8
#> 6    6 23.842360 11.99445  1.1244380 -3.54011558 8
```
