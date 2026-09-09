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
  window = c("none", "hann"),
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
  three channels, fit the three pairs separately and say so. `NA` marks
  a sample the record does not have; see "Gaps in the record".

- sfreq:

  Sampling rate in Hz. `freq` comes back in Hz; the default of 1 returns
  cycles per sample.

- segments:

  How many disjoint, non-overlapping blocks the record is cut into, at
  most. Each contributes one complex draw per frequency. A record with
  gaps in it usually supplies fewer, and `n` says how many it supplied.

- tapers:

  Number of orthogonal sine tapers applied within each segment. Each
  contributes one further draw. `1`, the default, means no taper at all
  rather than the first sine taper.

- smooth:

  Number of adjacent Fourier bins averaged together. The returned `freq`
  is the middle of each group. `tapers` and `smooth` may not both be
  above 1; see "How many degrees of freedom".

- window:

  The data window applied to each segment before its transform.
  `"none"`, the default, is a boxcar. `"hann"` is the raised cosine, for
  a spectrum steep enough that the untapered transform's high
  frequencies are leakage from its low ones. It buys no degrees of
  freedom and costs none, and it cannot be combined with `tapers` or
  `smooth`; see "How many degrees of freedom".

- frange:

  Optional `c(low, high)` in the same units as `freq`, applied after
  everything else.

## Value

A data frame with one row per retained frequency and columns `freq`,
`w11`, `w22`, `w12r`, `w12i` and `n`, plus `id` for the matrix form.
`w11` and `w22` are the two auto-spectra summed over draws, `w12r` and
`w12i` the real and imaginary parts of the summed cross-spectrum, and
`n` the degrees of freedom the record actually supplied. The matrix is
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

## Gaps in the record

`NA` in either signal marks a sample the pair does not have, which is
what artifact rejection leaves behind. Such a sample is not interpolated
and not skipped: the record is cut at it, and the segments are laid
inside the clean spans that remain, so no transform ever crosses a gap.
An artifact-rejected recording is then one call rather than one call per
surviving span and a hand-written sum.

A sample is usable only where BOTH signals have it, because the
cross-spectrum is a property of the pair. `NaN` counts as `NA`; an
infinity is refused, because it is a value the arithmetic cannot use
rather than a value the record is missing.

The segment LENGTH comes from the usable sample count divided by
`segments`, and each clean span then supplies as many whole segments as
fit in it, in record order, up to `segments` in total. The samples left
over at the end of each span are dropped, and a span shorter than one
segment supplies nothing at all.

So a record with ONE clean span always yields the full `segments`,
whatever was rejected from its ends, because the length divides down. A
record broken into SEVERAL usually yields fewer, since each span drops
its own remainder. On such a record, asking for FEWER segments can be
refused where more would be accepted: the length comes from the global
usable count, so a small `segments` can ask for a segment longer than
any surviving span. Three spans of 700, 1043 and 699 are refused at
`segments = 2`, which wants 1221 samples in a row, and give 3, 7 and 14
segments at 4, 8 and 16.

That is not hidden: `n` is the count the record supplied, it is per-row
data the density reads, and every standard error downstream is formed
from it. Nothing warns, because on a record with gaps the shortfall is
the normal case rather than a mistake.

In the matrix form each column is its own record, so a column with more
rejected samples than its neighbors gets a shorter segment, its own
frequency grid and its own `n`. That is correct and it is why `freq` is
a column of the frame rather than an attribute of it: a model reads
`freq` per row.

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

|                                 |            |
|---------------------------------|------------|
| configuration                   | measured n |
| `segments = 8`                  | 7.92       |
| `segments = 4, smooth = 2`      | 8.17       |
| `segments = 2, smooth = 4`      | 7.99       |
| `segments = 1, smooth = 8`      | 7.98       |
| `segments = 4, tapers = 2`      | 8.03       |
| `segments = 2, tapers = 4`      | 7.76       |
| `segments = 1, tapers = 8`      | 7.98       |
| `segments = 8, window = "hann"` | 8.02       |

Three combinations are refused rather than counted wrong.

**Tapers and smoothing together.** Both widen the same spectral window,
so their product is not the degrees of freedom: 4 tapers smoothed over 4
bins measures 5.80 against a nominal 16, and 2 tapers over 2 bins
measures 2.70 against 4. Raise `segments` instead.

**A window and smoothing together**, for the same reason and by the same
arithmetic. A Hann window makes each ordinate a weighted sum of its two
neighbors and itself, so adjacent bins are correlated and averaging them
adds less than it claims. Measured at a true coherence of zero, 3000
replicates a cell, nominal `n = 8` throughout: 4 segments smoothed over
2 bins delivers 5.66, 2 segments over 4 bins delivers 4.98 and 1 segment
over 8 bins delivers 4.68, against 8.19, 8.02 and 8.05 for the same
three without the window.

Declaring the shortfall instead of refusing it was considered and
rejected. There is no one number to declare: the same window loses a
different amount at each smoothing width, and frmtmb's own
[`whittle()`](https://aforren1.github.io/frmtmb/reference/whittle.html)
measured the same thing from the other side, where an honestly declared
equivalent degrees of freedom does not rescue a smoothed spectrum from
its raw-periodogram check at kernel widths of 7 and above
(`dev/reviews/2026-09-08-spectral2.md`). The correlation is a property
of the estimate, not of the number attached to it.

**A window and tapers together.** Both are a taper on the same segment,
and applying one over the other is neither.

**Overlapping segments** stay refused. Overlap raises the nominal count
without raising the independent one, so `n` would be a lie and every
standard error downstream would be too small.

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
is a modeling error rather than a degrees-of-freedom one and no count
will catch it.

## What a window does to the rows, which `n` does not say

`n` counts the draws behind ONE row and a window leaves that count
alone, which the table above measures. What a window and what `tapers`
both change is the relation BETWEEN rows, and `n` says nothing about
that.

Measured on the log of `w11` from this function, 400 replicates of white
noise at 1024 samples, at calls this function accepts. The last column
is the variance inflation factor of a mean over the rows, `1 + 2 * sum`
of the positive correlations, which is the number a smooth in frequency
actually pays. It is summed over twelve lags, which matters for one row:
at `tapers = 8` the correlation is still positive past lag 4, and
stopping there gives 6.45 instead of 8.02. Every other row has died by
lag 3 and is the same either way.

|                                 |        |        |        |           |
|---------------------------------|--------|--------|--------|-----------|
| configuration                   | lag 1  | lag 2  | lag 3  | inflation |
| `segments = 8`                  | -0.010 | -0.015 | -0.021 | 1.00      |
| `segments = 8, window = "hann"` | 0.397  | -0.009 | -0.027 | 1.79      |
| `segments = 8, tapers = 2`      | 0.525  | 0.104  | -0.035 | 2.24      |
| `segments = 4, tapers = 4`      | 0.745  | 0.502  | 0.271  | 4.07      |
| `segments = 1, tapers = 8`      | 0.876  | 0.744  | 0.615  | 8.02      |

So a Hann-windowed frame carries about one independent frequency in 1.8,
and `tapers` costs MORE rather than less: a `tapers = 4` frame carries
about one in four, and eight tapers cost a factor of eight, with
ordinates still correlated 0.6 three bins apart. An `s(freq)` smooth
fitted on either has fewer effective points than it has rows, so its
band is optimistic. By how much is not measured here: the inflation
factor says how much information the rows carry, not what a penalized
smooth does with it.

None of that touches `n`, which is right in every one of these
configurations: the per-row degrees of freedom and the correlation
between rows are different quantities, and `n` claims only the first.

Use a window or a taper where leakage is the larger error, which is a
steep spectrum, and not by default.

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
