# The cross-spectra of several channel pairs, stacked with a `pair` factor

[`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md)
models ONE channel pair. A recording with more channels has a pair for
every two of them, and the questions asked of it are about all of those
pairs at once. `frm_cross_pairs()` runs
[`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md)
on each pair and stacks the results, one block per pair, with a `pair`
factor that a model can use.

## Usage

``` r
frm_cross_pairs(X, pairs = NULL, ...)
```

## Arguments

- X:

  The recording: a numeric matrix or a data frame of numeric columns
  whose COLUMNS ARE CHANNELS and whose rows are samples, or a list of
  such matrices, one per epoch, each with the same channels. Channel
  names come from the column names, or are `ch1`, `ch2`, and so on when
  there are none.

- pairs:

  Which pairs to form. `NULL`, the default, forms every pair of two
  different channels, in column order. Otherwise a two-column matrix or
  data frame, one row per pair, naming each channel by name or by column
  number.

- ...:

  Passed to
  [`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md):
  `sfreq`, `segments`, `tapers`, `smooth`, `window`, `frange` and, for a
  list of epochs, `group`.

## Value

The rows of
[`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md)
for every pair, with the factors `pair` (`"a-b"` for channels `a` and
`b`), `ch1` and `ch2` in front. `ch1` is the signal
[`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md)
takes as `x`, so the fitted phase is the lead of `ch1` over `ch2`.

## Many pairs in one model

Fitting each pair on its own and correcting the p-values afterwards
treats the pairs as unrelated. One model over the stacked frame lets
them share information instead:

- `coh ~ 1 + (1 | pair)` gives each pair its own coherence, shrunk
  toward the pairs' common level by as much as the data say the pairs
  differ. A pair that stands out after shrinkage stands out against the
  others, which is the comparison a multiple-comparison correction tries
  to make after the fact.

- `coh ~ s(freq, by = pair)` with `pair` also a fixed effect gives each
  pair its own coherence spectrum, and
  [`frm_coherence()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_coherence.md)
  reads each one off with a band.

- `coh ~ 0 + pair`, with `mu`, `pow2` and `phase` also by `pair`, is the
  separate fits in one call: the likelihood is then a sum of one term
  per pair, so each pair's estimate is the one its own fit gives.

The power of a channel appears in every pair that includes it, so the
blocks share channels. Each block is still a correct likelihood for its
own pair; what the stacked likelihood does not model is the dependence
between two blocks that share a channel. Read the standard errors of a
stacked fit with that in mind.

## See also

[`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md),
[`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md).

## Examples

``` r
set.seed(2)
src <- rnorm(2048)
X <- cbind(Fz = src + rnorm(2048), Cz = 0.8 * src + rnorm(2048),
           Pz = rnorm(2048))
xp <- frm_cross_pairs(X, sfreq = 256, segments = 8)
table(xp$pair)
#> 
#> Fz-Cz Fz-Pz Cz-Pz 
#>   127   127   127 
```
