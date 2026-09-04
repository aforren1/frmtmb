# Posterior probabilities of the latent states

One matrix of posterior latent-state probabilities for a fit whose
family carries a
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
with a `latent_probs` slot: the probability of each state or class given
everything the model observed. Rows sum to one.

## Usage

``` r
latent_probs(fit, ...)
```

## Arguments

- fit:

  A `frmtmb_fit` whose family declares latent states.

- ...:

  Passed to methods; none of the built-in families take any.

## Value

A numeric matrix with one named column per latent state. Its rows are
observations for a family whose class belongs to the row (`lca()`, a
rowwise
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md))
and groups for one whose class belongs to a group
(`mixture(groups = )`).

## Details

This is one entry point for every structured family, and the family
decides what it returns.
[`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md),
`hmm_probs()` and `lca_probs()` are the same computation under
family-specific names and checks; use whichever names what you fitted.

## See also

[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md),
[`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md),
`hmm_probs()`, `lca_probs()`, `hmm_viterbi()` for the decoded path
rather than the probabilities

## Examples

``` r
set.seed(1)
dd <- data.frame(y = c(rnorm(80, 0, 1), rnorm(80, 5, 1)))
fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
head(latent_probs(fit))
#>         class1       class2
#> [1,] 0.9999999 7.623679e-08
#> [2,] 0.9999948 5.179722e-06
#> [3,] 1.0000000 2.623868e-08
#> [4,] 0.9888568 1.114321e-02
#> [5,] 0.9999888 1.123676e-05
#> [6,] 1.0000000 2.833862e-08
```
