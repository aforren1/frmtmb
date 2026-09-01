# Posterior class probabilities of a mixture fit

For an ordinary
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md) or
[`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
fit, one row per observation; for a group-level mixture (`groups = ~g`),
one row per group.

## Usage

``` r
mixture_probs(fit)
```

## Arguments

- fit:

  A `frmtmb_fit` with a mixture family.

## Value

A matrix of class probabilities (rows sum to one).
