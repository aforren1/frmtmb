# Sample from a frmtmb fit with tmbstan (NUTS)

Hands the fitted RTMB object to
[`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html);
from Stan's point of view the model is the (Laplace-free) joint density,
so all parameters including random effects are sampled unless
`laplace = TRUE` is passed through.

## Usage

``` r
as_tmbstan(fit, ...)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- ...:

  Passed to
  [`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html)
  (chains, iter, laplace, ...).

## Value

A `stanfit` object.
