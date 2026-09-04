# Expose a model's compiled functions

brms compiles Stan functions and `expose_functions()` makes them
callable from R. frmtmb compiles nothing: the model is an R closure
built by
[`build_objective()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.md)
and differentiated by RTMB, and a custom family's density is the plain R
function handed to
[`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md).
The method is defined so that a ported brms script gets that reason
rather than "could not find function".

## Usage

``` r
expose_functions(x, ...)

# S3 method for class 'frmtmb_fit'
expose_functions(x, ...)
```

## Arguments

- x:

  A `frmtmb_fit`.

- ...:

  Ignored; this method always stops.

## Value

This function never returns; it signals an error.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(40))
dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)
try(expose_functions(fit))
#> Error : expose_functions() has no Stan program to read on a frmtmb fit: a custom family's lpdf is the plain R function handed to custom_family(), callable as it is
```
