# Convert draws to a posterior draws object

Convert draws to a posterior draws object

## Usage

``` r
as_draws(x, ...)

# S3 method for class 'frmtmb_draws'
as_draws(x, ...)
```

## Arguments

- x:

  A `frmtmb_draws` object.

- ...:

  Unused.

## Value

A
[`posterior::draws_matrix`](https://mc-stan.org/posterior/reference/draws_matrix.html):
one column per sampled variable and one row per draw.
