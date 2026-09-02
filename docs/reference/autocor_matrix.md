# Estimated within-group residual correlation matrix

The correlation matrix `R` of the
[`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `arma()`, `cosy()` or
`unstr()` term of a fit, over the model's time levels. The residual
covariance of a group is `D R D` restricted to the time points that
group has, with `D` the diagonal matrix of that group's `sigma` values,
so `R` plus [`sigma()`](https://rdrr.io/r/stats/sigma.html) (or
`predict(dpar = "sigma")`) describes the whole residual structure.

## Usage

``` r
autocor_matrix(fit, resp = NULL)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- resp:

  Response name, for a multivariate fit.

## Value

A correlation matrix with the time levels as dimnames, or `NULL` when
the fit has no residual correlation term.

## See also

[`rescor_matrix()`](https://aforren1.github.io/frmtmb/reference/rescor_matrix.md)
for the ACROSS-response residual correlation of a `rescor = TRUE` fit,
which is a different structure and cannot be combined with this one.

## Examples

``` r
set.seed(1)
d <- expand.grid(week = 1:5, subj = factor(1:25))
e <- as.vector(apply(matrix(rnorm(125), 5, 25), 2, function(z) {
  as.vector(stats::filter(z, 0.6, "recursive"))
}))
d$x <- rnorm(125)
d$y <- 1 + 0.5 * d$x + e
fit <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
           data = d)
autocor_matrix(fit)
#>            1         2         3         4          5
#> 1 1.00000000 0.5480416 0.3003496 0.1646041 0.09020986
#> 2 0.54804158 1.0000000 0.5480416 0.3003496 0.16460405
#> 3 0.30034957 0.5480416 1.0000000 0.5480416 0.30034957
#> 4 0.16460405 0.3003496 0.5480416 1.0000000 0.54804158
#> 5 0.09020986 0.1646041 0.3003496 0.5480416 1.00000000
```
