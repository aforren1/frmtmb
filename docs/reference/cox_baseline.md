# The fitted baseline-hazard simplex of a `cox()` fit.

The fitted baseline-hazard simplex of a
[`cox()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
fit.

## Usage

``` r
cox_baseline(fit)
```

## Arguments

- fit:

  A `frmtmb_fit` with a
  [`cox()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  family.

## Value

The `Kbhaz` M-spline weights, summing to one.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(200))
dd$time <- rexp(200, exp(-0.5 + 0.7 * dd$x))
fit <- frm(bf(time ~ x), family = cox(), data = dd)
cox_baseline(fit)
#>         s1         s2         s3         s4         s5 
#> 0.01712166 0.21161600 0.11189645 0.06392626 0.59543964 
```
