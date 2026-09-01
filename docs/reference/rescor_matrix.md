# Estimated residual correlation matrix (rescor fits), else NULL

Estimated residual correlation matrix (rescor fits), else NULL

## Usage

``` r
rescor_matrix(fit)
```

## Arguments

- fit:

  A `frmtmb_fit`.

## Value

A correlation matrix or `NULL`.

## Examples

``` r
set.seed(2)
n <- 80
dd <- data.frame(x = rnorm(n))
# two responses that share a residual disturbance
e <- rnorm(n)
dd$y1 <- 1 + 0.5 * dd$x + e + rnorm(n, 0, 0.5)
dd$y2 <- 2 - 0.3 * dd$x + e + rnorm(n, 0, 0.5)
fit <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE) + gaussian(),
           data = dd)
rescor_matrix(fit)
#>           y1        y2
#> y1 1.0000000 0.8371743
#> y2 0.8371743 1.0000000

# a fit without rescor has no residual correlation to report
fit0 <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + gaussian(), data = dd)
rescor_matrix(fit0)
#> NULL
```
