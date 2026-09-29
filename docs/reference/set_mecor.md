# Control the correlation of noise-free latent variables

Several [me()](https://aforren1.github.io/frmtmb/reference/frmtmb-me.md)
terms model their latent (noise-free) values as correlated by default,
which is brms's default too. Add `set_mecor(FALSE)` to a formula to
model them as independent.

## Usage

``` r
set_mecor(mecor = TRUE)
```

## Arguments

- mecor:

  `TRUE` (the default) to estimate the correlation between the latent
  values of the
  [`me()`](https://aforren1.github.io/frmtmb/reference/frmtmb-me.md)
  terms that share a grouping, `FALSE` to fix it at zero.

## Value

An object of class `frmtmb_mecor`, which `+` adds to a formula built
with [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) or
[`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md).

## See also

[frmtmb-me](https://aforren1.github.io/frmtmb/reference/frmtmb-me.md)
for the model that
[`me()`](https://aforren1.github.io/frmtmb/reference/frmtmb-me.md)
specifies.

## Examples

``` r
set.seed(3)
n <- 120
tx <- rnorm(n)
tz <- 0.6 * tx + rnorm(n, 0, 0.8)
d <- data.frame(x = tx + rnorm(n, 0, 0.3), z = tz + rnorm(n, 0, 0.3),
                sx = 0.3, sz = 0.3)
d$y <- 1 + 0.5 * tx - 0.4 * tz + rnorm(n, 0, 0.5)

# independent latent values
f0 <- bf(y ~ me(x, sx) + me(z, sz)) + set_mecor(FALSE)
fit0 <- frm(f0, data = d)
```
