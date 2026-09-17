# Factor with numeric-coded levels for coordinate covariance structures

`ou()` and the spatial structures
([`exp()`](https://rdrr.io/r/base/Log.html), `gau()`, `mat()`) need the
positions of the term levels. `num_factor(x)` (one dimension) or
`num_factor(x, y)` (planar coordinates) encodes them in the level labels
the same way
[`glmmTMB::numFactor()`](https://rdrr.io/pkg/glmmTMB/man/numFactor.html)
does, so factors created by either function work.

## Usage

``` r
num_factor(x, y = NULL)
```

## Arguments

- x:

  Numeric positions (times, coordinates).

- y:

  Optional second coordinate.

## Value

A factor whose levels encode the unique positions.

## Examples

``` r
# unequally spaced observation times, kept as distances
num_factor(c(0, 1.5, 4))
#> [1] (0)   (1.5) (4)  
#> Levels: (0) (1.5) (4)

# planar coordinates for a spatial covariance
levels(num_factor(rep(1:3, 3), rep(1:3, each = 3)))
#> [1] "(1,1)" "(1,2)" "(1,3)" "(2,1)" "(2,2)" "(2,3)" "(3,1)" "(3,2)" "(3,3)"

# ou() reads the distances out of the level labels; a plain factor
# would only give it an ordering
set.seed(1)
tim <- c(0, 1, 1.5, 3)
n_g <- 40
S <- 0.9^2 * exp(-1.2 * abs(outer(tim, tim, "-")))
u <- matrix(rnorm(n_g * length(tim)), n_g) %*% chol(S)
dd <- data.frame(
  y = 1 + as.vector(t(u)) + rnorm(n_g * length(tim), 0, 0.4),
  g = factor(rep(seq_len(n_g), each = length(tim))),
  tim = num_factor(rep(tim, n_g))
)
fit <- frm(bf(y ~ 1 + ou(tim + 0 | g)) + gaussian(), data = dd)
round(VarCorr(fit)$g$cor[, "Estimate", ], 3)
#>         tim0  tim1 tim1.5  tim3
#> tim0   1.000 0.235  0.114 0.013
#> tim1   0.235 1.000  0.485 0.055
#> tim1.5 0.114 0.485  1.000 0.114
#> tim3   0.013 0.055  0.114 1.000
```
