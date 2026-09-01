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
round(VarCorr(fit)[[1]], 3)
#>          tim(0) tim(1) tim(1.5) tim(3)
#> tim(0)    0.823  0.193    0.094  0.011
#> tim(1)    0.193  0.823    0.399  0.045
#> tim(1.5)  0.094  0.399    0.823  0.094
#> tim(3)    0.011  0.045    0.094  0.823
```
