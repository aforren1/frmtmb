# Check a custom family's log-density for AD safety

Tapes the family's `lpdf` on test values and compares the AD gradient
against central finite differences. A mismatch usually means the lpdf
uses operations the tape cannot see (base
[`matrix()`](https://rdrr.io/r/base/matrix.html)/[`c()`](https://rdrr.io/r/base/c.html)
on advectors, branching on parameter values, `min`/`max`, clamping).

## Usage

``` r
check_custom_family(family, y, dpars, aterms = list(), tol = 1e-04)
```

## Arguments

- family:

  A `frmtmb_family` (from
  [`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
  /
  [`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)).

- y:

  A response vector of test data.

- dpars:

  Named list of numeric test values, one entry per dpar (each of length
  1 or `length(y)`).

- aterms:

  Named list of addition-term values (e.g. `trials`).

- tol:

  Maximum relative gradient error.

## Value

Invisibly `TRUE`; signals an error on failure.

## Examples

``` r
set.seed(1)
y <- rpois(50, 3)

# a hand-written poisson: check it before fitting anything with it
ok <- custom_family(
  "my_poisson", dpars = "mu", links = list(mu = "log"),
  lpdf = function(y, dpars, aterms) {
    y * log(dpars$mu) - dpars$mu - lgamma(y + 1)
  },
  type = "discrete"
)
check_custom_family(ok, y = y, dpars = list(mu = rep(2.5, 50)))

# base matrix() strips the advector class, so the tape sees constants
# and the gradient is silently wrong. The check catches it.
bad <- custom_family(
  "bad", dpars = "mu", links = list(mu = "log"),
  lpdf = function(y, dpars, aterms) {
    m <- matrix(dpars$mu, ncol = 1)
    y * log(m[, 1]) - m[, 1] - lgamma(y + 1)
  },
  type = "discrete"
)
try(check_custom_family(bad, y = y, dpars = list(mu = rep(2.5, 50))))
#> Error : Failed to tape the lpdf: Invalid argument to 'advector' (lost class attribute?). Typical cause: base matrix()/c() stripping the advector class, or branching on parameter values
```
