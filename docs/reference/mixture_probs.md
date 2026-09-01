# Posterior class probabilities of a mixture fit

For an ordinary
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md) or
[`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
fit, one row per observation; for a group-level mixture (`groups = ~g`),
one row per group.

## Usage

``` r
mixture_probs(fit)
```

## Arguments

- fit:

  A `frmtmb_fit` with a mixture family.

## Value

A matrix of class probabilities (rows sum to one).

## Examples

``` r
set.seed(3)
dd <- data.frame(y = c(rnorm(80, 0, 1), rnorm(80, 5, 1)))
fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)

p <- mixture_probs(fit)
head(p)
#>         class1       class2
#> [1,] 1.0000000 4.288987e-08
#> [2,] 0.9999990 1.000118e-06
#> [3,] 0.9999851 1.488386e-05
#> [4,] 1.0000000 1.798825e-08
#> [5,] 0.9999891 1.087896e-05
#> [6,] 0.9999952 4.800314e-06
rowSums(p)[1:3]                    # rows sum to one
#> [1] 1 1 1

# the hard assignment, and how well it recovers the truth
cl <- max.col(p)
table(cl, truth = rep(1:2, each = 80))
#>    truth
#> cl   1  2
#>   1 80  0
#>   2  0 80
```
