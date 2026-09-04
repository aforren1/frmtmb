# The fitted item profiles of an `lca()` fit

The class-conditional item-response probabilities `pi[j, c, k]`, which
are the main output of a latent class analysis: for each item, a
`K x C_j` table of the probability of each category given each class.
This is poLCA's `$probs`.

## Usage

``` r
lca_profiles(fit)
```

## Arguments

- fit:

  A `frmtmb_fit` with an
  [`lca()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.md)
  family.

## Value

A named list with one `K x C_j` probability matrix per item, of class
`frmtmb_lca_profiles`, carrying the estimated class sizes (the mean
prior class probabilities, poLCA's `$P`) in the `class_sizes` attribute.

## See also

[`lca()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.md),
[`lca_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_probs.md)

## Examples

``` r
set.seed(2)
n <- 200
cl <- rbinom(n, 1, 0.5) + 1
pr <- rbind(c(0.9, 0.85, 0.8), c(0.1, 0.15, 0.2))
Y <- matrix(0L, n, 3)
for (j in 1:3) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
dd <- data.frame(row = seq_len(n))
dd$Y <- Y
fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)

pf <- lca_profiles(fit)
pf
#> <lca profiles> 2 classes, 3 items
#> 
#> Estimated class sizes (mean prior probability):
#> class1 class2 
#>  0.437  0.563 
#> 
#> item1:
#>          cat1   cat2
#> class1 0.9465 0.0535
#> class2 0.1179 0.8821
#> 
#> item2:
#>          cat1   cat2
#> class1 0.8745 0.1255
#> class2 0.2005 0.7995
#> 
#> item3:
#>          cat1   cat2
#> class1 0.8209 0.1791
#> class2 0.2332 0.7668
pf[[1]]                       # item 1's K x C table
#>             cat1       cat2
#> class1 0.9464869 0.05351314
#> class2 0.1179051 0.88209491
attr(pf, "class_sizes")
#>    class1    class2 
#> 0.4370062 0.5629938 
```
