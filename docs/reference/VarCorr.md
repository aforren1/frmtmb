# Extract random-effect covariance matrices

Extract random-effect covariance matrices

## Usage

``` r
VarCorr(x, ...)

# S3 method for class 'frmtmb_fit'
VarCorr(x, ...)
```

## Arguments

- x:

  A `frmtmb_fit`.

- ...:

  Unused.

## Value

A named list of covariance matrices, one per random-effect term. The
names are the term labels, which can repeat when two blocks deparse the
same way (`(1 | gr(id, cov = A)) + (1 | id)`, the animal model's genetic
and permanent-environment terms). Index by position, not by name, when
that is possible in your model.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.4))
dd$y <- rnorm(200, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)

# the printed form shows SDs and correlations, as lme4 does
VarCorr(fit)
#>   x | g 
#>         Name Std.Dev. (Intercept)
#>  (Intercept)  0.88616            
#>            x  0.45121      0.0527
# the stored value is the covariance matrix itself
VarCorr(fit)[["x | g"]]
#>             (Intercept)          x
#> (Intercept)  0.78528800 0.02108648
#> x            0.02108648 0.20359360
# tidy shape for broom.mixed-style code
as.data.frame(VarCorr(fit))
#>     grp        var1 var2       vcov      sdcor
#> 1 x | g (Intercept) <NA> 0.78528800 0.88616477
#> 2 x | g           x <NA> 0.20359360 0.45121347
#> 3 x | g (Intercept)    x 0.02108648 0.05273605
```
