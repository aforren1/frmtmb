# Extract fixed effects

Extract fixed effects

## Usage

``` r
fixef(object, ...)

# S3 method for class 'frmtmb_fit'
fixef(object, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Unused.

## Value

A named list of coefficient vectors, one per dpar.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)

# one entry per distributional parameter, each on its link scale
fit <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), data = dd)
fixef(fit)
#> $mu
#> (Intercept)           x 
#>    1.236875    0.594386 
#> 
#> $sigma
#> (Intercept)           x 
#>  -0.0959919   0.2082866 
#> 
exp(fixef(fit)$sigma[["(Intercept)"]])   # sigma is modeled on the log
#> [1] 0.9084714

# flatten to the vector confint() and hypothesis() name their rows by
unlist(fixef(fit))
#>    mu.(Intercept)              mu.x sigma.(Intercept)           sigma.x 
#>         1.2368754         0.5943860        -0.0959919         0.2082866 
```
