# Build a family by name, with a link for any of its parameters

The frmtmb family constructors take `link` and a `link_<dpar>` for each
of their other distributional parameters, so
`student(link_sigma = "softplus")` needs nothing else and this is not
the usual way to reach a link.

## Usage

``` r
frm_family(family, ...)
```

## Arguments

- family:

  Family name, as a single string.

- ...:

  `link`, and any `link_<dpar>` argument the named family takes. Every
  one has to be named.

## Value

A `frmtmb_family` object, ready for
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md).

## Details

It exists for the families frmtmb has no constructor of its own for.
[`gaussian()`](https://rdrr.io/r/stats/family.html),
[`poisson()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html) and
[`Gamma()`](https://rdrr.io/r/stats/family.html) come from 'stats',
which knows nothing of a `link_sigma` and refuses several links brms
allows for the mean. frmtmb does NOT shadow them: a
[`gaussian()`](https://rdrr.io/r/stats/family.html) that answered a
`frmtmb_family` would break every
[`glm()`](https://rdrr.io/r/stats/glm.html) call in an attached session.
Naming the family instead reaches every argument, which is what brms
does with `brmsfamily()`.

## See also

[frmtmb-links](https://aforren1.github.io/frmtmb/reference/frmtmb-links.md)
for the links and which parameter takes which,
[frmtmb-families](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
for the constructors.

## Examples

``` r
# sigma on a softplus: gaussian() from 'stats' has no link_sigma
frm_family("gaussian", link_sigma = "softplus")
#> <frmtmb family> gaussian
#>   dpars: mu (identity), sigma (softplus)

# a mean link stats::poisson() refuses
frm_family("poisson", link = "softplus")
#> <frmtmb family> poisson
#>   dpars: mu (softplus)

set.seed(2)
d <- data.frame(x = rnorm(60))
d$y <- rnorm(60, 1 + d$x, exp(0.3 * d$x))
frm(bf(y ~ x, sigma ~ x),
    family = frm_family("gaussian", link_sigma = "softplus"), data = d)
#> frmtmb fit: y ~ x 
#> Family: gaussian   Method: ML 
#>  Links: mu = identity; sigma = softplus
#> 
#> logLik: -92.517  AIC: 193.034  nobs: 60 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)           x 
#>      0.9444      0.9301 
#>  sigma:
#> (Intercept)           x 
#>      0.7351      0.4798 
```
