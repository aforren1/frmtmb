# Build a family by name, with a link for any of its parameters

`brmsfamily()` is brms's generic family constructor. It names the family
with a string, takes the link for the mean as the second argument, and
takes a `link_<dpar>` argument for each of the family's other
distributional parameters.

## Usage

``` r
brmsfamily(family, link = NULL, ...)
```

## Arguments

- family:

  Family name, as a single string.

- link:

  The link for the mean, quoted or not. `NULL` gives the family's
  default. See
  [frmtmb-links](https://aforren1.github.io/frmtmb/reference/frmtmb-links.md)
  for the set each family takes.

- ...:

  Any `link_<dpar>` argument the named family takes, and its other
  constructor arguments. Every one has to be named.

## Value

A `frmtmb_family` object, ready for
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md).

## Details

The family constructors take the same arguments, so
`student(link_sigma = "softplus")` does not need it. They exist for the
families frmtmb has no constructor of its own for.
[`gaussian()`](https://rdrr.io/r/stats/family.html),
[`poisson()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html) and
[`Gamma()`](https://rdrr.io/r/stats/family.html) come from 'stats',
which knows nothing of a `link_sigma` and refuses several links brms
allows for the mean. frmtmb does NOT shadow them: a
[`gaussian()`](https://rdrr.io/r/stats/family.html) that answered a
`frmtmb_family` would break every
[`glm()`](https://rdrr.io/r/stats/glm.html) call in an attached session.

The name follows brms's spelling rules: case does not matter, `"normal"`
is `"gaussian"`, `"com_poisson"` is `"compois"`, and the prefixes
`"zi_"` and `"hu_"` stand for `"zero_inflated_"` and `"hurdle_"`. The
link may be unquoted, as in `brmsfamily("gaussian", inverse)`.

The family and its link as one character vector, `c("weibull", "log")`,
is the form brms accepts for `family =` in a model call, and
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) accepts it
there. `brmsfamily()` refuses it, as brms does.

## Differences from brms

An argument the family does not have is refused by name. brms ignores
it, so `brmsfamily("poisson", link_sigma = "log")` builds a poisson
family there. `threshold` and `refcat` are not arguments.

## See also

[frmtmb-links](https://aforren1.github.io/frmtmb/reference/frmtmb-links.md)
for the links and which parameter takes which,
[frmtmb-families](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
for the constructors.

## Examples

``` r
# sigma on a softplus: gaussian() from 'stats' has no link_sigma
brmsfamily("gaussian", link_sigma = "softplus")
#> 
#> Family: gaussian 
#> Link function: identity 
#> Parameters: mu (identity), sigma (softplus)
#> 

# a mean link stats::poisson() refuses, and brms's short names
brmsfamily("poisson", softplus)$link
#> [1] "softplus"
brmsfamily("zi_poisson")$link_zi
#> [1] "logit"

set.seed(2)
d <- data.frame(x = rnorm(60))
d$y <- rnorm(60, 1 + d$x, exp(0.3 * d$x))
frm(bf(y ~ x, sigma ~ x),
    family = brmsfamily("gaussian", link_sigma = "softplus"), data = d)
#>  Family: gaussian 
#>  Links: mu = identity; sigma = softplus
#> 
#> Formula: y ~ x 
#>    Data: d (Number of observations: 60) 
#>  Method: ML   logLik: -92.517   AIC: 193.034   BIC: 201.411 
#> 
#> Regression Coefficients:
#>                 Estimate Est.Error l-95% CI u-95% CI z value Pr(>|z|)
#> Intercept           0.94      0.15     0.64     1.25    6.15  7.8e-10
#> sigma_Intercept     0.74      0.16     0.43     1.04    4.73  2.2e-06
#> x                   0.93      0.12     0.70     1.16    7.77  8.0e-15
#> sigma_x             0.48      0.13     0.22     0.74    3.56  0.00037
```
