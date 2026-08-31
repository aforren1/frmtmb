# Additional response families

Family constructors without a
[stats::family](https://rdrr.io/r/stats/family.html) equivalent,
following brms naming.
[`gaussian()`](https://rdrr.io/r/stats/family.html),
[`poisson()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html), and
[`Gamma()`](https://rdrr.io/r/stats/family.html) from 'stats' are
accepted directly by [`frm()`](frm.md) and `+`.

## Usage

``` r
student(link = "identity")

lognormal(link = "identity")

negbinomial(link = "log")

nbinom1(link = "log")

Beta(link = "logit")

tweedie(link = "log")

compois(link = "log")

zero_inflated_poisson(link = "log")

zero_inflated_negbinomial(link = "log")

hurdle_poisson(link = "log")

multinomial(K)

cumulative(link = "logit")

beta_binomial(link = "logit")

skew_normal(link = "identity")

exgaussian(link = "identity")

bernoulli(link = "logit")

geometric(link = "log")

exponential(link = "log")

weibull(link = "log")

shifted_lognormal(link = "identity")

hurdle_gamma(link = "log")

hurdle_lognormal(link = "identity")

zero_inflated_binomial(link = "logit")

zero_inflated_beta(link = "logit")

asym_laplace(link = "identity")

sratio(link = "logit")

cratio(link = "logit")

acat(link = "logit")
```

## Arguments

- link:

  Link for `mu`.

- K:

  For `multinomial()`: number of response categories (columns of the
  count-matrix response); category 1 is the reference.

## Value

A `frmtmb_family` object.
