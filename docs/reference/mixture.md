# Finite mixture families

`mixture(fam1, fam2, ...)` builds a K-component mixture: each component
keeps its own distributional parameters, suffixed by the component index
(`mu1`, `sigma1`, `mu2`, ...), and the mixing proportions come from
`theta1 ... theta{K-1}` (multinomial-logit against the last component,
each with its own linear predictor - so mixing weights may depend on
covariates). The main model formula applies to every component mean;
override per component with `bf(y ~ x, mu2 ~ 1)`.

## Usage

``` r
mixture(..., groups = NULL)
```

## Arguments

- ...:

  Two or more component families.

- groups:

  Optional one-sided formula naming the latent-class grouping factor.

## Value

A `frmtmb_family`.

## Details

The likelihood is a parameter-branch-free logsumexp, so Laplace
machinery is untouched; the usual finite-mixture ML caveats apply
instead: the likelihood is invariant to component relabeling, so the
component means are initialized on spread-out response quantiles, and
multimodality is real (compare starts, or order the intercepts with
bounds: `set_prior("", class = "Intercept", dpar = ..., lb = )` per
component). Component families with extra parameters (ordinal) are not
supported.

With `groups = ~g` the mixture moves to the group level (latent
classes): every observation of a group shares one class draw, and the
marginal likelihood sums the class assignment per group. Continuous
random effects, smooths, and gp() terms are allowed in the component
formulas - the class sum happens conditional on the latent effects, so
one Laplace approximation integrates them (growth-mixture models).
Random effects written in a component formula are class-specific by
construction; the Laplace approximation of the class-mixture integrand
is not exact even for gaussian responses (a fraction of a log-likelihood
unit in typical well-separated problems). `quadrature = TRUE` makes the
integral numerically exact when the per-group integrand is univariate
(one scalar random intercept, in one class); with class-specific
intercepts in several classes the coordinates couple and quadrature
remains approximate - use
[`frmtmb.sample::check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.html)
to judge. Mixing-weight predictors are evaluated at each group's first
row (use group-constant covariates).
[`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)
returns the posterior class probabilities per group (or per observation
for ordinary mixtures), conditional on the random-effect modes.

## Examples

``` r
# two well-separated gaussian components
set.seed(3)
dd <- data.frame(y = c(rnorm(80, 0, 1), rnorm(80, 5, 1)),
                 x = rnorm(160))
fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
# one mu and sigma per component, plus the mixing weight theta1
fixef(fit)
#> $mu1
#>  (Intercept) 
#> -0.002346141 
#> 
#> $sigma1
#> (Intercept) 
#>  -0.1298411 
#> 
#> $mu2
#> (Intercept) 
#>    4.999987 
#> 
#> $sigma2
#> (Intercept) 
#>  0.01009985 
#> 
#> $theta1
#> (Intercept) 
#> 0.004758874 
#> 
# posterior class probability per observation
head(mixture_probs(fit))
#>         class1       class2
#> [1,] 1.0000000 4.288987e-08
#> [2,] 0.9999990 1.000118e-06
#> [3,] 0.9999851 1.488386e-05
#> [4,] 1.0000000 1.798825e-08
#> [5,] 0.9999891 1.087896e-05
#> [6,] 0.9999952 4.800314e-06

# the mixing weight can take its own predictor
frm(bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()), data = dd)
#> frmtmb fit: y ~ 1 
#> Family: mixture(gaussian, gaussian)   Method: ML 
#> logLik: -327.164  AIC: 666.328  nobs: 160 
#> 
#> Fixed effects:
#>  mu1:
#> (Intercept) 
#>   -0.002575 
#>  sigma1:
#> (Intercept) 
#>     -0.1302 
#>  mu2:
#> (Intercept) 
#>           5 
#>  sigma2:
#> (Intercept) 
#>     0.01038 
#>  theta1:
#> (Intercept)           x 
#>     0.01004    -0.07069 

# \donttest{
# latent classes: every observation of a group shares one class
set.seed(4)
n_g <- 40
cls <- rep(c(1, 2), each = n_g / 2)
dg <- data.frame(g = factor(rep(seq_len(n_g), each = 5)))
dg$y <- rnorm(nrow(dg), c(0, 4)[cls[as.integer(dg$g)]], 1)
fg <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
          data = dg)
head(mixture_probs(fg))   # one row per group, not per observation
#>   class1       class2
#> 1      1 1.262122e-12
#> 2      1 2.820791e-12
#> 3      1 8.436732e-16
#> 4      1 7.884398e-16
#> 5      1 4.031519e-09
#> 6      1 1.439555e-13
# }
```
