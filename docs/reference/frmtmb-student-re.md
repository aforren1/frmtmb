# Student-t distributed random effects

`(x | gr(g, dist = "student"))` gives a grouping term a Student-t latent
instead of a gaussian one, which is brms's spelling
([`brms::gr()`](https://paulbuerkner.com/brms/reference/gr.html),
argument `dist`). A group far from the others then costs the variance
component much less than it does under a gaussian latent, because the
t's tail can hold it.

## Value

`gr(g, dist = "student")` is a formula term, not a free-standing
function: [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md)
reads it at parse time, and the value it contributes is the heavy-tailed
random-effect block of the model, reachable through
[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md) and
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md).
This page itself documents the term grammar and returns nothing.

## Parameterization

A level's coefficients are drawn as `b_j = sqrt(nu * u_j) W z_j` with
`u_j ~ inv-chi2(nu)` and `z_j ~ N(0, I)`, which is the multivariate t
with `nu` degrees of freedom and scale matrix `Sigma = W W'`. The mixing
variable `u_j` is per LEVEL and shared across that level's coefficients,
so a correlated or a [`diag()`](https://rdrr.io/r/base/diag.html) block
is one multivariate t and not several independent univariate ones. brms
builds it the same way.

**`Sigma` is the SCALE, not the covariance.** The variance is
`Sigma * nu / (nu - 2)`, so a standard deviation is
`scale * sqrt(nu / (nu - 2))`.
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
stores the scale matrix, tags it with `nu`, and prints both columns;
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
and `frm_simulate(newparams = )` speak of it as `sd_<group>__<term>`,
which is the name brms gives the same quantity. The correlations are the
same either way.

## Why `nu` is fixed

brms estimates `nu` under a `gamma(2, 0.1)` prior truncated at 1, and
that prior is carrying the parameter. Maximum likelihood has no such
help: profiling a simulated fit with 20 groups leaves the whole grid
from 2.1 to 500 inside the 95% profile interval, and joint ML sends `nu`
to a boundary in 24% to 41% of replicates at 20 groups and still 5% to
14% at 100 (`dev/tre-feasibility.md`). So `nu` is a constant here, set
by `dist_nu` and defaulting to 5. That is the frequentist analogue of
brms's own `prior(constant(3), class = "df")`. To ask what the data say
about it, fit two or three values and compare
[`logLik()`](https://rdrr.io/r/stats/logLik.html).

`dist_nu` must exceed 2, so that the latent has a variance for
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md) to
convert and for new-level prediction to use.

## Accuracy

The t density is not log-concave, so the Laplace approximation is not
exact over it as it is over a gaussian latent. Measured against adaptive
quadrature (`dev/tre-feasibility.md`), the approximation pushes the
estimated latent SCALE UPWARD, and how much depends on how much the data
say about each level:

- With the latent scale near the residual SD, the bias is under 2% of
  one standard error at 8 observations per group and 0.2% at 25, at
  every `nu` from 2.5 up. It reaches 12% of a standard error in the
  worst case tested, `nu = 2.5` with 3 observations per group.

- It becomes material only where the variance component is small AND the
  groups are tiny: at 2 observations per group and a true scale a
  quarter of the residual SD, the Laplace estimate of the scale came
  back three times the exact one.

Two checks, both already in the package. `quadrature = TRUE`
marginalizes a scalar random intercept by Gauss-Kronrod quadrature
instead, which over a t latent is EXACT, not merely better; it is the
one to run when a t block's variance component matters and the groups
are small.
[`frmtmb.sample::check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.html)
measures the same thing without refitting, by NUTS on the objective: its
`z_shift` for `theta` reproduced the displacement above to within a
percentage point in the probe.

A last consequence of the constant:
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) carry roughly `G * c(nu)`
where `G` is the number of levels and
`c(nu) = lgamma((nu+1)/2) - lgamma(nu/2) + log(2/(nu+1))/2` (-0.226 at
`nu = 3`, -0.141 at `nu = 5`). Comparing two t fits with the same `nu`
and the same grouping is fine, because the offset cancels; comparing a t
fit against a gaussian one by AIC is not.

## What is refused

- `gr(cov = )` / `gr(prec = )`:

  A relationship matrix correlates the LEVELS, and the t's mixing
  variable is per level, so the joint density over the field is not a
  multivariate t and has no closed form to hand the Laplace machinery.
  brms writes this combination, but as a hierarchical construction Stan
  samples rather than a density.

- Other covariance structures:

  `us` (the default) and `diag` only. `ar1()`, `cs()`, `toep()` and the
  rest describe a covariance over the block's LEVELS, which the
  per-level mixing variable does not compose with.

- `mm()`:

  A multi-membership row loads several levels at once, so the per-level
  mixing variable has no single value on it.

- `|ID|` keys:

  Merged blocks are assembled as gaussian ones. Write the merged
  coefficients as one term instead -
  `(x1 + x2 | gr(g, dist = "student"))` is the same multivariate-t
  block.

## Downstream

[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md) and
its conditional variances, `sdreport()` standard errors, `REML = TRUE`
and
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
all work: the Laplace machinery does not care which density the latent
has. [`simulate()`](https://rdrr.io/r/stats/simulate.html) and
[`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)
draw a multivariate t with one chi-square per level.
`predict(allow_new_levels = TRUE)` inflates the unseen level's variance
by `nu / (nu - 2)`, so the interval has the right variance around a
heavier-tailed truth. It is still built as a gaussian interval, so it is
not the right quantile far into the tail. Use
[`simulate()`](https://rdrr.io/r/stats/simulate.html) for that.

## See also

[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
for the scale matrix,
[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
for what a `dist = "student"` block may be combined with, and
[`vignette("brms-migration")`](https://aforren1.github.io/frmtmb/articles/brms-migration.md).

## Examples

``` r
set.seed(1)
n <- 12
d <- data.frame(x = rnorm(20 * n), g = factor(rep(1:20, each = n)))
b <- rnorm(20)
b[20] <- b[20] + 6          # one outlying group
d$y <- 1 + 0.5 * d$x + b[d$g] + rnorm(20 * n)

fit_t <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))),
             family = gaussian(), data = d)
fit_n <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)

# the gaussian latent has to widen to cover the outlying group
VarCorr(fit_t)
#>   1 | g 
#>         Name   Scale Std.Dev.
#>  (Intercept) 0.84321   1.0886
#>    Student-t latent, nu = 5 (fixed); the stored matrix is the scale
VarCorr(fit_n)
#>   1 | g 
#>         Name Std.Dev.
#>  (Intercept)   1.5007

# heavier tails, at the cost of a fixed nu
frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 3))),
    family = gaussian(), data = d)
#> frmtmb fit: y ~ x + (1 | gr(g, dist = "student", dist_nu = 3)) 
#> Family: gaussian   Method: ML 
#> logLik: -384.804  AIC: 777.607  nobs: 240 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)           x 
#>      1.2809      0.6492 
#>  sigma:
#> (Intercept) 
#>     0.07106 
#> 
#> Random effects:
#>   1 | g 
#>         Name   Scale Std.Dev.
#>  (Intercept) 0.66974     1.16
#>    Student-t latent, nu = 3 (fixed); the stored matrix is the scale
```
