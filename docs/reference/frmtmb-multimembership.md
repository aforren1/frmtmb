# Multi-membership random effects

`(x | mm(g1, g2, ...))` says that one observation belongs to SEVERAL
levels of one grouping factor at once: a pupil taught in more than one
school, a fish caught in more than one net, a paper written by more than
one author. `mm()` is written where the grouping factor of a bar term
goes, and `mmc()` supplies the member-specific covariate of a random
slope over it. Both spellings follow brms.

## Value

`mm()` and `mmc()` are formula terms, not free-standing functions:
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) reads them
at parse time, and the value they contribute is the shared
multi-membership random-effect block of the model, reachable through
[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md) and
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md).
This page itself documents the term grammar and returns nothing.

## What mm() changes

Only the random-effect design matrix. The membership variables are
pooled into ONE grouping factor whose levels are every level named by
any member, and the block over them is an ordinary `us` (or `diag`)
block, exactly the block `(x | g)` would build. What differs is the
design row: instead of putting a 1 in one level's column it puts weight
`w_j` in each member level's column, so the effect the observation sees
is the weighted average of its members' effects.

Everything downstream is therefore unchanged. The covariance is the
usual one, the Laplace approximation is the usual one, and
[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md),
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md),
[`ngrps()`](https://aforren1.github.io/frmtmb/reference/ngrps.md),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) all read the block
with no multi-membership branch.

## Levels are pooled

The pooled level set is each membership variable's own levels,
concatenated in the order the variables were written and then
deduplicated - brms builds it the same way. So `mm(g1, g2)` and
`mm(g2, g1)` fit the same model but order the coefficients differently,
and a level that only one of the two variables carries still gets its
own coefficient. Because the members share one level set, the same label
in `g1` and in `g2` means the same school, which is the point of the
term; relabel one of them if it does not.

## Weights

`weights` is a matrix with one row per observation and one column per
membership variable, usually built with
[`cbind()`](https://rdrr.io/r/base/cbind.html). Without it every member
gets `1/J`, where `J` is the number of membership variables, and that
default is NOT rescaled. With it, `scale = TRUE` (the default) divides
each row by its sum, so the weights become proportions; `scale = FALSE`
uses the numbers as they are, negative ones included. Scaling refuses
negative weights and rows that sum to zero, because neither has a
proportion to be normalized to.

## Member-specific covariates

`mmc(x1, x2)` is ONE random-slope coefficient whose covariate value
differs by member: member 1 uses `x1`, member 2 uses `x2`. It takes one
variable per membership variable, they must be numeric, and it has to be
a term of its own on the left of the bar. A plain column on the left of
the bar - `(1 + z | mm(g1, g2))` - is the other case: one slope whose
covariate is the same for every member.

## What is refused

- Other covariance structures:

  `mm()` carries `us` and `diag` only. `ar1()`, `cs()`, `toep()`,
  [`exp()`](https://rdrr.io/r/base/Log.html), `gr(cov = )` and the rest
  all describe a covariance over the block's levels, and the pooled
  membership levels have no order, no coordinates and no relationship
  matrix for one to be defined on.

- `|ID|` keys:

  A merged block indexes one level set per observation row; an `mm()`
  row loads several at once.

- brms's other `mm()` arguments:

  `cor = FALSE` is `diag(x | mm(g1, g2))`, `id =` is the `|ID|` key, and
  `cov =` is `gr(g, cov = A)` over a single-membership factor, and
  `dist =` is `gr(g, dist = "student")` over one (see
  [frmtmb-student-re](https://aforren1.github.io/frmtmb/reference/frmtmb-student-re.md)),
  though not over a membership design, whose rows load several levels at
  once. `pw =` has no equivalent yet.

- Non-name members:

  `mm()` reads its membership variables as column names, as brms does.
  Build the column first.

On `newdata`, a membership level that was not in the fitted data needs
`allow_new_levels = TRUE`; that member then contributes the population
value while the row's remaining members still contribute their fitted
effects.

## One covariance per by-level

`mm(g1, g2, by = cbind(f1, f2))` is brms's by-split: `by` is a matrix
with one column per membership variable, and each pooled level gets the
covariance of its by-level. A pooled level must have one by-level
wherever it appears, in any member column; brms refuses the data
otherwise, and so does frmtmb, with brms's message. The block splits
into one block per by-level, and the standard deviations take brms's
names, `sd_mmg1g2__Intercept:cbind(f1, f2)1` (brms evaluates
[`cbind()`](https://rdrr.io/r/base/cbind.html) of two factors to their
integer codes). A new level on `newdata` takes the covariance of the
by-level its own member column names. The same split over a single
grouping factor is `gr(g, by = f)`.

## See also

[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md) and
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
for the fitted block,
[`vignette("brms-migration")`](https://aforren1.github.io/frmtmb/articles/brms-migration.md)
for the porting notes, and
[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
for what `mm()` may be combined with.

## Examples

``` r
set.seed(1)
n <- 200
d <- data.frame(
  x = rnorm(n),
  school1 = factor(sample(letters[1:8], n, TRUE)),
  school2 = factor(sample(letters[5:12], n, TRUE)),
  share1 = runif(n, 0.5, 1)
)
d$share2 <- 1 - d$share1
u <- rnorm(12, 0, 0.8)
names(u) <- letters[1:12]
d$y <- 1 + 0.5 * d$x +
  0.5 * u[as.character(d$school1)] +
  0.5 * u[as.character(d$school2)] + rnorm(n, 0, 0.5)

# equal membership: each pupil is half of each school
fit <- frm(bf(y ~ x + (1 | mm(school1, school2))) + gaussian(),
           data = d)
summary(fit)
#>  Family: gaussian 
#>  Links: mu = identity; sigma = log
#> 
#> Formula: y ~ x + (1 | mm(school1, school2)) 
#>    Data: d (Number of observations: 200) 
#>  Method: ML   logLik: -172.632   AIC: 353.264   BIC: 366.457 
#> 
#> Multilevel Hyperparameters:
#> ~mmschool1school2 (Number of levels: NA) 
#>               Estimate Est.Error l-95% CI u-95% CI
#> sd(Intercept)     0.92       0.2     0.61     1.41
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI z value Pr(>|z|)
#> Intercept     0.70      0.27     0.17     1.23    2.58   0.0099
#> x             0.49      0.04     0.41     0.58   11.94   <2e-16
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI
#> sigma     0.52      0.03     0.47     0.58
# one coefficient per pooled school level
ranef(fit)
#> $mm(school1, school2)   (1 | mm(school1, school2))
#>     Intercept
#> a  0.18370669
#> b -0.02184145
#> c -0.74232816
#> d  0.46518610
#> e  0.86556173
#> f  1.45758646
#> g -0.52474894
#> h  0.07022447
#> i  1.24095653
#> j -0.60886148
#> k -1.93404777
#> l -0.45139622
#> 

# the time each pupil spent in each school, as proportions
frm(bf(y ~ x + (1 | mm(school1, school2,
                       weights = cbind(share1, share2)))) + gaussian(),
    data = d)
#>  Family: gaussian 
#>  Links: mu = identity; sigma = log
#> 
#> Formula: y ~ x + (1 | mm(school1, school2, weights = cbind(share1, share2))) 
#>    Data: d (Number of observations: 200) 
#>  Method: ML   logLik: -208.589   AIC: 425.179   BIC: 438.372 
#> 
#> Multilevel Hyperparameters:
#> ~mmschool1school2weightsEQcbindshare1share2 (Number of levels: NA) 
#>               Estimate Est.Error l-95% CI u-95% CI
#> sd(Intercept)     0.96      0.23      0.6     1.52
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI z value Pr(>|z|)
#> Intercept     0.66      0.28     0.10     1.21    2.30    0.021
#> x             0.47      0.05     0.38     0.57    9.39   <2e-16
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI
#> sigma     0.63      0.03     0.57     0.69
```
