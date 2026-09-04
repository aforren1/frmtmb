# Prior objects, addressed by internal parameter name

The named-list prior spelling, as opposed to
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)'s
classes. Priors written this way apply on the INTERNAL parameter scale:
coefficients are on their link scale, and covariance parameters
(`theta_*`) are the unconstrained parameterization (log-SDs,
scaled-Cholesky terms), so `prior_normal(0, 1)` on `theta_1` is a
lognormal prior on that standard deviation.

## Usage

``` r
prior_normal(location = 0, scale = 1)

prior_t(df = 3, location = 0, scale = 1)

prior_lkj(eta = 1)
```

## Arguments

- location, scale, df:

  Prior parameters.

- eta:

  LKJ shape. `1` is uniform over correlation matrices, larger values
  concentrate toward the identity, and `0 < eta < 1` pushes toward the
  boundary.

## Value

A `frmtmb_prior` object.

## Details

[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) takes them
as a MAP penalty, and so does
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html),
where they take over exactly the parameters they name and leave the rest
of the prior stack in place.
[`par_template()`](https://aforren1.github.io/frmtmb/reference/par_template.md)
and
[`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.md)
name the addressable slots.

## See also

[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
for the class-based spelling, which is the one most models want.

## Examples

``` r
# the objects themselves are cheap descriptions
prior_normal(0, 2)
#> $kind
#> [1] "normal"
#> 
#> $location
#> [1] 0
#> 
#> $scale
#> [1] 2
#> 
#> attr(,"class")
#> [1] "frmtmb_prior"
prior_t(df = 3, location = 0, scale = 1)
#> $kind
#> [1] "t"
#> 
#> $df
#> [1] 3
#> 
#> $location
#> [1] 0
#> 
#> $scale
#> [1] 1
#> 
#> attr(,"class")
#> [1] "frmtmb_prior"

set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)

# names are internal parameter names, or whole components. theta_1
# is a log-SD, so a normal there is a lognormal on the SD.
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
           prior = list(beta = prior_normal(0, 5),
                        theta_1 = prior_t(3, 0, 1)))
prior_summary(fit)
#> $beta
#> $kind
#> [1] "normal"
#> 
#> $location
#> [1] 0
#> 
#> $scale
#> [1] 5
#> 
#> attr(,"class")
#> [1] "frmtmb_prior"
#> 
#> $theta_1
#> $kind
#> [1] "t"
#> 
#> $df
#> [1] 3
#> 
#> $location
#> [1] 0
#> 
#> $scale
#> [1] 1
#> 
#> attr(,"class")
#> [1] "frmtmb_prior"
#> 
```
