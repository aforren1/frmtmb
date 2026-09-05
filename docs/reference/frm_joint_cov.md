# The joint covariance of the fixed and random coefficients

The covariance of everything the fit estimates, `beta`, `betad`, `theta`
and the random-effect coefficients `b` together, in one matrix. It is
what a delta method over a fitted CURVE needs and what
[`vcov()`](https://rdrr.io/r/stats/vcov.html) cannot return: a penalized
smooth's wiggly part is a random-effect block even when the smooth is a
population term, so a covariance that stops at the fixed effects covers
none of it.

## Usage

``` r
frm_joint_cov(object)
```

## Arguments

- object:

  A fitted `frmtmb_fit`.

## Value

A list with

- `V`:

  the `p x p` joint covariance.

- `names`:

  length `p`; the PARAMETER COMPONENT each row belongs to (`"beta"`,
  `"betad"`, `"b"`, `"theta"`, ...), which is what
  [`frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.md)`$coef_pos`
  indexes.

- `labels`:

  length `p`; one label per row, `beta.<coefficient>` for a fixed effect
  and `b.<block>.<level>` for a random one.

A fit with no random effects has no joint precision, and `V` is then the
fixed-effect covariance `sdreport()$cov.fixed`; `names` says so.

## Details

`vcov(full = TRUE)` returns the OUTER parameter vector's covariance and
its row names are documented to be exactly
[`confint()`](https://rdrr.io/r/stats/confint.html)'s, so `b` is not in
it and cannot be added.
[`frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.md)
is the accessor for a design at `newdata`; this is the accessor for the
covariance those coefficients have.

The result is memoized on the fit, so the joint-precision solve is paid
once however many curves are drawn from it. It is also the ONLY route to
the covariance of an AUTOSCALED fit
(`frmtmb_control(autoscale = TRUE)`): a fresh
`RTMB::sdreport(getJointPrecision = TRUE)` goes round the
reparameterization and returns a covariance built on the unscaled
Hessian, which is a different matrix and is not marked as one.

## See also

[`frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.md),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, each = 10)))
dd$y <- rnorm(120, 1 + 2 * dd$x + rnorm(12, 0, 0.5)[dd$g], 0.4)
fit <- frm(bf(y ~ x + (1 | g)), data = dd)
jc <- frm_joint_cov(fit)
dim(jc$V)
#> [1] 16 16
table(jc$names)
#> 
#>     b  beta betad theta 
#>    12     2     1     1 
head(jc$labels)
#> [1] "beta.(Intercept)"        "beta.x"                 
#> [3] "betad.sigma_(Intercept)" "b.1 | g.1.(Intercept)"  
#> [5] "b.1 | g.2.(Intercept)"   "b.1 | g.3.(Intercept)"  
```
