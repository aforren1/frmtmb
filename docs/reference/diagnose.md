# Convergence diagnostics for a frmtmb fit

Reports the optimizer's own verdict plus six checks that a converged fit
can still fail: non-finite standard errors, flat directions, complete
separation in a binomial-type fit, a distributional parameter whose
maximum likelihood is outside its own range, predictor columns scaled
far from one, and variance components on the boundary of their parameter
space (lme4's `isSingular()`, read off the estimates rather than the
Hessian).

## Usage

``` r
diagnose(fit, quiet = FALSE)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- quiet:

  If `TRUE`, return the diagnostics without printing.

## Value

Invisibly, a list of diagnostics.

## Details

A DISTRIBUTIONAL PARAMETER AT THE END OF ITS LINK is one whose estimate
is far out on the link scale AND whose standard error is larger than the
estimate itself. The likelihood was still rising where the optimizer
stopped, so the number reported is the stopping point and not an
estimate.
[`student()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)'s
`nu` does this on any data with no heavy tails, because a student-t
reaches the gaussian only as `nu` goes to infinity: the maximum is never
attained, and two runs of the same model can report degrees of freedom
orders of magnitude apart while agreeing on every coefficient. Refit
with [`gaussian()`](https://rdrr.io/r/stats/family.html), or hold `nu`
somewhere finite with
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md).
The same parameter runs the other way, down to one, on data whose tails
are heavier than any identified `nu` can hold. The check names that too,
with a large NEGATIVE estimate and a natural-scale value of one, and
there [`gaussian()`](https://rdrr.io/r/stats/family.html) is the wrong
answer: the data is the message, and a prior is the way to hold `nu`
finite.

A FLAT DIRECTION is an outer parameter the likelihood does not depend
on: zero gradient and an empty Hessian row. It separates the two causes
of `NaN` standard errors. Parameters that trade off against each other
are over-parameterization, and the model is too big. Parameters the
likelihood is flat in are unidentified AT THIS POINT, and the remedy is
a starting value: a nonlinear term evaluated outside its own support (a
bump whose centre starts far from the data) is flat in several of its
parameters at once. One unusable direction makes EVERY standard error
`NaN`, so `bad_se` names the whole vector and `flat` names the cause.
The check is measured by perturbing each candidate and seeing whether
the gradient moves, and runs only when the covariance has already
failed.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
diagnose(fit)
#> Optimizer convergence code: 0 (relative convergence (4)) 
#> Max |gradient|: 0.0002949 at x 
#> Hessian positive definite: TRUE 
#> No convergence problems detected

# a random effect the data cannot support collapses to the boundary,
# which is a valid fit but a warning about the model
dd$h <- factor(rep(1:10, each = 10))
fit_s <- frm(bf(y ~ x + (1 | g) + (1 | h)) + gaussian(), data = dd)
d <- diagnose(fit_s, quiet = TRUE)
d$singular
#> NULL

# a predictor scaled far from one slows the optimizer down; the
# remedy is frmtmb_control(autoscale = TRUE)
dd$xbig <- dd$x * 1e5
diagnose(frm(bf(xbig ~ 1) + gaussian(), data = dd), quiet = TRUE)$scale
#> NULL
```
