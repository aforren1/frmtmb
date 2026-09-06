# Per-trial value estimates, prediction errors and choice probabilities

The quantities these papers plot. A learning model is fitted for its
latent trajectory as much as for its parameters, and the trajectory is
not in the fitted object: it is recomputed by replaying the recursion at
the estimates.

## Usage

``` r
frm_value_trace(fit)
```

## Arguments

- fit:

  A fit whose family came from this package.

## Value

A data frame with one row per row of the model frame, in the data's own
order.

## Details

The columns after `trial` are the family's own. Every family returns
`p`, the fitted probability of the choice that was actually made, and
`pe`, the prediction error the outcome produced. The value stores vary:
`q1` and `q2` for a two-armed family, `mu1` to `mu4` and `s1` to `s4`
for the Kalman filter, `qmf1` and `qmb1` for the two-step learner's
model-free and model-based stage-one values. Each value is the one the
choice on that trial was made ON, before the outcome moved it, which is
the ordering a plot of learning needs.

## What `p` is, and what it is not

`p` is the per-trial factor of the likelihood, so `sum(log(p))` is the
CONDITIONAL data log-likelihood given the fitted subject effects. In a
fit with no random effects that is `logLik(fit)` exactly, to machine
precision, and the test suite asserts it.

In a HIERARCHICAL fit the two differ, and by more than rounding:
[`logLik()`](https://rdrr.io/r/stats/logLik.html) reports the
Laplace-approximated MARGINAL likelihood, which adds the random-effect
density and the curvature term that a conditional quantity does not
carry. Measured on a 30-subject fit with `sd(id) = 1.07`, `sum(log(p))`
is -1507.8 against a [`logLik()`](https://rdrr.io/r/stats/logLik.html)
of -1542.5. Neither is wrong; they are different quantities, and only
the conditional one exists per trial.

This function REPLACES
[`stats::fitted()`](https://rdrr.io/r/stats/fitted.values.html) for
these families rather than supplementing it. They decline to declare a
mean on the response scale, because the response is a nominal option
code and `y - mu` would be arithmetic on a category, so
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
`predict(type = "response")` refuse. Everything a mean would have
carried is in this table, with more beside it.

## See also

[`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md),
[`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md)

## Examples

``` r
d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 30,
                     seed = 3)
d$choice <- frm_task_simulate(
  bandit2arm_delta(subject = id, trial = trial), d,
  pars = list(alpha = 0.4, tau = 3), seed = 3)[[1]]$choice
fit <- frmtmb::frm(frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
                   family = bandit2arm_delta(subject = id,
                                             trial = trial), data = d)
tr <- frm_value_trace(fit)
head(tr)
#>   subject trial        q1 q2         pe         p
#> 1       1     1 0.0000000  0  1.0000000 0.5000000
#> 2       1     2 0.5155390  0 -0.5155390 0.8426762
#> 3       1     3 0.2497585  0  0.7502415 0.6927611
#> 4       1     4 0.6365373  0  0.0000000 0.1118311
#> 5       1     5 0.6365373  0  0.3634627 0.8881689
#> 6       1     6 0.8239165  0  0.1760835 0.9359659
# the identity that ties the trace to the likelihood
c(from_trace = sum(log(tr$p)), logLik = as.numeric(stats::logLik(fit)))
#> from_trace     logLik 
#>   -81.9358   -81.9358 
```
