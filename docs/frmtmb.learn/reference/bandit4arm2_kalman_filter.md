# A Kalman filter over the restless four-armed bandit

The task of Daw and others (2006): four arms whose payoffs drift, so
that a subject has to keep exploring. A delta rule with a fixed learning
rate is the wrong model for it, because the right learning rate depends
on how uncertain the subject currently is about the arm it just pulled.
The Kalman filter is that model: it carries a mean AND a variance for
each arm, and its gain is the learning rate the uncertainty implies.

## Usage

``` r
bandit4arm2_kalman_filter(subject, trial = NULL, sigma_o = 4)
```

## Arguments

- subject:

  The column separating one learner's trial sequence from the next,
  given unquoted. It is carried by the family rather than by the
  formula, the way
  [`frmtmb.latent::hmm()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.html)
  carries its sequence grouping, because it orders the likelihood rather
  than entering a linear predictor.

- trial:

  The column giving trial order within a subject, given unquoted. `NULL`
  uses the order the rows appear in.

- sigma_o:

  The observation noise standard deviation, held FIXED rather than
  estimated. Daw and others use 4. It sets the scale the gain is
  measured against, and estimating it alongside `sigma0` and `sigmaD`
  asks the data to separate three quantities that enter through two
  ratios.

## Value

A `frmtmb_family` object, for `frm(family = )`.

## Details

    # choice, on the current posterior means
    P(arm k) <- softmax(tau * mu)
    # observation, chosen arm only
    gain     <- v[k] / (v[k] + sigma_o^2)
    mu[k]    <- mu[k] + gain * (payoff - mu[k])
    v[k]     <- v[k] * (1 - gain)
    # diffusion, every arm, ready for the next trial
    mu       <- lambda * mu + (1 - lambda) * center
    v        <- lambda^2 * v + sigmaD^2

hBayesDM calls this `bandit4arm2_kalman_filter` and spells the six
parameters `lambda`, `theta`, `beta`, `mu0`, `sigma0` and `sigmaD`.

## What the variance buys, and what it costs

The gain falls as a subject learns, so early trials move the mean
further than late ones without any parameter saying so. Nothing else in
this package has a learning rate that changes within a subject without a
covariate to change it.

The cost is identifiability. `sigma0`, `sigmaD` and the fixed
observation noise `sigma_o` enter the gain only through ratios, and
`tau` trades off against the scale of the payoffs. Fitting all six
freely on a short session is not advisable; fix `sigma_o` (it is an
argument, not a parameter, for that reason), and consider fixing `mu0`
and `center` at the task's own center when the design is known.

## Exploration bonus

Not implemented. Daw and others compare softmax choice with rules that
add a bonus for uncertainty, and the filter here carries the variance
those rules need, so it is a short addition; it is left out because a
bonus and `tau` are hard to separate on the data sizes these studies
run. The state is in
[`frm_value_trace()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_value_trace.md)
as `s1` to `s4` if you want to look.

## The Laplace caveat

frmtmb integrates the subject effects out with a Laplace approximation,
which is exact only when the conditional log-density is quadratic. For
binary choices it is not, and the fewer trials a subject has the less
quadratic it is. Measured on this family, the fixed effects survive
short sessions and the variance components do not: see the Laplace
section of
[`vignette("learning")`](https://aforren1.github.io/frmtmb/frmtmb.learn/articles/learning.md)
and
[`?frmtmb.learn`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frmtmb.learn-package.md)
for the numbers. `frm(importance =)`, which is the usual way to price
that error, is REFUSED for every family here; the refusal names the
seam.

## References

Daw, N. D., O'Doherty, J. P., Dayan, P., Seymour, B. and Dolan, R. J.
(2006). Cortical substrates for exploratory decisions in humans.
*Nature* 441, 876-879.

## See also

[`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md),
[`frm_task_design()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_design.md)

## Examples

``` r
d <- frm_task_design("bandit4arm_restless", n_subject = 5, n_trial = 40,
                     seed = 6)
d$choice <- frm_task_simulate(
  bandit4arm2_kalman_filter(subject = id, trial = trial), d,
  pars = list(lambda = 0.98, center = 50, tau = 0.15, mu0 = 50,
              sigma0 = 10, sigmaD = 3), seed = 6)[[1]]$choice
table(d$choice)
#> 
#>  1  2  3  4 
#> 46 56 58 40 
```
