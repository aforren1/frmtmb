# A Kalman filter over the restless four-armed bandit

The task of Daw and others (2006): four arms whose payoffs drift, so
that a subject has to keep exploring. A delta rule with a fixed learning
rate is the wrong model for it, because the right learning rate depends
on how uncertain the subject currently is about the arm it just pulled.
The Kalman filter is that model: it carries a mean AND a variance for
each arm, and its gain is the learning rate the uncertainty implies.

## Usage

``` r
bandit4arm2_kalman_filter(subject, trial = NULL, sigma_o = 4, bonus = FALSE)
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

- bonus:

  Add the exploration bonus to the choice rule, giving one more
  parameter `phi`. `FALSE`, the default, is the plain softmax over
  posterior means and is the model the first release fitted.

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

`bonus = TRUE` adds Daw and others' exploration bonus: an arm's utility
becomes

    util[j] <- tau * (mu[j] + phi * sqrt(s[j]))

so an arm the subject is UNCERTAIN about is worth more than its
posterior mean alone, by an amount `phi` the fit estimates. `phi` is an
ordinary distributional parameter, so it takes a formula like any other;
`phi = 0` is the plain softmax and is what `bonus = FALSE`, the default,
fits.

IT IS OFF BY DEFAULT so that every fit made before it existed is
unchanged, and the suite pins that by holding `phi` at zero with
`bf(phi = 0)` and requiring the two log-likelihoods to agree.

WHAT IT TRADES OFF AGAINST IS NOT `tau`, and this help said otherwise
before the study ran. The reasoning was that both control how far choice
departs from the current best arm. Measured at 30 subjects by 100 trials
with `center`, `mu0` and `sigma0` held at the task's own values, 60
replicates: `phi` recovers with a bias of +0.06 against a truth of 1.5,
a spread of 0.14 and coverage 0.93, and its estimates correlate with
`tau`'s at 0.12, which is nothing. What they correlate with is `sigmaD`,
at -0.81, and after the fact that is the obvious pair: `phi` multiplies
the posterior standard deviation and `sigmaD` sets how fast that
standard deviation grows, so the two scale the same term. Report them
together. `tau`'s own interval undercovers a little here, 0.85 against
the nominal 0.95.

Turn the bonus on when uncertainty-driven exploration is the question
the study asks. Leaving `center`, `mu0` and `sigma0` free as well is not
advisable on a session this size, for the reason the section above
gives.

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
for the numbers.

`frm(importance =)`, the usual way to price that error, now WORKS for
every family here. It was refused in the first release because the
correction needs one log-likelihood value per subject and the structured
protocol had no slot to put them in; the slots landed and the families
declare them. What the correction is worth depends on the design and not
on the family: at 40 subjects by 100 trials it moves a well-identified
subject-level standard deviation up by about 0.1 to 0.2 log units,
toward the truth, with good diagnostics, and at 20 trials it has nothing
to correct because the Laplace fit has usually collapsed the component
to zero already.
[`?frmtmb.learn`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frmtmb.learn-package.md)
carries the per-dataset table.

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
