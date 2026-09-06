# Counterfactual (fictitious) updating on two options

The unchosen option is updated too, in the opposite direction. In a
two-armed task where one option pays when the other does not, that is
what a subject who understands the task should do, and it is the model
of Glascher, Hampton and O'Doherty (2009) for probabilistic reversal
learning. hBayesDM calls it `prl_fictitious`.

## Usage

``` r
prl_fictitious(subject, trial = NULL)
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

## Value

A `frmtmb_family` object, for `frm(family = )`.

## Details

    EV[chosen]   <- EV[chosen]   + alpha * ( outcome - EV[chosen])
    EV[unchosen] <- EV[unchosen] + alpha * (-outcome - EV[unchosen])
    P(option 1)  <- plogis(tau * (EV1 - EV2) + bias)

The counterfactual outcome is the NEGATIVE of the realized one, which is
the model's assumption of an anticorrelated task rather than a reading
of the second `reward()` column. The second column is still what the
simulator pays a counterfactual choice with.

## What `bias` is, and why it does not port

`bias` is the indecision point: a constant added to the utility
difference, so a subject with `bias` above zero prefers option 1 at
equal value. It has an identity link.

**It is not hBayesDM's `alpha` renamed, and there is no fixed transform
between them.** hBayesDM's `prl_fictitious` writes the choice
probability as `inv_logit(beta * (alpha - (ev1 - ev2)))`, which is
`inv_logit(beta * (ev1 - ev2) - beta * alpha)` with the option labels
swapped. This family writes `inv_logit(tau * (ev1 - ev2) + bias)`.
Matching the two term by term gives

    bias = -tau * alpha_hBayesDM

The factor is `tau`, which is ESTIMATED, so the map depends on the fit
and no constant relates the two parameters. At the `tau = 3` this page's
example uses, an indecision point carried across unchanged is wrong by a
factor of three, and the resulting model converges quietly. Convert
through the equation above, using the same fit's `tau`, or refit.

What is implemented here is the equation printed at the top of this
page, and it is that equation the Stan identity in
`tests/testthat/test-stan-identity.R` checks to 8.5e-14. The relation
above is read off hBayesDM 2.0.0's published Stan source rather than
measured against a fit: `dev/hbayesdm-crosscheck.R` would measure it and
has not been run, for the toolchain reason recorded in its header.

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

Glascher, J., Hampton, A. N. and O'Doherty, J. P. (2009). Determining a
role for ventromedial prefrontal cortex in encoding action-based value
signals during reward-related decision making. *Cerebral Cortex* 19,
483-495.

## See also

[`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md),
[`bandit2arm_dual()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_dual.md)

## Examples

``` r
d <- frm_task_design("reversal", n_subject = 6, n_trial = 40, seed = 4)
# a reversal task pays -1 as often as +1
d$pay1 <- 2 * d$pay1 - 1
d$pay2 <- 2 * d$pay2 - 1
d$choice <- frm_task_simulate(
  prl_fictitious(subject = id, trial = trial), d,
  pars = list(alpha = 0.3, bias = 0, tau = 3), seed = 4)[[1]]$choice
fit <- frmtmb::frm(
  frmtmb::bf(choice | reward(pay1, pay2) ~ 1, bias ~ 1, tau ~ 1),
  family = prl_fictitious(subject = id, trial = trial), data = d)
frmtmb::fixef(fit)
#> $alpha
#> (Intercept) 
#>  -0.8286018 
#> 
#> $bias
#> (Intercept) 
#> -0.08634232 
#> 
#> $tau
#> (Intercept) 
#>   0.8673196 
#> 
```
