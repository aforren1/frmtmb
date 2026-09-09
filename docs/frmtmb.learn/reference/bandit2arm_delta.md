# Rescorla-Wagner delta learning on a two-armed bandit

The simplest value-learning model, and the one every other family here
is a variation on. A subject keeps a value estimate `Q` for each arm,
starts both at zero, and after each choice moves the chosen arm's
estimate toward what it just paid:

## Usage

``` r
bandit2arm_delta(subject, trial = NULL)
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

    Q[chosen] <- Q[chosen] + alpha * (reward - Q[chosen])

The choice is a softmax over the two values, which for two arms is a
logistic function of their difference:

    P(arm 1) <- plogis(tau * (Q1 - Q2))

This is the model hBayesDM calls `bandit2arm_delta`. Its two parameters
are named `A` and `tau` there and `alpha` and `tau` here; see
[`frm_learn_families()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_learn_families.md)
for the whole naming map.

## The data a trial carries

The response is the arm chosen, coded 1 or 2. `reward(pay1, pay2)`
carries what each arm WOULD have paid on this trial, in arm order. The
likelihood only ever reads the chosen arm's entry, so data that records
the received outcome alone can pass it twice: the fit is then EXACTLY
right rather than nearly so. The unchosen entry is multiplied by a zero
indicator before it reaches anything, so the objective is bitwise the
same function of the parameters whatever that column holds.

A DRAW is a different question. A simulated subject chooses for itself,
so paying it needs what the arm the real subject did not take would have
given, which a duplicated column does not record. With the two columns
identical on every trial the drawn choices carry no learning signal at
all, so [`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.html)
and
[`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md)
refuse that data by name rather than returning a draw from a task nobody
ran. Some rows equal is not the signature; every row is. `newdata` is
not a way round it, because the formula names one column twice and
`newdata` is read through that same formula; refit with a column per
arm.

This is true of every family in the package, not only this one. All
eight read the CHOSEN option's entry alone, which is measured rather
than asserted: replacing the unchosen entries with noise leaves the
log-likelihood bitwise unchanged in each of them.

## Reversal learning, and what the grammar buys

Nothing here is specific to a stationary bandit. A probabilistic
reversal-learning task is this family with a covariate on the learning
rate, because `alpha` is an ordinary distributional parameter with its
own linear predictor:

    frm(bf(choice | reward(pay1, pay2) ~ after_reversal + (1 | id),
           tau ~ 1 + (1 | id)),
        family = bandit2arm_delta(subject = id, trial = trial), data = d)

A separate `prl` family would fit one learning rate before the reversal
and one after. This fits the DIFFERENCE, with a standard error, and it
takes a factor with any number of levels, a smooth term or a random
slope in the same place.
[`vignette("learning")`](https://aforren1.github.io/frmtmb/frmtmb.learn/articles/learning.md)
works the example through.

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

Ahn, W.-Y., Haines, N. and Zhang, L. (2017). Revealing
neurocomputational mechanisms of reinforcement learning and
decision-making with the hBayesDM package. *Computational Psychiatry* 1,
24-57.

## See also

[`bandit2arm_dual()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_dual.md)
for separate rates on gains and losses,
[`prl_fictitious()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/prl_fictitious.md)
for counterfactual updating,
[`frm_value_trace()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_value_trace.md)
for the fitted values and prediction errors.

## Examples

``` r
d <- frm_task_design("bandit2arm", n_subject = 8, n_trial = 40,
                     seed = 1)
d$choice <- frm_task_simulate(
  bandit2arm_delta(subject = id, trial = trial), d,
  pars = list(alpha = 0.4, tau = 3), seed = 1)[[1]]$choice
fit <- frmtmb::frm(frmtmb::bf(choice | reward(pay1, pay2) ~ 1,
                              tau ~ 1),
                   family = bandit2arm_delta(subject = id,
                                             trial = trial), data = d)
frmtmb::fixef(fit)
#> $alpha
#> (Intercept) 
#>  -0.6021218 
#> 
#> $tau
#> (Intercept) 
#>    1.188803 
#> 
head(frm_value_trace(fit))
#>   subject trial        q1 q2         pe         p
#> 1       1     1 0.0000000  0  1.0000000 0.5000000
#> 2       1     2 0.3538584  0  0.6461416 0.7616541
#> 3       1     3 0.5825010  0  0.4174990 0.8712927
#> 4       1     4 0.7302365  0 -0.7302365 0.9166345
#> 5       1     5 0.4718362  0  0.5281638 0.8247849
#> 6       1     6 0.6587314  0 -0.6587314 0.8968508
```
