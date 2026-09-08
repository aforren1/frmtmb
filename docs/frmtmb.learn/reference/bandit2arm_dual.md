# Two-armed delta learning with separate rates for gains and losses

[`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md)
with one learning rate replaced by two, so that good news and bad news
are absorbed at different speeds. This asymmetry is the parameter most
of the clinical literature is about, and it is the model hBayesDM calls
`prl_rp` when the split is on the sign of the OUTCOME.

## Usage

``` r
bandit2arm_dual(subject, trial = NULL, split = c("pe", "outcome"))
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

- split:

  `"pe"` splits on the sign of the prediction error, `"outcome"` on the
  sign of the outcome itself.

## Value

A `frmtmb_family` object, for `frm(family = )`.

## Details

    pe  <- reward - Q[chosen]
    Q[chosen] <- Q[chosen] + (pe > 0 ? Arew : Apun) * pe

## Which sign, and why it matters here more than elsewhere

Two models are in circulation and they are not the same model. One
splits on the sign of the PREDICTION ERROR: an outcome better than
expected is a gain even when it is a small reward. The other splits on
the sign of the OUTCOME, which is what `prl_rp` does. `split` chooses,
and the default is `"pe"`.

They also cost different things to tape. The outcome's sign is DATA, so
`split = "outcome"` selects the rate with a column of zeros and ones
computed once and adds nothing to the tape. A prediction error is an AD
quantity, and RTMB refuses a comparison on one outright, so
`split = "pe"` selects with
[`sign()`](https://rdrr.io/r/base/sign.html). That is exact rather than
smoothed, and the kink it introduces at `pe == 0` is the model's own:
the update is zero there under both rates, so the two branches agree at
the crossing and the derivative from each side is the right one.

With `split = "outcome"`, an outcome of exactly zero counts as
punishment, which is what `prl_rp` does and what the usual `0`/`1`
payoff coding needs. Code punishment as `-1` if you want the symmetric
reading.

One consequence of the kink is worth knowing before fitting `"pe"`
hierarchically on GRADED payoffs. The joint log density is then
non-smooth in the random effects wherever a prediction error crosses
zero, and the Laplace inner solve stops short of the conditional mode:
measured on the identity fixture, the largest gradient on the subject
effects is 4.0e-02 under `"pe"` against 1.4e-15 under `"outcome"` on the
same data. The fit and its log likelihood are unaffected to machine
precision, but a convergence warning on a `"pe"` fit with graded payoffs
is expected rather than alarming. On binary payoffs neither split has a
live kink, and the two are the same model anyway.

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

[`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md),
[`frm_learn_families()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_learn_families.md)

## Examples

``` r
d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 40, seed = 2)
d$choice <- frm_task_simulate(
  bandit2arm_dual(subject = id, trial = trial), d,
  pars = list(Arew = 0.5, Apun = 0.15, tau = 3),
  seed = 2)[[1]]$choice
fit <- frmtmb::frm(
  frmtmb::bf(choice | reward(pay1, pay2) ~ 1, Apun ~ 1, tau ~ 1),
  family = bandit2arm_dual(subject = id, trial = trial), data = d)
frmtmb::fixef(fit)
#> $Arew
#> (Intercept) 
#>   0.8827643 
#> 
#> $Apun
#> (Intercept) 
#>   -1.837557 
#> 
#> $tau
#> (Intercept) 
#>   0.9807577 
#> 
```
