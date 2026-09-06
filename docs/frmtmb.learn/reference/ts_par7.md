# The two-step task: a model-based and model-free hybrid

The task of Daw and others (2011), and the model that made it famous. A
trial has two choices. The first leads, usually but not always, to one
of two second-stage states; the second is a choice within that state,
and it pays. A purely model-free learner repeats what was rewarded. A
model-based learner knows the transition structure and works backward
through it, so a reward that arrived after a RARE transition makes it
switch rather than stay. `w` is how much of each.

## Usage

``` r
ts_par7(subject, trial = NULL, p_common = 0.7)
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

- p_common:

  The probability that a first-stage choice leads to its common
  second-stage state. 0.7 in the original task. Known, not estimated:
  the subject's BELIEF about it is not separable from `w` on the data
  these studies collect.

## Value

A `frmtmb_family` object, for `frm(family = )`.

## Details

    QMB[a] <- P(common) * max(Q2[state a leads to]) +
              P(rare)   * max(Q2[the other state])
    Q1[a]  <- w * QMB[a] + (1 - w) * QMF[a]
    P(a)   <- softmax(tau1 * (Q1 + pers * repeated))
    # after the outcome
    delta1 <- Q2[s2, a2] - QMF[a1];  QMF[a1] <- QMF[a1] + alpha1 * delta1
    delta2 <- reward   - Q2[s2, a2]; Q2[s2, a2] <- Q2[s2, a2] + alpha2 * delta2
    QMF[a1] <- QMF[a1] + alpha1 * lambda * delta2

`lambda` is the eligibility trace: how much of the second-stage
prediction error reaches the first-stage value directly. `pers` is
perseveration, a constant bonus for repeating the previous trial's
first-stage choice. hBayesDM calls the seven-parameter version `ts_par7`
and spells them `a1`, `beta1`, `a2`, `beta2`, `lambda`, `w` and `pi`.

## Two ways this differs from hBayesDM's `ts_par7`

Five of the seven parameters are a plain rename. Two things are not, and
both are read off hBayesDM 2.0.0's published Stan source.

**`pers` does not port, because perseveration sits on the other side of
the temperature.** hBayesDM adds it OUTSIDE `beta1`:
`inv_logit(beta1 * (v2 - v1) + pi * (...))`. This family adds it INSIDE
the bracket `tau1` multiplies, so its contribution to the same logit is
`tau1 * pers * (...)`. Matching term by term,

    pers = pi_hBayesDM / tau1

`tau1` is estimated, so as with
[`prl_fictitious()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/prl_fictitious.md)'s
`bias` there is no constant transform: convert through a fit's own
`tau1`, or refit.

**The eligibility trace uses a different prediction error, and this
family follows the paper.** The trace term here is formed from the
PRE-update `delta2`, which is the equation Daw and others (2011) print.
hBayesDM's code updates the stage-two value first and then forms the
trace from the already-updated value, so its effective trace is
`lambda * a1 * (1 - a2) * delta2`. The two agree only at `a2 = 0`. This
is a difference in the reference implementation rather than an error
here, and it is recorded because a user comparing fits will otherwise
meet it with no stated cause.

## The data a trial carries

The response is the FIRST-stage choice, 1 or 2. `stage2(state, choice)`
carries the second-stage state (1 or 2, for the two states) and the
choice made in it (1 or 2). `payoff(p1, p2, p3, p4)` carries what each
of the four second-stage options would have paid, indexed
`2 * (state - 1) + choice`.

The transition probability is a known constant of the TASK rather than a
parameter, so `P(state | first choice)` contributes a constant to the
log likelihood and is dropped. It still enters the model-based value,
which is what `p_common` is for.

## Why [`simulate()`](https://rdrr.io/r/stats/simulate.html) is refused

One trial's draw is three numbers: the first-stage choice, the state the
environment answers with, and the second-stage choice.
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) and
[`frmtmb::frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.html)
return one response vector, and there is nowhere to put the other two,
so both refuse by name.
[`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md)
returns whole data frames and does the draw properly; it is also what
every other family in this package uses for a recovery study.

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

Daw, N. D., Gershman, S. J., Seymour, B., Dayan, P. and Dolan, R. J.
(2011). Model-based influences on humans' choices and striatal
prediction errors. *Neuron* 69, 1204-1215.

## See also

[`frm_task_design()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_design.md),
[`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md)

## Examples

``` r
d <- frm_task_design("twostep", n_subject = 5, n_trial = 40, seed = 8)
d <- frm_task_simulate(
  ts_par7(subject = id, trial = trial), d,
  pars = list(alpha1 = 0.4, tau1 = 3, alpha2 = 0.4, tau2 = 3,
              lambda = 0.6, w = 0.5, pers = 0.2), seed = 8)[[1]]
head(d[, c("id", "trial", "choice", "state2", "choice2")])
#>   id trial choice state2 choice2
#> 1  1     1      1      2       1
#> 2  1     2      2      2       2
#> 3  1     3      1      1       2
#> 4  1     4      1      1       2
#> 5  1     5      2      2       1
#> 6  1     6      2      2       1
```
