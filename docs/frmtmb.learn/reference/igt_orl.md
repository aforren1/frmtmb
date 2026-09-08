# Outcome-representation learning for the Iowa gambling task

The model Haines, Vassileva and Ahn (2018) built to fix what
[`igt_pvl_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/igt_pvl_delta.md)
cannot see. A prospect-theory learner carries one number per deck, the
expected value, so it cannot distinguish a deck that pays a little often
from one that pays a lot rarely once the two average out. ORL carries
three:

## Usage

``` r
igt_orl(subject, trial = NULL)
```

## Arguments

- subject:

  The column separating one learner's trial sequence from the next,
  given unquoted.

- trial:

  The column giving trial order within a subject, given unquoted. `NULL`
  uses the order the rows appear in.

## Value

A `frmtmb_family` object, for `frm(family = )`.

## Details

- `EV`, the expected value, updated by a delta rule with SEPARATE rates
  for a gain and a loss;

- `EF`, the expected frequency of a gain, updated by a delta rule on the
  SIGN of the outcome, and updated counterfactually on the decks that
  were not played;

- `PS`, perseverance, which is 1 for the deck just played and decays
  toward zero on the others.

The three combine linearly into the choice utility,

    util[j] <- EV[j] + betaF * EF[j] + betaP * PS[j]

with no inverse temperature: the softmax sensitivity is fixed at 1, and
the scale of the utilities is carried by the payoffs and by `betaF` and
`betaP` instead. That is hBayesDM's parameterization and it is kept, so
the payoff SCALE matters here in a way it does not for a family with a
free `tau`; see the section below.

This is the model hBayesDM calls `igt_orl`.

## How well it is identified, measured

The first release left this family out on the grounds that it is weakly
identified. Measured, that claim is half right and the half that is
wrong is the half that would have mattered.

Every parameter RECOVERS. At 30 subjects by 100 trials with a random
intercept on `Arew`, 60 replicates, every bias is at or inside its Monte
Carlo error and coverage runs 0.92 to 0.98; 60 of 60 fits gave a usable
interval. So a point estimate from this family is not the hazard a
variance component from a twenty-trial session is.

What is weak is the SEPARATION of the three choice-rule parameters, and
it shows up in the correlations of the estimates rather than in their
biases. Across replicates `betaF` and `betaP` correlate at -0.77, `k`
with `betaP` at -0.53 and `k` with `betaF` at +0.45. The two learning
rates are comparatively clean (`Arew` with `Apun`, 0.38) and neither
correlates with the choice-rule three above 0.21.

The practical consequence is specific. A study comparing groups on
`Arew` or `Apun` is on solid ground. A study comparing them on `betaF`
alone is at risk of attributing to outcome frequency what belongs to
perseverance, because the two are estimated against each other; report
them together, or test them jointly.
[`?frmtmb.learn`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frmtmb.learn-package.md)
carries the whole table.

## The payoff scale

hBayesDM divides the Iowa gambling task's payoffs by 100 before fitting,
and the reason is not cosmetic: with the softmax sensitivity fixed at 1,
doubling the payoffs doubles every utility and makes the model twice as
deterministic. `frm_task_design("igt")` already returns payoffs on that
scale. Data brought in raw dollars must be divided the same way, or
`betaF` and `betaP` are estimated against a `EV` a hundred times larger
than they were built for.

## What differs from hBayesDM

One transform. hBayesDM puts `K` on `(0, 5)` through a probit and decays
perseverance by `3^K`; this family estimates `k` on `(0, Inf)` with a
log link and decays by the same `3^k`. So the two `k` values are the
same number on the same scale, and only the bound differs: a fit here
can report a decay faster than hBayesDM's prior allows, and a fit that
wants to must be reported as such rather than silently pinned at 5.
Everything else is a rename: `Arew` and `Apun` keep their names, and
`betaF` and `betaP` keep theirs.

## The data a trial carries

`payoff(pay1, pay2, pay3, pay4)`, one column per deck, holding what each
deck WOULD have paid on that trial. The likelihood reads the played
deck's entry; the other three make the simulator coherent, and the
counterfactual frequency update reads only the SIGN the played deck
produced, not the others' values.

## References

Haines, N., Vassileva, J. and Ahn, W.-Y. (2018). The
outcome-representation learning model: a novel reinforcement learning
model of the Iowa gambling task. *Cognitive Science* 42, 2534-2561.

## See also

[`igt_pvl_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/igt_pvl_delta.md)
for the prospect-theory model of the same task,
[`frm_value_trace()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_value_trace.md)
for the three value stores per trial.

## Examples

``` r
d <- frm_task_design("igt", n_subject = 6, n_trial = 60, seed = 8)
d$choice <- frm_task_simulate(
  igt_orl(subject = id, trial = trial), d,
  pars = list(Arew = 0.3, Apun = 0.1, k = 0.5, betaF = 1,
              betaP = 1), seed = 8)[[1]]$choice
table(d$choice)
#> 
#>   1   2   3   4 
#>  47 106  82 125 
```
