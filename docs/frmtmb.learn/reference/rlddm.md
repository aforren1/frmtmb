# Reinforcement learning with a drift-diffusion choice rule

A delta learning rule feeding the drift rate of a two-boundary Wiener
diffusion, so that the choice and the response time are one likelihood
rather than two. The learning half is
[`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md)'s
exactly:

## Usage

``` r
rlddm(subject, trial = NULL, max_ndt = NULL)
```

## Arguments

- subject:

  The column separating one learner's trial sequence from the next,
  given unquoted.

- trial:

  The column giving trial order within a subject, given unquoted. `NULL`
  uses the order the rows appear in.

- max_ndt:

  Upper bound for the non-decision time. `NULL`, the default, uses the
  fastest response in the data. A value above it is refused, because it
  admits parameters at which the fastest trial has no likelihood.

## Value

A `frmtmb_family` object, for `frm(family = )`.

## Details

    Q[chosen] <- Q[chosen] + alpha * (reward - Q[chosen])

and the choice half replaces the softmax with a first-passage density.
Evidence accumulates at a rate proportional to how far apart the two
value estimates have grown,

    v <- drift * (Q[upper] - Q[lower])

between boundaries `bs` apart, starting a fraction `bias` of the way up,
and the response time is `ndt` plus the time the accumulator takes to
touch a boundary. A trial's contribution is then the joint density of
WHICH boundary was touched and WHEN, which is what makes this a model of
response times rather than a choice model with a timing footnote: a
learning rate estimated from choices alone cannot tell a subject who
deliberated from one who guessed, and this one can.

The model is Pedersen, Frank and Biele (2017). hBayesDM does not carry
it; see the section below for what it does carry.

## What it corresponds to in hBayesDM

Nothing, and that is worth stating plainly because every other family in
this package has a counterpart. hBayesDM ships `choiceRT_ddm`, which is
this family's choice rule with NO learning: one drift rate per subject
and no value store. It also ships `bandit2arm_delta`, which is this
family's learning rule with a softmax instead of a diffusion. `rlddm()`
is the two joined, and the join is the model rather than a renaming of
either.

Where the names do map, they map as hBayesDM's `choiceRT_ddm` spells the
diffusion: its `alpha` is the boundary separation, which is `bs` here,
its `beta` is the start point, which is `bias`, its `delta` is the drift
rate, which is `drift` times the value difference, and its `tau` is the
non-decision time, which is `ndt`. `alpha` here is the LEARNING RATE, as
it is in every other family in this package. That collision is why the
diffusion parameters carry `frmtmb.eam`'s names rather than hBayesDM's:
one convention across this package's families is worth more than an
inherited one that would make `alpha` mean two different things two
families apart.

## The data a trial carries

The response is the response TIME. The boundary reached travels in
`dec()`, coded 0 for the lower boundary and 1 for the upper, which is
the term and the coding
[`frmtmb.eam::wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.html)
uses. Arm 1 is the lower boundary and arm 2 the upper, so
`reward(pay1, pay2)` is in the same arm order as every other family
here.

    frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1 + (1 | id),
           drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1),
        family = rlddm(subject = id, trial = trial), data = d)

## Holding the start point at a half

Most of the reinforcement-learning literature fixes `bias` at 0.5,
because with two ARMS rather than a correct and an error response there
is no reason for the accumulator to start nearer one boundary. It is an
ordinary distributional parameter here, so it is estimated by default
and held with
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.html)'s constant
form:

    frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1 + (1 | id),
           drift ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5), ...)

Estimating it is the more general model and costs one parameter; fixing
it is what makes `drift` and `bs` easier to separate on a short session.

## The non-decision time is bounded, not logged

The density is zero at and below `ndt`, so the likelihood has a hard
edge at `ndt = min(rt)` and a log link would let the optimizer walk over
it. `ndt` therefore gets a logit scaled onto `(0, max_ndt)`, with
`max_ndt` defaulting to the fastest response in the data, which makes
the constraint structural. This is the same construction
[`frmtmb.eam::wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.html)
uses and it is written again here rather than borrowed, because a link
is not what that package exports.

## Recovery, measured

`dev/learn-recovery-round2.R` runs it and
[`?frmtmb.learn`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frmtmb.learn-package.md)
carries the table. All five parameters recover at 30 subjects by 100
trials with a random intercept on the learning rate: every bias is
inside its Monte Carlo error, and 60 of 60 replicates produced a usable
interval. Coverage is at the nominal rate for `drift`, `bs`, `ndt` and
`bias`, and 0.87 for `alpha`, so a learning rate's Wald interval is
slightly optimistic here and should be read as such.

THE PAIR THAT TRADES OFF IS NOT THE ONE TO EXPECT, and this section said
otherwise before the study was run. A first-passage density is driven
largely by the RATIO of the drift to the boundary, so `drift` and `bs`
look like the pair at risk; measured, their estimates correlate at 0.17
across replicates, which is nothing. What does co-vary is `bs` with
`ndt`, at -0.57, and `drift` with `bias`, at +0.56. Both make sense
after the fact: the boundary separation and the non-decision time are
the two ways to make responses slower, and the drift and the start point
are the two ways to favor a boundary. Read those two pairs together.

## Drawing data needs RWiener, and fitting does not

[`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md)
on this family draws a boundary and a time JOINTLY from the diffusion,
through
[`frmtmb.eam::ddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ddm_simulate.html),
and that route is exact only when `RWiener` is installed. Without it the
draw is refused rather than approximated: `frmtmb.eam`'s own fallback
draws the boundary first and then rejects paths until one reaches it,
which fails outright on the rows where one boundary is strongly favored,
and a simulator that works on most rows and errors on the rest is worse
than one that says what it needs. Nothing about FITTING needs `RWiener`;
the likelihood is
[`frmtmb.eam::wiener_lpdf()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_lpdf.html)
and is exact either way.

## Where the density comes from

[`frmtmb.eam::wiener_lpdf()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_lpdf.html),
through the export rather than through a colon. The Navarro and Fuss
(2009) series pair, its smooth blend and its tape safety all belong to
that package and are not repeated here. This is the only dependency one
extension of frmtmb has on another, and it is one function.

## References

Pedersen, M. L., Frank, M. J. and Biele, G. (2017). The drift diffusion
model as the choice rule in reinforcement learning. *Psychonomic
Bulletin and Review* 24, 1234-1251.

Navarro, D. J. and Fuss, I. G. (2009). Fast and accurate calculations
for first-passage times in Wiener diffusion models. *Journal of
Mathematical Psychology* 53, 222-230.

## See also

[`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md)
for the same learning rule under a softmax,
[`frmtmb.eam::wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.html)
for the diffusion without learning,
[`frm_value_trace()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_value_trace.md)
for the per-trial drift rates.

## Examples

``` r
d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 40,
                     seed = 1)
s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
                       pars = list(alpha = 0.4, drift = 3, bs = 1.6,
                                   ndt = 0.2, bias = 0.5), seed = 1)
fit <- frmtmb::frm(
  frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) ~ 1,
             drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1),
  family = rlddm(subject = id, trial = trial), data = s[[1]])
frmtmb::fixef(fit)
#> $alpha
#> (Intercept) 
#>  -0.1298849 
#> 
#> $drift
#> (Intercept) 
#>    3.003392 
#> 
#> $bs
#> (Intercept) 
#>   0.4032789 
#> 
#> $ndt
#> (Intercept) 
#>    1.303546 
#> 
#> $bias
#> (Intercept) 
#>  0.05694329 
#> 
head(frm_value_trace(fit))
#>   subject trial        q1 q2   drift_t         pe       dens
#> 1       1     1 0.0000000  0  0.000000  0.0000000 0.76659084
#> 2       1     2 0.0000000  0  0.000000  1.0000000 0.31815357
#> 3       1     3 0.4675743  0 -1.404309  0.5324257 0.04651865
#> 4       1     4 0.7165229  0 -2.151999 -0.7165229 2.11932154
#> 5       1     5 0.3814952  0 -1.145779  0.6185048 1.63379607
#> 6       1     6 0.6706922  0 -2.014351 -0.6706922 2.42811813
```
