# Draw datasets from a learning family's generative process

Walks each subject forward through the design, drawing a choice from the
current value store and learning from the payoff that choice earns. It
is the same recursion the likelihood tapes, at `mode = "simulate"`, so a
draw and the density that scores it cannot drift apart.

## Usage

``` r
frm_task_simulate(
  family,
  data,
  pars,
  nsim = 1L,
  response = "choice",
  seed = NULL
)
```

## Arguments

- family:

  A family from this package, such as
  [`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md).

- data:

  A design from
  [`frm_task_design()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_design.md),
  or any data frame with the same canonical columns.

- pars:

  Named list of parameters on their NATURAL scales (a learning rate in
  `(0, 1)`, not its logit). Each element is one value, recycled
  everywhere; or one value per subject, in the order
  `levels(factor(data$id))` gives; or one value per ROW. The per-row
  form is what a covariate that varies WITHIN a subject needs, which is
  what a reversal is: see the example.

- nsim:

  Datasets to draw.

- response:

  Name of the column the drawn choice goes into.

- seed:

  Passed to [`set.seed()`](https://rdrr.io/r/base/Random.html) when not
  `NULL`.

## Value

A list of `nsim` data frames, each `data` with the drawn columns
replaced.

## Details

This route takes parameters directly, per subject or per row, rather
than through a formula and a random-effect covariance. That is what a
recovery study wants, and it makes the simulator independent of the
fitting machinery, so agreement between a fit and its truth is evidence
about both.
[`frmtmb::frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.html)
is the other route and goes through the fitted grammar; the two agree by
construction for every family except
[`ts_par7()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/ts_par7.md),
which
[`frmtmb::frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.html)
refuses because one trial's draw is three numbers.

## See also

[`frm_task_design()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_design.md),
[`frmtmb::frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.html)

## Examples

``` r
d <- frm_task_design("bandit2arm", n_subject = 5, n_trial = 30,
                     seed = 1)
sims <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial),
                          d, pars = list(alpha = 0.4, tau = 3),
                          nsim = 2, seed = 1)
table(sims[[1]]$choice)
#> 
#>   1   2 
#> 105  45 

# a learning rate that changes within a subject, one value per row
r <- frm_task_design("reversal", n_subject = 5, n_trial = 30, seed = 2)
a <- ifelse(r$after_reversal == "after", 0.6, 0.2)
rs <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial),
                        r, pars = list(alpha = a, tau = 3), seed = 2)
table(rs[[1]]$choice)
#> 
#>  1  2 
#> 74 76 
```
