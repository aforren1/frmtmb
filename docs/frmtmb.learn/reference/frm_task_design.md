# Build a trial-level design for a learning task

One row per subject and trial, with the payoff schedule of EVERY option
fixed in advance and a placeholder response. Fixing the whole schedule
ahead of the choices is what makes a draw from these models coherent: a
simulated subject that takes the arm the real one did not still needs to
be paid.

## Usage

``` r
frm_task_design(
  task = c("bandit2arm", "reversal", "bandit4arm_restless", "igt", "twostep"),
  n_subject = 30L,
  n_trial = 100L,
  p = c(0.7, 0.3),
  seed = NULL,
  ...
)
```

## Arguments

- task:

  One of `"bandit2arm"` (a stationary two-armed bandit), `"reversal"`
  (the same with the contingency reversed at the halfway point, and an
  `after_reversal` factor to model it with), `"bandit4arm_restless"`
  (four arms whose payoffs follow the decaying Gaussian random walk of
  Daw et al. 2006), `"igt"` (the Iowa gambling task's four decks) or
  `"twostep"` (the two-stage Markov task of Daw et al. 2011).

- n_subject, n_trial:

  Subjects, and trials per subject.

- p:

  Two payoff probabilities, for the two-armed tasks.

- seed:

  Passed to [`set.seed()`](https://rdrr.io/r/base/Random.html) when not
  `NULL`.

- ...:

  Task-specific settings: `decay` and `center` and `sd_walk` and
  `sd_obs` for the restless bandit, `p_reward` for the two-step task's
  stage-two options.

## Value

A data frame, ready for
[`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md).

## Details

The column names are canonical, because
[`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md)
reads them: `id`, `trial`, `choice`, and `pay1` to `pay2` or `pay4`. A
fit to your own data needs none of this; the family reads whatever
columns your formula names.

## See also

[`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md),
[`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md)

## Examples

``` r
d <- frm_task_design("reversal", n_subject = 4, n_trial = 20, seed = 1)
head(d)
#>   id trial choice after_reversal pay1 pay2
#> 1  1     1      1         before    1    0
#> 2  1     2      1         before    1    1
#> 3  1     3      1         before    1    0
#> 4  1     4      1         before    0    0
#> 5  1     5      1         before    1    1
#> 6  1     6      1         before    0    0
table(d$after_reversal)
#> 
#> before  after 
#>     40     40 
```
