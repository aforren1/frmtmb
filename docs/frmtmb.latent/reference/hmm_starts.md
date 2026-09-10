# Refits of an `hmm()` model from jittered starting values

A hidden Markov likelihood has more than one mode, and no convergence
diagnostic reports which one a fit reached. `hmm_starts()` refits the
model `n` times from starting values scattered around its estimates,
returns the best fit it found, and reports the spread of the OPTIMA the
refits reached, so that "the optimizer converged" and "the optimizer
converged to the answer" stop being the same sentence.

## Usage

``` r
hmm_starts(
  fit,
  n = 8L,
  jitter = 2,
  seed = NULL,
  grad_tol = 0.001,
  keep = FALSE
)
```

## Arguments

- fit:

  A `frmtmb_fit` with an
  [`hmm()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.md)
  family.

- n:

  Number of refits. Each one costs a whole fit.

- jitter:

  Perturbation size, in standard errors of the fit's own estimates.
  Larger values explore further and converge less often.

- seed:

  Optional integer. Given, the starts are reproducible and the session's
  random stream is restored afterwards.

- grad_tol:

  Largest `max|gradient| / |logLik|` a refit may have and still count as
  converged.

- keep:

  If `TRUE`, every refit that returned a fit at all is kept in `$fits`,
  in the order they ran, so it is shorter than `n` when a refit errored.
  These are whole fitted models with their tapes; `n` of them at a
  realistic size is gigabytes, which is why the default discards all but
  the best.

## Value

An object of class `frmtmb_hmm_starts`: `$best` (the highest-likelihood
fit, which is the original when no refit beat it, and whose own `call`
carries the winning start by value so that re-evaluating it reproduces
the fit), `$table` (one row per start: `logLik`, `grad_rel`, `status`,
`seconds`, and the `fn_evals` / `gr_evals` the optimizer took, which is
a cost measure that does not move with the machine's load), `$spread`
and `$spread_all` (`NA` when only one fit is in the set, since a
comparison that was never made is not a spread of zero), `$modes`,
`$mode_tol` (the one threshold above), and the counts.

## What the numbers mean

- the spread:

  The range of the log-likelihood over the refits that CONVERGED. A
  second spread is printed over every refit that finished, converged or
  not, because a summary that silently drops the runs that disagreed
  reports a spread that is tight for the wrong reason.

- converged:

  A refit counts as converged when it returns without error and its
  largest absolute gradient, divided by the absolute log-likelihood it
  reached, is at most `grad_tol`. The ratio rather than the gradient
  itself is what is tested, because on a long chain the objective's own
  magnitude sets the scale the optimizer stops at.

- the modes:

  Refits whose log-likelihoods agree to within `grad_tol^2` times the
  original fit's log-likelihood are counted as the same optimum. Two
  modes with a wide gap between them is the finding; one mode found `n`
  times is the reassurance.

That last threshold is returned as `$mode_tol`, and the same value
decides all three questions this function answers: whether a refit
displaces `$best`, whether the summary says the original found a local
optimum, and how the modes are merged. **It moves as the SQUARE of
`grad_tol`**, so loosening the convergence test by a factor of 100
loosens the merge threshold by 10 000. That is deliberate and it has a
range: on the probe `hmm_starts()` exists for, the 8.099-unit detection
survives every `grad_tol` up to about 0.08 and is gone at 0.086, where
the threshold has grown past the gap itself. The default of 1e-3 puts
the threshold at about 1e-03 log-likelihood units on a fit of -1096,
which is seven orders of magnitude below that gap.

## What it costs, and what `n` to use

A refit is a whole fit: the frame, the tape and the optimization again.
One call therefore costs `n` fits plus one `sdreport()`, which is paid
once for the jitter scale and is the same one
[`summary()`](https://rdrr.io/r/base/summary.html) would have paid.

Measured at the 25 000-row design of `dev/scale-findings.md` (50
sequences of 500, `K = 3` gaussian, `tr12 ~ (1 | id)`), alone on the
machine, three blocks of the same work from the same seed: one refit is
**175.7 s and 26.5 gradient evaluations**, against 183.5 s and 16
gradient evaluations for the plain fit it started from, and the
`sdreport()` is 37.2 s. So `n = 4` is 12.3 minutes there, `n = 8` is
24.0 and `n = 16` is 47.5. The instrument's control (the largest block
over the smallest) reported 1.023 and the gradient counts were identical
across blocks, which is the check a clock cannot make.
`dev/latent-2p3-starts-cost.R` is the script.

Twenty-four minutes is not something to do by accident on a model that
took three to fit, so on anything that size set `n` from the `seconds`
line of a first, small call rather than from the default.

`n` cannot have a default that is right at every size, and the one here
is a starting point rather than a recommendation. What it is chosen
from: on the probe in `dev/hmm-feasibility.md` a single jittered refit
at `jitter = 2` reaches the better optimum about 6 times in 10, so
`n = 8` misses it with probability about 4e-4 and `n = 4` with
probability about 2e-2. That rate is a property of that surface, not a
constant, so treat `n` as something to raise when the printed summary
shows more than one optimum and to price from the `seconds` line when
the model is large.

## Where the starts come from

Each start is the fit's own estimates plus a normal draw per outer
parameter, with standard deviation `jitter` times that parameter's
standard error. Standard errors make `jitter` a distance in the answer
rather than in an internal coordinate. Where the fit reports no usable
standard error, the median of the ones it does report is used, and where
it reports none at all every parameter is moved by `jitter` on the
internal scale; the printed summary names which of the three applied.
Random effects start from the original fit's conditional modes, since
the inner problem is re-solved either way.

## See also

[`hmm()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.md),
[`frmtmb::frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.html),
which varies the OPTIMIZER from one fixed start and answers a different
question.

## Examples

``` r
set.seed(4)
n_seq <- 12; len <- 20
G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
mu <- c(0, 3); sg <- c(0.6, 0.6)
dd <- do.call(rbind, lapply(seq_len(n_seq), function(id) {
  s <- integer(len); s[1] <- 1L
  for (t in 2:len) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
  data.frame(id = id, t = seq_len(len), y = rnorm(len, mu[s], sg[s]))
}))
fit <- frm(bf(y ~ 1),
           family = hmm(K = 2, gaussian(), time = t, group = id),
           data = dd)
ms <- hmm_starts(fit, n = 3, seed = 1)
ms
#> <hmm_starts> 3 refits, jitter 2 (standard errors)
#>   converged      3 of 3
#> 
#>   best logLik    -330.8885276
#>   original       -330.8885276   (the original fit was the best found)
#> 
#>   logLik spread, the original and the 3 converged refits: 7.310120509e-09
#>   logLik spread, every one of the 4 that finished    : 7.310120509e-09
#> 
#>   distinct optima among the original and the converged: 1
#>        logLik found
#>  -330.8885276     4
#> 
#>   seconds: 0.05641 total, 0.0179 median per refit
#>   gradient evaluations: 86 total, 28 median per refit
logLik(ms$best)
#> 'log Lik.' -330.8885 (df=6)
```
