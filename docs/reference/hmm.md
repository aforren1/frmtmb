# Hidden Markov models

`hmm(K, family)` fits a `K`-state hidden Markov model: the response at
each time point is drawn from one of `K` state-dependent copies of
`family`, and the unobserved state follows a first-order Markov chain
along `time` within `group`. The state sequence is summed out exactly by
the forward algorithm, which is evaluated on the same AD tape as
everything else, so random effects, smooths and distributional
predictors compose with it unchanged.

## Usage

``` r
hmm(
  K,
  family = stats::gaussian(),
  time = NULL,
  group = NULL,
  init = c("stationary", "estimated", "uniform"),
  trans = ~1
)
```

## Arguments

- K:

  Number of hidden states (at least 2, at most 9 - beyond that the
  `tr{i}{j}` dpar names stop being unambiguous).

- family:

  The state-dependent (emission) family. Any univariate family with a
  `mu` parameter works, plus
  [`multinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  for categorical emissions. Ordinal families, mixtures and nested
  `hmm()` components are refused.

- time:

  Variable giving the order of observations within a sequence, as a bare
  name or a one-sided formula. Omit it to use each row's position in the
  data.

- group:

  Variable identifying the sequences, as a bare name or a one-sided
  formula. Omit it to treat the whole data set as one sequence.

- init:

  Initial-state distribution: `"stationary"` (the default),
  `"estimated"`, or `"uniform"`.

- trans:

  Default one-sided formula for every transition cell.

## Value

A `frmtmb_family`.

## Parameters

Each of the wrapped family's distributional parameters is copied once
per state and suffixed with the state index, exactly as
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
does: `hmm(2, gaussian())` has `mu1`, `mu2`, `sigma1`, `sigma2`. The
main model formula applies to every state's location parameter; override
one state with `bf(y ~ x, mu2 ~ x + (1 | id))`. Every dpar takes the
full formula grammar, random effects included.

Transition probabilities are a row-wise multinomial logit with **state 1
as the reference cell in every row**, named `tr{i}{j}` for the move from
state `i` to state `j` (`j >= 2`): a two-state chain has `tr12` and
`tr22`. This is `depmixS4`'s parameterization, and the coefficients
agree with it to four decimals (`dev/hmm-feasibility.md`, probe B1).
`hmmTMB` and `moveHMM` instead reference the diagonal; the two span the
same model and only the reported coefficients differ. `trans` gives
every transition cell the same default formula; a single cell is
overridden the ordinary way, `bf(y ~ 1, tr12 ~ x)`. The covariate value
at time `t` drives the transition from `t` to `t + 1`, so a sequence's
last row never contributes a transition - again `depmixS4`'s convention.

## Initial distribution

- `"stationary"`:

  The stationary distribution of the transition matrix, solved on the
  tape. Costs no parameters and is the right default for long sequences
  and for many short ones. It needs a CONSTANT transition matrix, so it
  is refused when any transition dpar carries a predictor.

- `"estimated"`:

  `K - 1` free logits (state 1 the reference), reported as
  `hmm_ldel_1 ...` and counted in `df`.

- `"uniform"`:

  Fixed at `1 / K`.

## Decoding and the fitted values

`E[y_t]` under an HMM is `sum_k P(S_t = k | y) mu_k(x_t)`, which needs
the smoothed state probability and therefore a backward pass over the
whole sequence.
[`hmm_probs()`](https://aforren1.github.io/frmtmb/reference/hmm_probs.md)
returns those probabilities and
[`hmm_viterbi()`](https://aforren1.github.io/frmtmb/reference/hmm_viterbi.md)
the maximum-a-posteriori state path;
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
`predict(type = "response")` and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) all route
through
[`hmm_probs()`](https://aforren1.github.io/frmtmb/reference/hmm_probs.md),
so they report the occupancy-weighted mean rather than any single
state's.

## Label switching and local optima

The likelihood is invariant to permuting the states, so the state means
start at spread response quantiles (see
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md));
a start with every state mean equal sits on the symmetry axis and the
optimizer never leaves it. Relabeling between runs is expected and is
not fought. Multimodality is real and is not signalled by any
convergence diagnostic: on a random-effect model the default cold start
has been measured converging 8.1 log-likelihood units below the global
optimum with `convergence == 0` and a positive-definite Hessian. Compare
several starts
([`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md))
before reporting.

## Missing responses

An `NA` response is a time point the chain passes THROUGH without
emitting: the row is kept, its emission factor is masked out, and the
transition into and out of it is unchanged. This departs from the
package-wide `na.action`, deliberately - dropping the row would shorten
the chain and make one transition stand in for several, which biases the
transition matrix (measurably: with 3 of 20 points missing the fitted
off-diagonal moved from 0.115 to 0.128 against a true 0.10).
[`nobs()`](https://rdrr.io/r/stats/nobs.html) therefore counts every
row. [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) are `NA` at a
masked row, while
[`hmm_probs()`](https://aforren1.github.io/frmtmb/reference/hmm_probs.md)
is not: the neighbouring observations still say where the chain was.

## Cost

Evaluating the likelihood is linear in the number of ROWS and free in
the number of sequences: a thousand sequences of five cost what one of
five thousand does. Both the tape and the decoding passes are organized
by time step, so what they actually cost is
`length of the LONGEST sequence` times `K^2`. At 20 000 rows a fit takes
about 3 s cut into a thousand short sequences and about the same as one
long chain, but
[`hmm_probs()`](https://aforren1.github.io/frmtmb/reference/hmm_probs.md)
takes 0.06 s in the first shape and 1.3 s in the second. Below 5 000
rows nothing is noticeable. On a long chain the objective's own
magnitude also makes the optimizer's relative convergence test bite
before its gradient test does; judge `max|grad|` per observation rather
than in absolute terms.

## Boundaries

An unpenalized multinomial logit will send an emission or transition
probability to 0 whenever a category is rare inside a state, and the
optimizer then reports singular convergence at a perfectly good optimum.
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
on the affected logit is the remedy.

## Random effects and the Laplace approximation

A random effect in a state's linear predictor is integrated by the
Laplace approximation OUTSIDE the exact state sum. That integrand is a
mixture over state sequences and is not Gaussian, so the approximation
is genuinely approximate here even for a gaussian response. Measured
against adaptive Gauss-Hermite quadrature on a 40-sequence,
25-observation model with one scalar random intercept, the bias was
**-0.126 in the log-likelihood** (8.9e-5 relative) and 4.4e-4 absolute
in the parameters (`dev/hmm-feasibility.md`, probe D1).
`quadrature = TRUE` is refused: its rule integrates a random effect
against a PRODUCT of per-row densities, and the forward algorithm is not
one.

## What is refused

`REML` (the restricted likelihood would integrate out only one state's
location coefficients, which matches no standard definition),
`quadrature`, `frmtmb_control(profile = TRUE)`,
[`weights()`](https://rdrr.io/r/stats/weights.html), `cens()`,
[`trunc()`](https://rdrr.io/r/base/Round.html), `se()` and `mi()` on the
response, multivariate models and `rescor`, `residuals(type = "osa")`,
`predict(se.fit = TRUE)` on the response scale, and
[`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md).
A grouping in which every sequence has length 1 is refused too: the
chain is then unidentified and the model is a
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md).

## See also

[`hmm_probs()`](https://aforren1.github.io/frmtmb/reference/hmm_probs.md),
[`hmm_viterbi()`](https://aforren1.github.io/frmtmb/reference/hmm_viterbi.md),
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)

## Examples

``` r
set.seed(11)
n_seq <- 20; len <- 25
G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
mu <- c(0, 3); sg <- c(0.6, 0.6)
dd <- do.call(rbind, lapply(seq_len(n_seq), function(id) {
  s <- integer(len); s[1] <- 1L
  for (t in 2:len) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
  data.frame(id = id, t = seq_len(len),
             y = rnorm(len, mu[s], sg[s]), state = s)
}))
fit <- frm(bf(y ~ 1),
           family = hmm(K = 2, gaussian(), time = t, group = id),
           data = dd)
fixef(fit)
#> $mu1
#> (Intercept) 
#> -0.01572669 
#> 
#> $sigma1
#> (Intercept) 
#>  -0.5834065 
#> 
#> $mu2
#> (Intercept) 
#>    3.053871 
#> 
#> $sigma2
#> (Intercept) 
#>  -0.5550471 
#> 
#> $tr12
#> (Intercept) 
#>   -2.274024 
#> 
#> $tr22
#> (Intercept) 
#>    1.117463 
#> 

# smoothed state probabilities and the MAP path
head(hmm_probs(fit))
#>         state1       state2
#> [1,] 1.0000000 1.420375e-11
#> [2,] 1.0000000 4.670620e-09
#> [3,] 1.0000000 3.274867e-11
#> [4,] 1.0000000 2.084125e-08
#> [5,] 1.0000000 5.615109e-09
#> [6,] 0.9999152 8.478153e-05
mean(hmm_viterbi(fit) == dd$state)
#> [1] 1

# fitted() is the occupancy-weighted mean, not state 1's
cor(fitted(fit), dd$y)
#> [1] 0.9299445

# \donttest{
# one state's mean takes its own predictor, random effects included
frm(bf(y ~ 1, mu2 ~ 1 + (1 | id)),
    family = hmm(K = 2, gaussian(), time = t, group = id),
    data = dd)
#> frmtmb fit: y ~ 1 
#> Family: hmm(2, gaussian)   Method: ML 
#> logLik: -617.932  AIC: 1249.86  nobs: 500 
#> 
#> Fixed effects:
#>  mu1:
#> (Intercept) 
#>    -0.01573 
#>  sigma1:
#> (Intercept) 
#>     -0.5834 
#>  mu2:
#> (Intercept) 
#>       3.054 
#>  sigma2:
#> (Intercept) 
#>      -0.555 
#>  tr12:
#> (Intercept) 
#>      -2.274 
#>  tr22:
#> (Intercept) 
#>       1.117 
#> 
#> Random effects:
#>   mu2: 1 | id 
#>         Name  Std.Dev.
#>  (Intercept) 7.362e-05

# covariate-dependent transitions: trans = sets every cell's default
dd$x <- rnorm(nrow(dd))
frm(bf(y ~ 1),
    family = hmm(K = 2, gaussian(), time = t, group = id,
                 init = "estimated", trans = ~x),
    data = dd)
#> frmtmb fit: y ~ 1 
#> Family: hmm(2, gaussian)   Method: ML 
#> logLik: -610.619  AIC: 1239.24  nobs: 500 
#> 
#> Fixed effects:
#>  mu1:
#> (Intercept) 
#>    -0.01582 
#>  sigma1:
#> (Intercept) 
#>     -0.5836 
#>  mu2:
#> (Intercept) 
#>       3.054 
#>  sigma2:
#> (Intercept) 
#>     -0.5546 
#>  tr12:
#> (Intercept)           x 
#>    -2.10651     0.08341 
#>  tr22:
#> (Intercept)           x 
#>     1.29178    -0.05649 
# }
```
