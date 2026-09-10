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
  [`frmtmb::multinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.html)
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
[`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)
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

`"estimated"` estimates ONE initial distribution shared across every
sequence. This is worth knowing before comparing with `hmmTMB`, whose
`initial_state = "estimated"` estimates a SEPARATE one per sequence. On
20 sequences of 100 rows that is 38 extra parameters on hmmTMB's side
and **19.678 log-likelihood units**, with nothing in either package's
output saying the two fits are of different models. (Against
`init = "uniform"` here the same hmmTMB fit is 21.231 units higher,
because uniform costs frmtmb its own two initial-distribution parameters
as well; the two figures differ by that 1.553 and it is the first that
belongs to the sentence above.)

Pin the initial distribution on both sides before comparing.
`init = "uniform"` here against `fixpar = list(delta0 = <all NA>)` there
makes the two likelihoods the same function, and they then agree to
1e-11 relative (`dev/latent-2p3-hmmtmb-probe4.R`).

## Decoding and the fitted values

`E[y_t]` under an HMM is `sum_k P(S_t = k | y) mu_k(x_t)`, which needs
the smoothed state probability and therefore a backward pass over the
whole sequence.
[`hmm_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_probs.md)
returns those probabilities and
[`hmm_viterbi()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_viterbi.md)
the maximum-a-posteriori state path;
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
`predict(type = "response")` and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) all route
through
[`hmm_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_probs.md),
so they report the occupancy-weighted mean rather than any single
state's.

## Label switching and local optima

The likelihood is invariant to permuting the states, so the state means
start at spread response quantiles (see
[`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html));
a start with every state mean equal sits on the symmetry axis and the
optimizer never leaves it. Relabeling between runs is expected and is
not fought.

Multimodality is real and NO convergence diagnostic reports it. On a
two-state model with a random intercept in each state mean, 25 sequences
of 30, the default cold start converges to -1096.09575602 where the
optimum is -1087.99646521: **8.099 log-likelihood units below it**, with
`convergence == 0`, a positive definite Hessian, no non-finite standard
error, a largest gradient of 7.00028e-04 (6.4e-07 of the
log-likelihood), and `diagnose()` printing "No convergence problems
detected". That last clause is true of this fit by a factor of 1.43:
`diagnose()`'s clean verdict needs `max|grad|` under 1e-3 ABSOLUTE, and
a construction whose gradient landed just the other side of that would
get a warning instead. What is 8.099 units wrong here is the answer, not
the gradient, and no threshold on a gradient can see that. The optimum
is confirmed independently by a hand-rolled `MakeADFun(random =)` and by
hmmTMB to the last digit (`dev/hmm-feasibility.md`, probe D4;
`dev/latent-2p3-repro81.R` reproduces it against this family).

[`hmm_starts()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_starts.md)
is the remedy: it refits from jittered starting values, returns the
best, and reports the spread of the OPTIMA the refits reached. On that
probe, `hmm_starts(fit, n = 8, jitter = 2)` recovered the better optimum
on 5 of 5 seeds, and so did `jitter = 1`, `4` and `8`; `jitter = 0.5`
recovered it on 1 of 5, which is what the default of 2 is chosen from.

[`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.html)
is a different check, and on this probe it is a reassuring one: it
re-runs every optimizer from the SAME start, so it tests the optimizer
rather than the surface. All four (nlminb, optim, bobyqa, NLopt L-BFGS)
reach -1096.096 with `convergence == 0` and a log-likelihood spread of
7.8e-07, which is agreement on the wrong answer
(`dev/latent-2p3-allfit.R`).

## Recovery at a realistic scale

Measured at the design named in the realistic-scale table of
`dev/extension-gaps-plan.md`: 50 sequences of 500 steps, 25 000 rows,
`K = 3` gaussian with state means -2, 0 and 3 and a common standard
deviation of 0.7, and a transition matrix whose free logits are `tr12`
-1.5, `tr13` -2.5, `tr22` 2.0, `tr23` -1.0, `tr32` -1.0 and `tr33` 2.0.
Scripts and tables: `dev/latent-2p3-hmm.R` and `dev/latent-findings.md`.

**Fixed transitions, against depmixS4.** Recovery over 40 replicates
(seeds 20260910 to 20260949): the largest bias on any of the twelve
parameters is 0.014 on a transition logit whose own standard error is
0.083, the spread across replicates matches the standard error each fit
reports, and Wald coverage is 95.6 percent (459 of 480, Monte Carlo
interval 93.8 to 97.5). On the response scale the state means come back
within 0.021, the state standard deviations within 0.017, and every cell
of the transition matrix within 0.014. The fitted labels were the true
labels in 40 of 40, and the Hessian was positive definite in 40 of 40.
Against depmixS4 on the same data, best of four random EM starts at
`tol = 1e-12`, over 3 replicates: the two log-likelihoods agree to
between 3.1e-14 and 4.9e-12 relative, and the means, standard deviations
and transition matrix to 1.5e-06, which is 1.6e-04 of one of this fit's
own standard errors.

**`tr12 ~ (1 | id)`, against hmmTMB.** Recovery over 15 replicates, with
a true `sd(tr12 | id)` of 0.6: it comes back at 0.615 on average with a
spread of 0.067, and Wald coverage over all thirteen parameters is 94.4
percent (184 of 195, Monte Carlo interval 91.1 to 97.6). Against hmmTMB
fitting the same model on the same data, over 6 replicates, the two
Laplace marginal log-likelihoods agree to between 6.1e-12 and 6.9e-11
relative and the estimates to 1.3e-04, which is 1.5e-02 of one standard
error. Fifteen replicates put about ten points of Monte Carlo error on
any single parameter's coverage, so read the overall figure rather than
the rows.

No replicate of either arm converged to a local optimum. That is a
statement about these two designs and not about the family: see "Label
switching and local optima" for a design where the cold start does
exactly that, with every diagnostic clean.

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
[`hmm_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_probs.md)
is not: the neighbouring observations still say where the chain was.

## Cost

Evaluating the likelihood is linear in the number of ROWS and free in
the number of sequences: a thousand sequences of five cost what one of
five thousand does. Both the tape and the decoding passes are organized
by time step, so what they actually cost is
`length of the LONGEST sequence` times `K^2`. At 20 000 rows a fit takes
about 3 s cut into a thousand short sequences and about the same as one
long chain, but
[`hmm_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_probs.md)
takes 0.06 s in the first shape and 1.3 s in the second. Below 5 000
rows nothing is noticeable. On a long chain the objective's own
magnitude also makes the optimizer's relative convergence test bite
before its gradient test does; judge `max|grad|` per observation rather
than in absolute terms.

## Boundaries

An unpenalized multinomial logit will send an emission or transition
probability to 0 whenever a category is rare inside a state, and the
optimizer then reports singular convergence at a perfectly good optimum.
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.html)
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
`conditional_effects()`. A grouping in which every sequence has length 1
is refused too: the chain is then unidentified and the model is a
[`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html).

## See also

[`hmm_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_probs.md),
[`hmm_viterbi()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_viterbi.md),
[`hmm_starts()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_starts.md),
[`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)

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
#>  Links: mu1 = identity; sigma1 = log; mu2 = identity; sigma2 = log;
#>         tr12 = identity; tr22 = identity
#> 
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
#>  Links: mu1 = identity; sigma1 = log; mu2 = identity; sigma2 = log;
#>         tr12 = identity; tr22 = identity
#> 
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
