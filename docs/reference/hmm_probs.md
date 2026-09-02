# Posterior state probabilities of an hmm fit

The smoothed occupancy `P(S_t = k | y)` for every row of the data, from
a forward-backward pass at the estimates. This is the
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
analog
[`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md),
one rung up: an HMM's per-row state probability conditions on the WHOLE
sequence, not on that row alone, which is why it needs a backward pass
and cannot come out of the family's density.

## Usage

``` r
hmm_probs(fit)
```

## Arguments

- fit:

  A `frmtmb_fit` with an
  [`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md) family.

## Value

An `n x K` matrix of probabilities whose rows sum to one, with columns
named `state1 ... stateK`.

## Details

Under a random-effect model the probabilities are conditional on the
random-effect modes, as every other post-fit quantity in the package is.
Rows whose response is `NA` still get a probability: the chain passes
through them and the neighboring observations inform the state.

## See also

[`hmm_viterbi()`](https://aforren1.github.io/frmtmb/reference/hmm_viterbi.md),
[`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md)

## Examples

``` r
set.seed(12)
dd <- data.frame(id = 1, t = 1:120)
s <- integer(120); s[1] <- 1L
G <- matrix(c(0.92, 0.08, 0.15, 0.85), 2, 2, byrow = TRUE)
for (t in 2:120) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
dd$y <- rnorm(120, c(0, 3)[s], 0.5)
fit <- frm(bf(y ~ 1),
           family = hmm(K = 2, gaussian(), time = t, group = id),
           data = dd)
p <- hmm_probs(fit)
head(p)
#>            state1       state2
#> [1,] 9.999989e-01 1.121901e-06
#> [2,] 1.000000e+00 1.786469e-08
#> [3,] 9.999676e-01 3.235686e-05
#> [4,] 4.209560e-08 1.000000e+00
#> [5,] 2.937353e-10 1.000000e+00
#> [6,] 8.655297e-10 1.000000e+00
rowSums(p)[1:3]
#> [1] 1 1 1
mean(max.col(p) == s)
#> [1] 1
```
