# Most likely state path of an hmm fit (Viterbi)

Global decoding: the single state sequence with the highest posterior
probability, per sequence, by the Viterbi recursion. This is not the
same as taking the most likely state at each row separately
(`max.col(hmm_probs(fit))`, local decoding), and the two can disagree

- only the Viterbi path is guaranteed to be a possible path under the
  transition matrix.

## Usage

``` r
hmm_viterbi(fit)
```

## Arguments

- fit:

  A `frmtmb_fit` with an
  [`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md) family.

## Value

An integer vector, one state index per row of the data.

## Details

State LABELS are arbitrary: the likelihood is invariant to permuting the
states, so which run calls a state "1" is not meaningful across fits.
Compare paths after matching states by their fitted means.

## See also

[`hmm_probs()`](https://aforren1.github.io/frmtmb/reference/hmm_probs.md),
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
v <- hmm_viterbi(fit)
table(v, truth = s)
#>    truth
#> v    1  2
#>   1 58  0
#>   2  0 62
```
