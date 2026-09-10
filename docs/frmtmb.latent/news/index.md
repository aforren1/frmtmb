# Changelog

## frmtmb.latent 0.3.0

- **[`hmm_starts()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_starts.md)
  refits from jittered starting values and reports the spread of the
  optima.** The default cold start has been measured converging **8.099
  log-likelihood units below the optimum** with `convergence` 0, a
  positive definite Hessian, and `diagnose()` printing “No convergence
  problems detected”; `frm_allfit()` reaches that same wrong optimum on
  all four optimizers with a spread of 7.8e-07, so agreement between
  optimizers is not evidence here.
  [`hmm_starts()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_starts.md)
  recovers from it on 5 of 5 seeds at jitter 1 and above. It costs about
  176 seconds per refit at 50 sequences by 500 steps, so
  [`?hmm`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.md)
  deliberately recommends no default `n`.

- **BREAKING:
  [`lca()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.md)’s
  starting values changed, and a fit can move.** The old rule cut
  subjects on the mean of their item codes, which is blind to any design
  whose classes differ in WHICH items they endorse rather than how many.
  On the realistic-scale design that score explains 0.000814 of its own
  variance between the true classes, and the old rule reached a local
  optimum 243 to 284 log-likelihood units below `poLCA(nrep = 10)` on 8
  of 200 replicates, with a positive definite Hessian on seven of the
  eight. On at least two of those the cause was the start and not the
  surface: all 10 of 10 poLCA single random starts found the global
  optimum there.

  The rule now clusters on the response pattern and shrinks each class
  profile nine tenths of the way toward the pooled one. It reaches
  poLCA’s optimum on 200 of 200 replicates of the block its one constant
  was tuned on, and on 200 of 200 of a second block of seeds it never
  saw.

  A fit therefore moves where the two starts differ, and it moves to the
  better optimum. Class 1 is still the low-score end wherever that score
  means anything; where it does not, the labeling is deterministic but
  need not match 0.2.2. **A 0.2.2 fit can be reproduced**:
  [`?lca`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.md)
  carries the old rule as ten lines of R, checked to reproduce the old
  starting values and the old log-likelihood at a relative difference of
  exactly 0.

- **Two coefficients under-cover, and there is no remedy yet.** Pooled
  over 400 replicates of the realistic-scale design, `class3:x2` covers
  91.75 percent and `class4:x2` 90.75 against a nominal 95, while the
  third binary gating slope and the other six coefficients contain 95. A
  profile interval was measured and does not help, 103 of 118 both ways,
  because the standardized error has no heavy tail and no curvature: the
  shortfall is in the magnitude of the standard error rather than the
  shape of the likelihood.
  [`?lca`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.md)
  says so rather than recommending a remedy nobody has measured.

## frmtmb.latent 0.2.2

- Requires frmtmb 0.55.1. A scale row in the gated measurement tier,
  which finds the hidden Markov post-fit passes are 0.25 and 0.09
  seconds at 50 sequences by 500 steps, so the R loops there are not a
  cost worth removing.

## frmtmb.latent 0.2.1

Requires frmtmb 0.55.0, for the hazard-container lint that now runs in
this package’s own check.

## frmtmb.latent 0.2.0

[`hmm()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.md)
declares the per-sequence log-likelihood, so `frm(importance = )` and
`loo()` on draws reach it. Requires frmtmb 0.53.0 for the factorization
slot.

- [`hmm()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.md)
  declares `frmtmb_structure(loglik_group = )`: the log-likelihood of
  every sequence, which the forward recursion already produces on its
  way to the total. `frm(importance = )` can correct an hmm() fit whose
  grouping factor is its own sequence variable, and `loo()` on sampled
  draws has one column per sequence instead of a refusal. There is no
  `loglik_row`, and its absence is the answer rather than an omission: a
  row’s emission density is not its contribution to the likelihood,
  because the state that emitted it was reached through every earlier
  row. `residuals(type = "deviance")` stays refused for exactly that
  reason, in the words it already used.

- [`lca()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.md)
  declares neither slot, and says why in the source: its likelihood IS
  rowwise, because one row is one subject’s whole item response pattern,
  so `loo()` and `frm(importance = )` reach it through the ordinary path
  and a slot would be a second definition of the same numbers.

## frmtmb.latent 0.1.0

First release.
[`hmm()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.md)
and
[`lca()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.md)
were part of frmtmb through v0.47.0 and move here unchanged, as step 10
(the last) of dev/structured-family-protocol.md.

- `hmm(K, family)` fits hidden Markov models: covariate-dependent
  transitions, a choice of initial distribution, the forward algorithm
  for the exact likelihood, and
  [`hmm_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_probs.md)
  and
  [`hmm_viterbi()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm_viterbi.md)
  for the smoothed occupancies and the decoded path.

- `lca(K)` fits latent class analysis in the manner of poLCA:
  conditionally independent polytomous items, class-conditional item
  profiles as family extra parameters, and a multinomial logit on the
  class membership so that latent class REGRESSION is the ordinary
  fixed-effect machinery.
  [`lca_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_probs.md)
  and
  [`lca_profiles()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_profiles.md)
  read the fit.

- Both reach frmtmb only through its exported extension API, and the
  suite asserts it: `test-structure-latent.R` scans every function in
  this namespace and fails if one reaches a frmtmb internal.

- The compatibility rules that name either family register from
  `.onLoad()`, so `frm_compat()` gains their rows when this package is
  loaded and carries no dangling reference when it is not.

- No behavior changed in the move. Every test came over with its
  assertions intact.
