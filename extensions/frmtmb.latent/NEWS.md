# frmtmb.latent 0.2.0

`hmm()` declares the per-sequence log-likelihood, so `frm(importance
= )` and `loo()` on draws reach it. Requires frmtmb 0.53.0 for the
factorization slot.

* `hmm()` declares `frmtmb_structure(loglik_group = )`: the
  log-likelihood of every sequence, which the forward recursion
  already produces on its way to the total. `frm(importance = )` can
  correct an hmm() fit whose grouping factor is its own sequence
  variable, and `loo()` on sampled draws has one column per sequence
  instead of a refusal. There is no `loglik_row`, and its absence is
  the answer rather than an omission: a row's emission density is not
  its contribution to the likelihood, because the state that emitted
  it was reached through every earlier row. `residuals(type =
  "deviance")` stays refused for exactly that reason, in the words it
  already used.

* `lca()` declares neither slot, and says why in the source: its
  likelihood IS rowwise, because one row is one subject's whole item
  response pattern, so `loo()` and `frm(importance = )` reach it
  through the ordinary path and a slot would be a second definition of
  the same numbers.

# frmtmb.latent 0.1.0

First release. `hmm()` and `lca()` were part of frmtmb through
v0.47.0 and move here unchanged, as step 10 (the last) of
dev/structured-family-protocol.md.

* `hmm(K, family)` fits hidden Markov models: covariate-dependent
  transitions, a choice of initial distribution, the forward algorithm
  for the exact likelihood, and `hmm_probs()` and `hmm_viterbi()` for
  the smoothed occupancies and the decoded path.

* `lca(K)` fits latent class analysis in the manner of poLCA:
  conditionally independent polytomous items, class-conditional item
  profiles as family extra parameters, and a multinomial logit on the
  class membership so that latent class REGRESSION is the ordinary
  fixed-effect machinery. `lca_probs()` and `lca_profiles()` read the
  fit.

* Both reach frmtmb only through its exported extension API, and the
  suite asserts it: `test-structure-latent.R` scans every function in
  this namespace and fails if one reaches a frmtmb internal.

* The compatibility rules that name either family register from
  `.onLoad()`, so `frm_compat()` gains their rows when this package is
  loaded and carries no dangling reference when it is not.

* No behavior changed in the move. Every test came over with its
  assertions intact.
