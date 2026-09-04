# Changelog

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
  `.onLoad()`, so
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.html)
  gains their rows when this package is loaded and carries no dangling
  reference when it is not.

- No behavior changed in the move. Every test came over with its
  assertions intact.
