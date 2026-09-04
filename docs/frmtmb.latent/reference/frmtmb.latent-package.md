# frmtmb.latent: Latent-State Families for 'frmtmb' Models

Adds two families whose latent variable is a discrete state to the model
grammar of 'frmtmb'. 'hmm()' fits hidden Markov models with
covariate-dependent transitions, summing the state sequence out exactly
with the forward algorithm and decoding it with the forward-backward and
Viterbi passes. 'lca()' fits latent class analysis in the manner of
'poLCA', with conditionally independent polytomous items and a
multinomial logit on the class membership. Both are written entirely
against the structured-family protocol that 'frmtmb' exports, so neither
needs a branch inside the core package.

## See also

Useful links:

- <https://aforren1.github.io/frmtmb/>

- <https://github.com/aforren1/frmtmb>

- Report bugs at <https://github.com/aforren1/frmtmb/issues>

## Author

**Maintainer**: Alex Forrence <alex.forrence@gmail.com>
([ORCID](https://orcid.org/0000-0002-9728-6337))

Authors:

- Alex Forrence <alex.forrence@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-9728-6337))
