Adds two families whose latent variable is a discrete state to the model
grammar of 'frmtmb'. 'hmm()' fits hidden Markov models with
covariate-dependent transitions, summing the state sequence out exactly
with the forward algorithm and decoding it with the forward-backward and
Viterbi passes. 'lca()' fits latent class analysis in the manner of
'poLCA', with conditionally independent polytomous items and a
multinomial logit on the class membership. Both are written entirely
against the structured-family protocol that 'frmtmb' exports, so neither
needs a branch inside the core package.
