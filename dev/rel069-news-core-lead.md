Four lanes, each with an adversarial review and punch rounds:
`ciharden` (tests and CI made independent of platform rounding, CI
pinned to `ubuntu-24.04` with a Stan job on Ubuntu 26.04 beside it,
less memory in the `gp()` kriging covariance), `surface` (the brms
post-processing calls a ported script meets on a maximum-likelihood
fit now answer or refuse by name), `setier` (which standard errors a
fit can report, a boundary fit said as lme4 says it, separation named)
and `optima` (fits that stopped short of the maximum likelihood:
`mo()`, the probit tails, `cs()` ordinal mixtures, nlminb's rejected
trial point). Each lane's `dev/<lane>-findings.md` has the validation,
the numbers and the scripts, and `dev/reviews/2026-10-07-<lane>.md`
has the review; `dev/round-20261007.md` is the round's record.
In brms's own ported suite, bin 1 now passes 390 of 494 assertions,
up from 389 at 0.68.1.
