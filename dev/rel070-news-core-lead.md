Three lanes, each with an adversarial review and punch rounds:
`moparity` (brms features that were refused or diverged: `mo()`
crossed with a factor and in group-level terms, Dirichlet priors on a
simplex, `update(prior =)` merging as brms merges, the `cox()`
baseline), `optrap` (fits that stopped short of the maximum because of
the response's units, a restart that raised, or a badly conditioned
design) and `seflags` (which standard errors a fit reports: an exact
design alias, curved ridges, an sd the likelihood does not read, raw
polynomials). Each lane's `dev/<lane>-findings.md` has the validation,
the numbers and the scripts, and `dev/reviews/2026-10-08-<lane>.md`
has the review; `dev/round-20261013.md` is the round's record.
In brms's own ported suite, bin 1 passes 389 of 494 assertions (390 at
0.69.0): `brmsfit-methods:326` asks a fixture that does not converge
for finite standard errors, and it held at 0.69.0 by the luck of where
the optimizer stopped.
