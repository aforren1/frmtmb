INTERIM, six of the seven lanes (illse not merged). Seven lanes, each with an adversarial review and fix rounds:
`dblbar` (every `||` term read as brms reads it), `famsimplex` (brms's
simplex families, `dirichlet()` and `logistic_normal()`), `famuni`
(`frechet()`, `gen_extreme_value()` and `discrete_weibull()`, and
frmtmb.eam's `wiener()` links), `postfix` (the post-processing gaps of
the brms port: `frm_multiple()` methods, `pp_check(resp = )`, and the
default-start error of a nonlinear fit), `optstart` (least-squares
starting values for location-scale fits, and a restart that neither
gains nor moves), `illse` (accurate standard errors on ill-conditioned
designs) and `remeasure` (the vignette port and bin 2 of brms's suite,
measured; no package code). Each lane's `dev/<lane>-findings.md` has
the validation, the numbers and the scripts, and
`dev/reviews/2026-10-13-<lane>.md` has the review;
`dev/round-20261014.md` is the round's record.
