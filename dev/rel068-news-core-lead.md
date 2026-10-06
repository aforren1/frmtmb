Five lanes, each with an adversarial review and punch rounds: `fixes`
(`emmeans()` on transformed predictors, `trunc()` and `se()` in
`conditional_effects()`, the ordinal leftovers of 0.67.0, brms's `s()`
null-space basis, the nonlinear flat-direction warning), `gpby`
(`gp(x, by = )` with brms's arguments and prior classes, and the exact
`gp()` kriging covariance at new positions), `ordmix` (ordinal
mixtures, `hurdle_cumulative()` with `thres(gr = )` and `cs()`),
`nanse` (a fit says which parameters have no standard error, and why,
and `frm_allfit()` starts where `lme4::allFit()` does) and `vigport`
(the brms vignette port measured again, records only). The
consolidation lifted the refusal of `cs()` on `cumulative()` and made
`conditional_effects(method = "predict")` report the median, both the
user's decisions of 2026-10-06. Each lane's `dev/<lane>-findings.md`
has the validation, the numbers and the scripts, and
`dev/reviews/2026-10-05-<lane>.md` (`2026-10-06-nanse.md` for lane
nanse) has the review; `dev/round-20261005.md` is the round's record.
In brms's own ported suite, bin 1 now passes 389 of 494
assertions, up from 387 at 0.67.0.
