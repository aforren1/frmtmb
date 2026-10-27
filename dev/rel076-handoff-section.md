## The 0.76.0 round, 2026-10-27

`dev/round-20261027.md` is this round's record. Each lane has its
findings (`dev/mocode-findings.md`, `dev/acprior-findings.md`,
`dev/loglik-findings.md`) and its review with the fix-round checks
(`dev/reviews/2026-10-27-<lane>.md`). Three lanes off `70476dc3`
(0.75.0, pushed), merged onto `70476dc3` in this order:

- `mocode`: `mo()` codes an integer from its minimum, as brms does
  (BREAKING since 0.15.0); fractional, missing, factor and out-of-range
  new data refused; `conditional_effects()` and emmeans hold `mo()` as
  brms does; prediction and marginaleffects standard errors include the
  simplex, `NaN` at an empty category; a saved 0.75.0 fit coded the old
  way is refused.
- `acprior`: frmtmb.sample's autocorrelation defaults flat on the
  natural parameter; a prior on class `sigma` under a first-order
  `cov = TRUE` `ar/ma/arma` read at brms's innovation sd; the cosy chart
  on brms's [0, 1]; declared bounds kept under a user prior; one prior
  slot per residual correlation parameter (BREAKING). New in the
  sampling API: `autocor_flat_priors()`, `autocor_bound_chart()`,
  `autocor_innov_sd()`.
- `loglik`: `log_lik()` of draws for residual covariance terms (brms's
  leave-one-out conditional) through `autocor_row_loglik()` and the
  seam slot `frmtmb_register_autocor(jacobian =)`, which `sar()` fills;
  the 13 bin-2 `sar()` and `fcor()` rows ported, 8 of them with the new
  ledger class "adapted-input".

Versions: frmtmb **0.76.0** (minor, breaking); frmtmb.sample
**0.23.0** and frmtmb.spatial **0.3.0** (minor), each needing frmtmb
0.76.0, measured by `dev/rel076-floor.sh`; survival 0.2.0, latent
0.7.1, eam 0.13.0, coupling 0.6.1, ode 0.7.2 and learn 0.7.3 unchanged.
The release library is `C:/Users/adf44/source/r/rellib-r14`;
`rellib-r13` keeps 0.75.0.

### Decisions made at the consolidation, for the coordinator

- **The sampling-API page keeps both lanes' names** (acprior's three,
  then loglik's one); the only textual conflicts outside NEWS.
- **`brmsfit-methods:924` is a divergence by policy**, as mocode's
  review asked: the block's newdata hold `Exp = 4` on every row, and
  brms builds that model only under testmode. Three neighbouring rows
  of the same block (923, 927, 937) keep their divergence outcome but
  now record the refused refit as their observation.
- **The bins 1 and 2 figure is 619 of 823**, not the 620 the brief
  expected: 607 + 13 (the bin-2 port) - 1 (row 924).
- **Two consolidation tests** in frmtmb.sample's `test-loglik-autocor.R`:
  `cosy()` with `mo()` and `ar(cov = TRUE)` through `frm_sample()` and
  `log_lik()`, against the dense conditional of each draw. The `sar()`
  prior tables with both changes are checked by a probe
  (`dev/rel076-overlap.R`), since no test file loads frmtmb.sample and
  frmtmb.spatial together.
- **The floors are measured without an install** (`dev/rel076-floor.sh`
  loads copies of rellib-r14's builds over rellib-r13), because this
  consolidation installs only into rellib-r14.
- **The lanes' dev/ scripts are tracked, the reviewers' are not** in
  the release tree, as at 0.75.0; the loglik reviewer's three output
  folders of ledger tables are left out of its lane commit too.

### The release review and its fixes

`dev/reviews/2026-10-27-release.md`: FIX, no blocker or major. The
coordinator asked for findings 1 to 5 and 7 as release fixes
(`dev/rel076-release-fix2.R`, sourced by the resolve step); finding 6
(three ledger reasons the run no longer reaches) is filed:

1. NEWS lists the refusal of a two-valued integer `mo()`.
2. The record's comparison of the two tmbstan objectives was invalid
   (tmbstan's process-global pointer, `dev/rtmb-pitfalls.md` item 23);
   the reviewer's measurement replaces it, and the backlog item is
   resolved: one ulp in the summation order, amplified by the sampler.
3. `ac_bound_theta()` takes a bound within 64 ulps of a computed end as
   that end, so `ub = 1` on a row-standardized W is accepted, and its
   refusal prints the ends in 15 digits.
4. `?frm_sample` says the sigma intercept's prior under `sigma ~ x` is
   read at the log marginal sd.
5. `log_lik(newdata =)` names `posterior_predict(newdata =)` only where
   it answers (`draws_pp_structured()`).
6. The release commit carries the review and its scripts.

### What the 0.76.0 round is evidence for

**Two tmbstan fits in one process share one objective.** tmbstan's
`sampling()` sets process-global pointers to the fit it runs, so
`log_prob()` of an earlier stanfit evaluates the latest one. The
consolidation compared two runs that way, found them bitwise equal at
800 points and called the cause of their different draws unknown; the
release review evaluated each right after its own run and found the
ulp. A comparison that cannot fail is not evidence: build the case
where the two must differ first.

**A lane's log can be older than its last fix round.** Three of the
lanes' logs that the equality rerun compares against (mocode's
`edge-lane.txt`, acprior's `agree.txt` and `agree-up.txt`) were written
before those lanes' fix rounds and differ from the merged build. Rerun
on each lane's own final library, the merged build's output is the
same line for line (`dev/rel076-equality-lanelib.sh`). Compare with the
lane's final build, not with the newest log of that name.

**An expected figure in a brief is arithmetic, not a measurement.**
The brief added the 13 bin-2 rows to 607 and did not take
off the bin-1 row the same round turns into a divergence. Regenerate,
then report what the ledger says.

**`cd dir && ( ... ) &` runs the cd in the background.** One
consolidation command started the bin-2 pipeline that way and then ran
`dev/ciharden-scan.sh` in the session's own directory, the main
checkout, which made an empty run there (two gitignored files, removed;
record, "The round's process"). Put the cd inside the background group
or use absolute paths.
