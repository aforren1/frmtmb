## The 0.77.0 round, 2026-10-27 to 2026-10-28

`dev/round-20261028.md` is this round's record. Each lane has its
findings (`dev/sigmainn-findings.md`, `dev/mointerop-findings.md`,
`dev/priordraws-findings.md`, `dev/reslot-findings.md`) and its review
with the fix-round checks (`dev/reviews/2026-10-27-<lane>.md`). Four
lanes off `435e7ce6` (0.76.0, committed on main, not pushed), merged
onto `435e7ce6` in this order (reslot before priordraws, which reads
its `(flat)` rows):

- `sigmainn`: under the covariance form of `ar()`, `ma()` and `arma()`,
  first order and higher, `sigma` is brms's innovation sd (BREAKING);
  the likelihood and the simulator read the marginal sd through
  `autocor_marg_sigma()`; 0.76.0's prior layer for the innovation sd is
  removed, the export `autocor_innov_sd()` with it; saved 0.76.0 fits
  are refused where an answer reads sigma.
- `mointerop`: `mo()` on one or two values, as brms (BREAKING);
  marginaleffects `vcov =`, gap-aware standard errors and brms's
  prediction arguments; the rescor bin-2 rows.
- `reslot`: one prior slot per random-effect and `rescor` parameter
  (BREAKING); `prior = "flat"` is brms's `""` in every slot; the
  named-list prior spelling refused; a `route = "fit"` table refused by
  `frm_sample()`.
- `priordraws`: `frm_sample(sample_prior = )`, `prior_draws()`, the
  stored warmup, `update()` of draws, `Evid.Ratio` from prior draws
  (BREAKING), through core's new `prior_only` field of
  `build_objective()`.

Versions: frmtmb **0.77.0** (minor, breaking); frmtmb.sample
**0.24.0** (minor, breaking), needing frmtmb 0.77.0, measured by
`dev/rel077-floor.sh`; spatial 0.3.0, survival 0.2.0, latent 0.7.1,
eam 0.13.0, coupling 0.6.1, ode 0.7.2 and learn 0.7.3 unchanged. The
release library is `C:/Users/adf44/source/r/rellib-r15`; `rellib-r14`
keeps 0.76.0.

### Decisions made at the consolidation, for the coordinator

- **The coerced fit-route table (reslot's review, R1).** Keeping the
  route through `as.data.frame()` would break the pinned rule that the
  coerced table equals a hand-written data frame
  (`test-get-prior-route.R`). Instead `as.brmsprior()` marks a data
  frame that holds frmtmb's printed `(flat)` on a `source = "default"`
  row ("table"), and frm_sample() refuses it where such a row replaces
  one of brms's proper defaults; a coerced route = "sample" table passes.
  It can refuse an intended `(flat)` edit of a coerced sample-route
  table (filed).
- **(flat) on a `mo()` simplex** stopped frm_sample() with "Invalid
  Dirichlet prior ... got 0" (found by the mointerop x reslot probe, on
  reslot's lane too). Fixed as release fix C1: the simplex keeps
  dirichlet(1), brms's reading of `""` there.
- **The 18 rows of brms's own numeric helpers** are divergences by
  policy (the coordinator's decision); the bin-2 ledger reads a policy
  verdict for an unported row from the adapter spec.
- **Dead code removed**: `flat_spec()`'s `innov_sd` store, and the
  `overrides` checks in `prior_summary.frmtmb_draws()` and
  `sample_prior_plan()` (`er_entries()` had gone with priordraws'
  rewrite).
- **frmtmb.spatial keeps 0.3.0**: only its copy of a test helper
  changed, and its tests pass on 0.76.0 with the same counts.

### What the 0.77.0 round is evidence for

**A cross-lane probe finds what neither lane could.** `(flat)` on a
`mo()` simplex failed on reslot's lane already, but no test of either
lane built it; the probe of the pair mointerop x reslot did. Probe each
pair of lanes on the merged build before the tiers run.

**`sh` with Rtools first on PATH segfaults an R child at exit.** The
bin-2 pipeline, started as `sh dev/brmsport-bin2-all.sh` with Rtools'
directories first on PATH, stopped twice on a segfault of
`dev/brmsport-bin2-survey.R` after its output was complete; the same
script exits 0 on its own, and the pipeline passed under `bash` with
the default PATH. Run the pipeline with `bash`, PATH as Git gives it.

### The release review and its fixes

`dev/reviews/2026-10-28-release.md`: FIX, one major. The fixes are
`dev/rel077-release-fix2.R`, sourced by the resolve step:

1. `[.frmtmb_prior_rows` keeps a prior table's route through a row or
   column selection and `subset()`, so frm_sample() still refuses
   frm()'s table cut down (before, `subset()` and `pf[, cols]` sampled
   `(flat)` in place of brms's defaults in silence); a coerced table is
   marked by any printed `(flat)`.
2. `?frmtmb-marginaleffects` and NEWS say the sandwich spellings are
   refused under a residual correlation term.
3. `c()` takes a prior table with its route, as brms's `c()` does.
4. The commit script exports `GIT_OPTIONAL_LOCKS=0`, so its dry run
   writes nothing.
5. Nits: a stale comment, a long line, `qskew_normal()` in the policy
   reason (the ledger regenerated once more: 633 of 823).

**A guard keyed on an attribute is only as good as the operations that
keep it.** R1 and reslot's S1 both tested the route mark through the
edits a user makes on purpose, and neither through `subset()` or a
column selection, which drop attributes in base R. When a guard rides on
an attribute, list the base generics that copy the object and test each.
