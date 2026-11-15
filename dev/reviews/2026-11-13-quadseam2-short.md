# Short check of lane quadseam2's fix round (round of 2026-11-13, base ec8178d3)

Checks the fixes for F1, F2 and the minor item of
`dev/reviews/2026-11-13-quadseam2.md` and the "Fix round" section of
`dev/quadseam2-findings.md`. It does not review the lane again. Scripts
`dev/reviews/quadseam2-short-*`, logs `dev/reviews/quadseam2-short-log/`
(gitignored). Reference-BLAS runs are pinned to one P-core mask
(`quadseam2-review-aff.cmd`). At most 2 R processes at a time. File
compares used `git show ec8178d3:<path>`, `git hash-object` and `diff`
on copies.

## Verdict: SHIP

Build: the lane library's installed code is the fix-round source
(`quadseam2-short-build.R`, `build.txt`): `R/fit.R functions 93 differ
0`, `R/autocor-latent.R functions 57 differ 0`, `R/se-check.R functions
84 differ 0`. frmtmb.latent loads from rellib-r21. Its `R/` and
DESCRIPTION hash the same as ec8178d3, and only its test file changed.

## 1. F1: the rule's escape (R/fit.R). Fixed

The escape call site in `fit_assembled()` is the same as at ec8178d3,
apart from a 4-line comment: `esc_lap`, an unconditional
`escape_stationary()`, then `stationary_escape_keep(opt, esc_lap)` when
the records differ. The body of `stationary_escape_keep()` is identical
to ec8178d3. The quadrature sd restart fills `stationary_escape` and
`stationary_escape_laplace` from `fit$opt` when the restart's opt lacks
them (`R/fit.R:1931-1933`).

esc2 on the lane library (`esc2-lib.txt`) joined with the reviewer's
`esc2-base.txt`:

    fits 45 lane below base by >1e-8: 0 differ beyond 1e-10: 0
    ALPHAX/11   rule 6.667 ll -84.5783086103 alpha 951780,-2644632   (= base)
    SYMSMALL/1  rule 2.291 ll -84.9253951576 alpha 9359630           (= base)
    WEAK/4      rule 1.453 ll -74.1377898690 alpha 3409589           (= base)

5222 gains kept (`dev/quadseam2-hmm.R` via `quadseam2-short-hmm.R`,
`hmm-lib.txt`): HMM seeds 1 to 3 and MIX seeds 1 to 3 have stall =
free to 4 decimals (-998.2551, -960.6004, -953.5887, -635.4970,
-579.4192, -614.0458), with escape starts 6 on each. These are the
same numbers as the lane's `hmm-lib-fr.txt`.

Differences that remain in the escape area against ec8178d3, each
justified:
- `escape_stationary()`: the `repeat` loop. It runs again from the kept
  optimum while it gains, for at most `n_on` passes. The record's
  `starts` value is the total over the passes. This is item 5222, and
  without it the gains above are lost. The loop ends because
  `max_passes` is fixed from the first call, and each pass is finite.
- `stationary_escapes()`: `attr(out, "n_on")`, the count of predictors
  on their point. The loop above needs it, and it does not change the
  starts list.
- The `quad_sd_up_escape()` block after the escape (item 4861), with the
  copy of the records (the third point of F1). The restart starts from
  the kept parameters with only the small sds moved up, so an escaped
  alpha stays off its point.
- `sd_up_escape()` is split into `sd_up_start()` plus the restart, so
  the quadrature route can share the probe. On the Laplace route the
  behavior is the same: the reviewer's 18 `ident` fits were identical.
- Merge note, kept for the consolidation: priors2's
  `check_map_unbounded()` must also be called on `esc$opt` in the
  quadrature restart branch. It is in the findings.

## 2. New MEDIUM filing (OpenBLAS). Not a regression

`dev/ciharden-openblas.sh 0.3.26 lapack` (`ob-build.txt`: Rblas 74 to
OpenBLAS, Rlapack 616), `OPENBLAS_CORETYPE=Haswell`, both arms pinned to
mask 0C00 (`quadseam2-short-ob.cmd`, `dev/quadseam2-symsmall.R`,
SYMSMALL seed 1):

    base  ESC laplace starts 2 in 87.216780 out 86.497099 alpha_out 6.727
          ESC rule    starts 2 in 87.216780 out 87.216780 alpha_in 7.668e-15
          SEED 1 ll -87.2167801711 alpha 7.66803e-15 warn frmtmb_warning
    lane  (the same three lines, character for character)
    base, second process: the same

The rule's optimum walks back to alpha = 0 after the Laplace escape. The
rule's escape runs on both builds, and both of its starts come back.
Base and lane are bitwise identical, so 0.82.0 already has this
behavior. The filing is correct as MEDIUM. `dev/ciharden-out/` was
deleted after the run.

## 3. F2: latres_quadratic() (R/autocor-latent.R). Fixed

Code: `h <- pmax(h, 1)`. Rows not finite at the fit are left out
(`read`). If any read row has a non-finite third or fourth difference,
the answer is FALSE. `suppressWarnings()` encloses only `row_lpdf()`.

Attack (`quadseam2-short-latres.R`, base and lane, 33 fits, `ar(t, g,
cov = TRUE)`, 25 groups of 6, seeds 901 to 903). The families have a
support where a step of 1 leaves the support or is not quadratic: Beta
logit with a mean near 0.03 and near 0.97, Beta identity at 0.5, Beta
log at 0.9, von_mises tan_half, von_mises identity near 2.6, Gamma
inverse at eta 2:

    gate FALSE on 21 of 21; gate TRUE on 12 of 12 lognormal; no warning
    from inside the gate (0 "gatewarn")
    rows 33 differ (ll, latres count, other warnings) 0   (base vs lane)

7 of the 21 bounded-support fits warn `frmtmb_latres_laplace` on both
builds (BETAlog seeds 1 to 3, BETAid 2, VMtan 2 and 3, VMid 3).
On FALSE, the caller falls through to the shift test of 0.82.0. Thus
FALSE gives 0.82.0's answer. It can keep a spurious warning that 0.82.0
also gave, and it can never drop a real warning. The default is correct.
A family whose support on the link scale is bounded cannot have a log
density that is quadratic over the full real line, so TRUE is never the
correct answer for it.

FEWBERN (`quadseam2-review-fewbern.R`, `fewbern-lib.txt`): seed 703
(marginal sd 7.557e-05) warns ("no curvature to read its first-order
correction by"), gate FALSE. Seeds 701, 702 and 704 warn, gate FALSE.
The step sweep reads FALSE at each value, because the override now goes
through `pmax(h, 1)`.

Lognormal at y x 1e-8 and 1e8, seeds 40401 to 40403, `mi()` and
weight-0 spellings: gate TRUE and latres 0 on 12 of 12. The logLik is
the same for both spellings, on base and on the lane.

## 4. NEWS and backlog. Correct

The new section has 4 bullets, each BREAKING. The "escaped once" bullet
is gone. The latres bullet is BREAKING. A search of the section for
version numbers, `dev/`, lane names and em dashes finds nothing. The
only "once per" is in the 5222 bullet, and it is correct there.

The 5408 CLOSED line says "not a defect" and gives the three fits as the
reason. It also says that the second escape stays and returns at once
on a fit off its point. This is accurate: `escape_stationary()` returns
at once when `stationary_escapes()` is empty. The reviewer suggested
that the line be removed. This wording records the result instead, and
I accept it.

## 5. Tests. Correct

- `test-stationary-passes.R` accepts `starts` of 4 or 6. Both are valid
  endings. With `max_passes = 2`, 4 means one pass that gained and left
  no predictor on its point. 6 means a second pass of 2 starts. A pass
  that gained nothing also gives 4, but the test then fails
  `all(abs(a) > 1)` and the logLik check, so the wrong ending is still
  caught.
- `test-quad-escape-record.R`: the record test is back as at ec8178d3,
  and there is a regression test on ALPHAX seed 11. The data draw
  matches esc2 (seed 1011, the two earlier designs' draws are made
  first). It asserts a Laplace record, a rule gain > 1e-3 x |logLik|,
  and `all(abs(alpha) > 1)`. On the first-round build, the reviewer's
  `esc2-lib.txt` shows ALPHAX/11 alpha `4.801e-16,-5.186e-15`, which
  fails the alpha assertion. The lane's `fr2-qer-firstround.txt` has
  pass=7 fail=2.
- LOW, optional: the regression test catches warnings with its own
  handler and accepts `class(w)[1]` of `frmtmb_conv_code` or
  `frmtmb_warning`. The round rule asks for `allow_warnings()` by
  message. The generic `frmtmb_warning` class admits any frmtmb warning
  that has no subclass, not only the standard-error one. A match on the
  standard-error message is narrower. This does not block the ship.

## 6. Test runs on the lane library (reference BLAS, one file per process)

`quadseam2-short-tests.sh`. Each `lib:` line shows frmtmb 0.82.0 from
`wt-quadseam2-lib`:

| file | pass | fail | err | skip | warn |
|---|---|---|---|---|---|
| test-quad-escape-record.R | 9 | 0 | 0 | 0 | 0 |
| test-stationary-passes.R | 4 | 0 | 0 | 0 | 0 |
| test-latres-quadratic.R | 6 | 0 | 0 | 0 | 0 |
| test-latres.R | 305 | 0 | 0 | 0 | 0 |
| test-quadrature-defects.R | 89 | 0 | 0 | 0 | 0 |
| test-se-tier.R | 114 | 0 | 0 | 0 | 0 |
| frmtmb.latent test-hmm-stationary.R | 12 | 0 | 0 | 0 | 0 |

## Housekeeping

Free space on C: was 146.06 GB at the start and 146.07 GB at the end.
The OpenBLAS build was made and deleted. No run failed on disk or allocation. No git writes.
No lane files were edited. The new files are `quadseam2-short-*` and
this record.
