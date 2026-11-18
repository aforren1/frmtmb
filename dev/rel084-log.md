# Running log of the 0.84.0 consolidation (round of 2026-11-16)

Release tree `C:/Users/adf44/source/r/frmtmb-wt-rel084` (branch `wt-rel084`, off
`b255350e`), library `C:/Users/adf44/source/r/rellib-r23`, base `rellib-r22` (0.83.0).
Machine output in `dev/rel084-log/` (gitignored). The record is `dev/round-20261116.md`.

## Lane records read (2026-11-17)

Per-lane digests by four read-only subagents, `dev/rel084-log/digest-<lane>.md`
(gitignored). The last verdict on disk in each lane tree's `dev/reviews/`:

- SHIP on file: ordstarts (`ordstarts-short.md`), sparsese, drawdir, gpselimit
  (`<lane>-short.md`).
- FIX on file, a fix round after it, and no later short check in the lane tree,
  the release tree or `round-briefs/`: cdfname (review only), alapse
  (`alapse-short.md`, then fix round 2), nbrobust (review only; records-only fix
  round), slopetrig (review only; fix round changed `R/autoscale.R`), squant
  (`squant-short.md` "FIX (one item)", then fix round 2: compat rows and docs),
  nlband (`nlband-short.md`, then fix round 2: `R/fit-restart.R`,
  `R/par-template.R`), quadesc (review only; fix round changed `R/fit.R`).
- The brief says all eleven are SHIP; the coordinator's round notes say the same.
  Taken as the brief says; the missing short-check files are reported.

## Merge (2026-11-17)

`bash dev/rel084-merge.sh` (`dev/rel084-log/merge/`): sparsese 48 paths, gpselimit
56, slopetrig 41, alapse 60 (1 conflict), nlband 65 (1), quadesc 37, nbrobust 40,
squant 65 (1), drawdir 47, ordstarts 103, cdfname 71, samplers 29, codesize 15,
cran 11. Skipped as extracted copies: `dev/reviews/sparsese-review-base/` (6),
`dev/reviews/ordstarts-build-review-base/` (10). No stan-cache program was missing.

Conflicts, resolved by hand (Edit):
- `R/se-check.R` `se_lost_clauses()` order vector (alapse): union, slopetrig's
  `"shortcor",` line kept, then `"nonfinite", "nocurv", ald_se_reason)) {`.
  gpselimit's `"gpwhite"` merged cleanly on the `"boundary"` line.
- `R/par-template.R` `nl_fit_probed()` pole path (nlband): both comments, filter
  `c("warning", "frmtmb_nl_local", "frmtmb_nl_restart", "frmtmb_boundary_fit",
  "frmtmb_gp_basis")` (merge step 2).
- `R/compat.R` (squant): alapse's restated `gp_pred()` x asym_laplace note, then
  squant's `elf_note` block and its 12 rows (merge step 6).

Merge steps applied (2026-11-17):
- Step 5 (alapse x squant): `quantile_elf()` first in `ald_se_routes`
  (`R/ald-fixed.R`, squant's wording, which names `?quantile_elf`), a sentence
  in `ald_laplace_message()` (`R/ald-laplace.R`, with "use smoothing = 1"),
  `?frmtmb-families` "Quantile regression inference" names `[quantile_elf()]`
  first; NEWS (`dev/rel084-news.R`). The `@noRd` comment of `ald_se_routes` no
  longer names a lane.
- Step 1 (cdfname): `Rscript dev/cdfname-rename.R --dry-run`: "0 line(s) in 1025
  candidate file(s)"; the run: "RENAME done: 0 line(s)". Acceptance grep
  (`grep -rIn '\blcdf\b' --exclude-dir=.git --exclude-dir=dev .`): NEWS.md 8,
  R/families.R 1, test-cdf-slot-name.R 13, eam NEWS 7, survival NEWS 2.
- NEWS (`dev/rel084-news.R`): core 24 bullets, 15 BREAKING, 7 sections; eam
  0.14.0, survival 0.4.0, coupling 0.6.2, one bullet each, "Needs frmtmb 0.84.0".
- Backlog: the merged heading renamed "## Filed at the 0.84.0 release
  (2026-11-16)"; 10 lane subsections (cdfname filed none); 4 HIGH, 14 MEDIUM,
  44 LOW; 18 CLOSED lines of the round (3+1+2+1+2+2+3+1+2+1). Upstream bugs:
  no id collides; R (stats) R-1, RTMB-3 "Also", RTMB-4 measurements, samplers'
  four entries present.
- Versions (`dev/rel084-release-fix.sh`): frmtmb 0.84.0; eam 0.14.0 and
  survival 0.4.0, floor `frmtmb (>= 0.84.0)`; coupling 0.6.2 (floor 0.61.0
  kept); run-tests.R and run-session.R name rellib-r23. codemeta.json
  regenerated (`codemetar::write_codemeta()`, which read the README URL and
  topics from the GitHub API, as at 0.83.0).
- Install (`dev/rel084-install.sh`, `install-1.txt`): frmtmb 0.84.0, eam 0.14.0,
  latent 0.7.5, ode 0.7.4, survival 0.4.0, spatial 0.5.0, coupling 0.6.2,
  sample 0.30.0, learn 0.7.4, drmTMB 0.7.0, all in rellib-r23.

Focus tier (40 changed test files, every gate on, `dev/rel084-verify.sh focus 3`):
rellib-r23 ran 40 of 40, pass 6733, fail 0, err 0, warn 0; 40 lib: lines name
rellib-r23. rellib-r22 (`base 2`): 24 of 40 files fail or err, 351 failures and
65 errors (pass 5899).

Key tables: a subagent copied the lanes' scripts to `dev/rel084-kt/` (rellib-r23
alone) and ran cdfname and alapse (MATCH), then handed back early with the
ordstarts driver still running; `dev/rel084-kt/run-queue.sh` runs the rest one R
process at a time, then merge step 6 (`compat-ext.R none|all`). Its write of a
summary file was refused by the tool; its report is in `dev/round-20261116.md`.
`bash dev/rel084-chain.sh` started 2026-11-17 (2 R processes; phases in its header).

## Key tables on rellib-r23 (2026-11-17, `dev/rel084-kt/`, logs `dev/rel084-log/kt/`)

Compared by a read-only subagent against the lanes' logs and findings; every log's
lib line names rellib-r23.
- cdfname MATCH (route probe 7 routes, rds-b, eam probe). alapse MATCH (rq table
  11 lines and 179 logLiks identical; design B 46 seeds, max gap 8.86e-11; D/E
  silent 9/9; bounded, over, ties 22/22; short-ties 160/160, 0 of 365, 37/60, 60/60).
- ordstarts MATCH: sum6 "silent shortfalls before 19 after 3 of 320 ; D down from
  a proper code-0 A0: 0"; A0 313 of 313 equal 0.83.0; weights "pass 21 of 21";
  short-p9 byte-identical. Costs moved slightly (D/A0 sum 1.79 vs lane 1.85).
- nbrobust MATCH (reference BLAS; jacobian bands, 3805 of 3888 max rel 4.610e-12;
  nb200 21/21 at 17.2212, nb300 42/42 at 18.9741; start at 30 21/21).
- slopetrig MATCH (5033: fits 50 engaged 50 short 0; lme4 0 of 20; gr(by =)
  engaged 10 of 10, D short 7; calibration counts identical; cor and guard).
- squant MATCH (qgam, ulps, tape gradients identical; HELD 21/21/21).
- nlband MATCH (fix1, studg, sweep identical; survey "ALL rows 2555 differ 0").
- sparsese MATCH (3370ref s1 to s3; fix1 and cum4500 identical).
- drawdir: repro and fires12 identical (11 warn, slpw:12 reweights, Pareto k 0.8);
  survey 286 fits, 0 summed fires, per-direction fires as the lane. mix40
  cumulative MOVED: no verdict 11 -> 18 of 40. `dev/rel084-kt/mix40-warn.R`:
  ordstarts' swapped start reaches higher maxima (seed 8 logLik -375.940 vs
  -378.694 on rellib-r22; seed 10 -396.461 vs -399.826; seed 27 -375.430 vs
  -378.658), all code 0, where predict() says the covariance could not be
  recovered from the Hessian (frmtmb_warning) instead of a draw verdict.
  Explained by ordstarts; not a draw-check change.
- gpselimit MATCH (30 fits, 13 white rows message, plateau, gp(by =) NaN 4).
- quadesc MATCH on the reference arm (63 of 189; alphax -85.1020..-84.9742; veto
  PASS 21 of 21, alpha 3.44123, gain 11.309, evals 43). OpenBLAS arms of quadesc
  and nbrobust not rerun here (the chain's OpenBLAS phase covers the test files).
- Step 8 probe (13 white gp fits): the draw check fires on none, no reweight, no
  warning; largest dir 2.719. Step 4 probe: 50 D/E fits, all frmtmb_ald_local and
  frmtmb_ald_no_se, 0 with a gp, boundary, limit or converged_at_max class.
- Step 6: test-compat.R with no extension pass 878, with all eight pass 1158,
  fail 0, err 0, warn 0.
- Plain tier (chain phase 1): ran 553 of 553, pass 36271, fail 0, err 0, warn 0,
  skip 293; 553 lib: lines rellib-r23.
- Gated tier (chain phase 2, every gate on): ran 553 of 553, pass 40090, fail 0,
  err 0, warn 0, skip 0; 553 lib: lines rellib-r23. Per suite (files, plain pass /
  skip, gated pass): frmtmb 353, 26730/243, 29975; coupling 11, 551/5, 557; eam
  33, 1855/5, 1879; latent 13, 700/2, 706; learn 16, 474/13, 559; ode 11, 549/1,
  551; sample 92, 4286/18, 4654; spatial 8, 551/6, 634; survival 16, 575/0, 575.
- Ported brms suite (chain phase 3): bin 1 14 files, ledger, gen, guards exit 0;
  `dev/rel084-ledger-diff.R`: pass 398, cannot transfer 48, divergence 48 of 494,
  "Rows whose outcome, class or reason moved: 0". Bin 2: 636 of 823 (77.3%), 636 of
  723 (88.0%); `dev/rel084-bin2-diff.R`: 0 rows, 0 classes, 0 value checks moved.
  Regenerated files that moved in text only: `brmsport-bin2-groups.tsv` (G1 names
  `$cdf`, cdfname's rename of the bin-2 scripts), `brmsport-bin2-valuechecks.tsv`
  (3 HOLDS rows' max diffs, e.g. log_lik:417 5.77e-15 -> 3.55e-15: nbrobust's
  density), and the raw record `dev/brmsport-log/rec-frmtmb-methods.tsv` row
  brmsfit-methods:169 (divergence -> cannot transfer, a STALE spaghetti row; the
  ledger keeps it as adjudicated at 0.83.0, row 169).
- Docs check (chain phase 3): roxygen again over nine packages, 0 files moved;
  installed code against the tree: 0 functions differ in all nine (frmtmb 2707
  compared, sample 335, eam 194, latent 67, survival 51, learn 47, spatial 46, ode
  38, coupling 24); check_pkgdown OK for all nine; build-docs 9 of 9 sites.
- Survey (chain phase 4, `dev/rel084-survey.sh all`, RECORD=1, rellib-r23;rellib-r22):
  172 pairs, 172 lib: lines rellib-r23 (frmtmb 0.84.0); per-pair wall sum 892 s.
  Headline **89.3, Wald [84.8, 93.9]** (bootstrap [84.5, 93.5]), no defect 93.0
  [89.2, 96.8]: both as 0.83.0. `dev/rel084-equiv-delta.R`: 0 verdicts changed
  (cannot measure 15, divergence by decision 7, equivalent 131, refusal correct 19,
  as before); l1_const moved on 0 pairs; l1_grad on 1, p131 (cs_pred() x mixture,
  25.870 -> 23.637, l1_const 0, equivalent L1). Compat frames
  (`frame-diff-r22-r23.txt`): 12012 -> 12116 pairs, 104 only after (quantile_elf);
  3 status moves, asym_laplace x {homcs, homtoep, toep} works -> conditional
  (alapse); 20 note-only changes (alapse's restated notes, cdfname's `cdf`).
- Floors (chain phase 5, `dev/rel084-floor.sh 2`; each extension as built in rellib-r23
  over rellib-r22's 0.83.0 core): eam (floor lowered to load) 5 of 33 files fail or err
  (contaminant 3 fail, ndt-bound 1 err, rdm-gng 1 fail 2 err, surface 1 fail,
  wiener-cdf 9 fail): needs 0.84.0. survival (lowered) 9 of 16 files err or fail (cox
  6 err, cox-bhaz 6, rp-floored 5 fail 12 err, ...): needs 0.84.0. coupling (floor
  0.61.0 kept) 11 of 11 pass; sample (floor 0.83.0 kept) 92 of 92 pass.
- Translator audit (chain phase 6, `dev/rel084-audit.sh`): 132 pairs (the record's L1
  pairs and p069; p072 left L1 at 0.83.0), 132 with a result, 66 with latent blocks,
  130 blocks, 2 untestable, 0 BLIND, 0 seen by less than 100 x tol; 132 lib: lines
  rellib-r23.
- R CMD check --as-cran (chain phase 7, `dev/rel084-check.sh`, R_LIBS rellib-r23):
  frmtmb 0.84.0 1 NOTE, frmtmb.eam 0.14.0 1 NOTE, frmtmb.survival 0.4.0 1 NOTE (each
  "Skipping checking math rendering: package 'V8' unavailable"), frmtmb.coupling 0.6.2
  OK. Logs copied to `dev/rel084-log/checklogs/`.
- OpenBLAS (chain phase 8, `dev/rel084-ob.sh`, 23 files of `dev/rel084-log/ob.jobs`,
  every gate on, 4 threads): 0.3.26 and 0.3.32 with LAPACK, `OPENBLAS_CORETYPE`
  unset and Haswell: each ran 23 of 23, every file fail 0, err 0, skip 0, warn 0;
  23 lib: lines rellib-r23 per run; blas lines name Rob0.3.26-lapack and
  Rob0.3.32-lapack. The R copies were deleted after (`dev/ciharden-out/` gone).
- Sessions plain (chain phase 9, `run-sessions.sh rel0.84.0 2`): 27 of 27, every SESSION
  line fail 0 err 0 warn 0; 243 of 243 lib: lines rellib-r23; LEAK 70 (46 package
  load, 24 extension session state, 0 untagged). frmtmb 26730, 26898, 27004 (forward,
  reverse, preload); sample 4286; coupling 551, eam 1855, latent 700, learn 474, ode
  549, survival 575, spatial 551.
- Sessions gated (chain phase 10, to 05:11): 27 of 27, every SESSION line fail 0 err 0
  warn 0; 243 of 243 lib: lines rellib-r23; LEAK 221 (194 package load, 142 Stan model
  loads; 27 extension session state; 0 untagged). frmtmb 29975, 30143, 30249; sample 4654.
- CHAIN END 2026-11-18 05:11. Scratch removed: `dev/rel084-check/` (logs kept in
  `dev/rel084-log/checklogs/`), `dev/docsci-site-rel084/`, `dev/ciharden-out/` (by the
  OpenBLAS script), the floors' scratch library (by its script), `dev/rel084-progress.md`.

## After the release review (2026-11-18, `dev/reviews/2026-11-18-rel084.md`, FIX, four text items)

- F1: the mix40 attribution corrected in `dev/round-20261116.md` (Semantic conflicts,
  Key tables drawdir, What is not done) and `dev/round-handoff.md`: C0, not the swap.
  The entry above under "Key tables" (ordstarts' swapped start) is superseded.
- F2: `R/families-quantile.R` (so `man/quantile_elf.Rd`) and the `gp_pred()` x
  quantile_elf note in `R/compat.R` cite `dev/reviews/squant-short-gp.R` and
  `squant-short-gpprof.R`; both added to CONS_EDITS of `dev/rel084-commit.sh`.
- F3 and F4 b (NEWS, through `dev/rel084-news.R`): "Before," for the two
  negative-binomial bullets; "another link of `mu`"; the bound sentence qualified to a
  maximum inside the bounds.
- F4: step 4 of the record carries the review's measurement; the backlog's 0.84.0
  section gains "### From the 0.84.0 release review" (1 LOW, 1 MEDIUM).
- `bash dev/rel084-review-fix.sh` (`dev/rel084-log/review-fix.txt`): roxygen, checkRd,
  core install, three test files, test-compat.R both arms, check_pkgdown and build-docs
  for core.
