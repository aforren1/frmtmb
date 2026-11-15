# Running log of the 0.83.0 consolidation (round of 2026-11-13)

Release tree `C:/Users/adf44/source/r/frmtmb-wt-rel083` (branch `wt-rel083`, off
`8159568e`, 0.82.1), library `C:/Users/adf44/source/r/rellib-r22`, base `rellib-r21`
(0.82.0). Machine output in `dev/rel083-log/` and `dev/ciharden-log/` (gitignored).
The record is `dev/round-20261113.md`.

## Lane records read (2026-11-14)

The last verdict of each lane's review and short check, read against the brief's "all
fourteen: SHIP" (per-lane digests by five read-only subagents, checked against the
files' modification times):

- SHIP, no edit after it: famend, refits.
- SHIP, doc or record edits after it that the check asked for: quadseam2 (a test now
  uses `allow_warnings()`, 8 expectations; a backlog remark), slopescale (roxygen and
  NEWS wording of the autoscale trigger), nlstop2 ("FIX, record only", the record edits
  made; NEWS and `?huber` wording), gpfield (n1 to n3 wording in NEWS and
  `?predict.frmtmb_fit`).
- FIX, answered by a fix round no reviewer has read: selimit2 (F4b: predictions of
  whitened REML fits above `joint_root_max` get NaN standard errors, R/se-check.R,
  R/predict.R), priors2 (fix round 3: `open_end_probe()` and the warning
  `frmtmb_prior_unbounded_local`, R/fit-open-ends.R), silent2 (F2: factor responses
  compared by label, R/frame.R, R/confint.R), arspat2 (N1 to N3 NEWS; the margin slot
  checked in full, R/autocor-seam.R), parity2 (F3a: `ce_grids_build(re_formula = )`,
  R/conditional-effects.R), wcheck2 (S1: no refusal beside a Student-t or another
  package's block, R/ald-laplace.R), predraw2 (dilution limit documented; the no-loo
  allowance conditional), hygiene (fix round 3: a named environment is never merged,
  R/registries.R).

Paths the merge did not take (untracked, not ignored, machine output of a reviewer):
`.github/scripts/__pycache__/coverage-plan.cpython-314.pyc` and
`dev/reviews/hygiene-review-cm-base.json` (lane hygiene).

## Done (2026-11-14)

- Merge (`dev/rel083-merge.sh`, `dev/rel083-log/merge/`): selimit2 66 paths, quadseam2
  62, slopescale 34, nlstop2 51, priors2 103, famend 82 (1 conflict), silent2 57 (1 file,
  3 conflicts), arspat2 66, refits 85, parity2 85, wcheck2 43, gpfield 84 (1 conflict),
  predraw2 52, hygiene 72; 2 skipped.
- Conflicts: `dev/rel083-fix-at0.R` (merge step 2), `dev/rel083-fix-families.R` (step 1),
  `dev/rel083-fix-fit.R` (refits x gpfield, beyond the brief's list). Step 3:
  `dev/rel083-fix-quadmap.R`. Step 5: `dev/rel083-nl-starts.R`,
  `tests/testthat/test-nl-starts-control.R`.
- NEWS (`dev/rel083-news.R`): core 89 lane blocks in 18 sections plus the `nl_starts`
  bullet, 50 BREAKING; frmtmb.sample 8 bullets, frmtmb.spatial 2, frmtmb.latent a "Needs"
  section. Backlog: the one merged heading renamed "Filed at the 0.83.0 release
  (2026-11-13)"; 108 CLOSED lines of the round kept (the lanes' 6+1+7+10+4+7+8+6+3+14+12
  +10+4+16). dev/upstream-bugs.md: no id collides.
- Versions (`dev/rel083-release-fix.sh`), codemeta.json regenerated
  (`codemetar::write_codemeta()`, which read the repository's README URL and topics from
  the GitHub API), install (`dev/rel083-install.sh`, `install-1.txt`).
- Focus tier, 93 changed test files, every gate on: rellib-r22 fail 0, err 0 (one file
  with 2 warnings, my new test, fixed and rerun: pass 21); rellib-r21: 72 files fail or
  err (513 failures, 110 errors).
- Key tables: subagent, `dev/rel083-log/kt/summary.md`, scripts `dev/rel083-kt/`.
- `bash dev/rel083-chain.sh` from 11:23 (phases in its header; `dev/rel083-log/chain.txt`).

## After the chain (2026-11-14, 19:35 to 20:40)

- Chain phases 1 to 8 (11:23 to 19:35): plain 549 of 549 and gated 549 of 549, fail 0,
  err 0, warn 0; docs check clean (0 roxygen moves, 0 functions differ, check_pkgdown OK,
  9 of 9 sites); survey 172 of 172; floors (each needed); audit 133 pairs, 0 BLIND; R CMD
  check frmtmb, sample, spatial, latent each 1 NOTE (V8); OpenBLAS 0.3.26 and 0.3.32,
  CORETYPE Haswell and unset, 22 of 22 clean in each; sessions 27 of 27 in both tiers.
- Step 5 wording (coordinator): `dev/rel083-nl-starts-doc.R`, `dev/rel083-news.R` rerun,
  core roxygen (`man/frmtmb_control.Rd` only), core reinstalled (`install-2.txt`);
  `test-nl-starts-control.R` 21 pass; checkRd 0; check_pkgdown OK; installed code 0
  differ.
- Survey: p072's measured verdict moved, but its 0.80.0 adjudication row held it at
  "equivalent"; `dev/rel083-p072.R` makes it a divergence by decision, and
  `dev/equiv-summary.R dev/rel083-equiv-out --record` regenerated the record: headline
  89.3 [84.8, 93.9].
- Ledger diffs: 0 rows moved in bin 1 and bin 2.
- SE comparison by test identity (`dev/rel083-se-*.{R,sh}`, `se-cmp.txt`): fired 185,
  already NaN 179, newly NaN 6, no base 0.
- One brief overlap: the R call that inserted the round-handoff section ran while the
  SE run held 3 R processes (4 for a few seconds).
- Commit script dry run (`commit-dry-run.txt`). Its first two dry runs stopped on the
  drift guard: `dev/release/run-session.R` (the rellib-r22 line, now a consolidation
  edit) and two multi-lane pairs: quadseam2's `R/fit.R` (merge step 3 sits next to its
  lines; exempt) and selimit2's `R/se-check.R` (its hunks beside quadseam2's and
  slopescale's; a replay of the merge gives the release copy, and the guard now replays).
- Scratch removed: `dev/rel083-check/` (logs kept in `dev/rel083-log/checklogs/`),
  `dev/docsci-site-rel083/`, `dev/ciharden-out/` (by the OpenBLAS script), the floors'
  scratch library.

## After the release review (2026-11-14, 21:00 to 21:45)

F1 to F3 and the test-sugar.R seed fixed (`dev/round-20261113.md`, "After the release
review"); core roxygenised and reinstalled; `dev/rel083-fix-rerun.sh` clean; ledger 1
row moved (row 169, re-adjudicated); dry run clean, rel083 998 paths, 16 index sha256
unchanged. Nothing running.
