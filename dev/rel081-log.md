# Consolidation log, round of 2026-11-09 (frmtmb 0.81.0)

Running log of the consolidator. Brief:
`C:/Users/adf44/source/r/round-briefs/rel081-brief.md`. Release tree
`frmtmb-wt-rel081` (wt-rel081, off e0cce906); library `rellib-r20`.

## 2026-11-10

- Start. Release tree clean on e0cce906, stan-cache 535 programs, rellib-r20
  empty. C: 38 GB free.
- Lane trees all on e0cce906 on wt-r81-<lane>.
- Merge: `bash dev/rel081-merge.sh`, order translator, compat, code8,
  hollow, support (the brief's). 0 conflicts in any lane (33, 55, 86, 48,
  51 paths). 66 paths skipped as machine output (review-*-out*,
  compat-m-race*). The first run wrote LF: Git Bash's grep drops CR, so
  the CRLF test read none. All 71 modified files were rewritten CRLF after
  (`git ls-files --eol -m`: 71 i/lf w/crlf), and the script now tests
  with tr.
- Release fixes: `Rscript dev/rel081-release-fix.R` (skip_on_cran in the
  GEV importance/quadrature test, support record names von_mises `ub`,
  versions frmtmb 0.81.0, sample 0.28.0 with floor 0.81.0, survival
  0.2.3, run-tests.R/run-session.R on rellib-r20, adjudication edits: 7
  rows changed, ids now p002 p021 p072 p069).
- Install: `bash dev/rel081-install.sh <rellib-r20>` (log
  `dev/rel081-log/install-1.txt`): roxygen of core and eight extensions,
  all ten packages in rellib-r20 (drmTMB copied from r19). Roxygen added
  no diff beyond the lanes' own man pages.
- Merge check loo.R (hollow x support): `dev/rel081-loo-exact.R` (log
  `dev/rel081-log/loo-exact-r20.txt`). GEV `y ~ x`, seed 349, 200 draws;
  three draws moved so row 21 is off their support (t -0.19, -0.065,
  -0.029). `log_lik(draws)`, `log_lik(draws, newdata = draws$data)` and
  `newdata = data` are identical, -Inf in the moved cells (62 cells -Inf,
  0 other non-finite); `row_lpdf(exact = FALSE)` there is -2.1e82,
  -2.8e251, -1.0e304. `loo()` both ways stops with the loo package's "All
  input values must be finite or -Inf" (same message): a -Inf log_lik
  cell is +Inf in its ratios. Observation for the record, not a conflict.
- Pre-survey pairs (`bash dev/rel081-pre.sh`, log `dev/rel081-log/pre.txt`,
  1 pair at a time): p021, p047, p048, p131, p136 (translator) and p154,
  p010, p146, p167, p169 (compat) on rellib-r20: every verdict, level,
  `l1_const`, `l1_grad` and every shared column identical to the lane's
  own outdir (translator-noreg3-out, compat-noreg2-edit-out). 10 of 10
  `lib:` lines rellib-r20 0.81.0.
- NEWS: `Rscript dev/rel081-news.R` (rerunnable: released sections read
  from e0cce906). One `# frmtmb 0.81.0` section, 19 bullets: Convergence
  and standard errors (5), Families (3), Hypothesis tests (2),
  Compatibility table (1), Bug fixes (2: stage2, and make_conditions() of
  draws, which hollow's NEWS did not name), Limits (6). frmtmb.sample
  0.28.0 (5 bullets, floor line), frmtmb.survival 0.2.3 (1). The von_mises
  `ub` case named; the families bound bullet marked BREAKING (support's
  record classes it so).
- Backlog: `Rscript dev/rel081-backlog.R`: one section "Filed at the
  0.81.0 release (2026-11-10)", the three lanes' items as written (4),
  then 5 items from the reviews, the lanes' records and this
  consolidation (live-row count line, reader floors, translator gaps,
  code-8 message on ubuntu needs CI, loo() on a -Inf cell).
- upstream-bugs: the translator's brms-49 collides with the existing
  brms-49 (predictor.bprepnl) and brms-50 exists: renumbered brms-51 by
  `dev/rel081-release-fix2.R` in upstream-bugs.md, adjudication.tsv
  (p021's reason) and translator-findings.md (with a note); the review
  keeps its text.
- Test chain started 07:44 (`bash dev/rel081-chain.sh`: focus on r20,
  base on r19, plain, ref; 3 processes).
- Focus tier on rellib-r20 (25 changed or new test files, every gate;
  `dev/ciharden-log/r081-focus.sum`): 25 of 25 ran, 0 fail/err/warn but
  the expected `test-brms-suite-methods.R:419` (brmsfit-methods:326 "now
  HOLDS but is recorded as 'cannot transfer'"), which the ledger
  regeneration clears. 25 of 25 `lib:` lines rellib-r20.
- Base arm, the same files on rellib-r19 (`r081-focus-r19.sum`): 12 of 25
  files fail or err, 76 failures and 12 errors: compat 46, map 6,
  se-conv 4+2, se-tier 10, mo-coding 2, se-check 1, cens-trunc 2+1,
  families-extreme 0+3, brms-utilities 0+1, loo-newdata-fitted 5+4,
  royston-parmar 0+1. 25 of 25 `lib:` lines rellib-r19.
- Scripts written for the later phases: rel081-chain2.sh (survey with
  frames and key tables; audit with floors; ledger with docs; checks;
  sessions), rel081-keytables.sh, rel081-floor.sh, rel081-audit.sh,
  rel081-survey.sh, rel081-commit.sh, rel081-rerun.R and the adapted
  check, docs, counts, installed, sessions-sum, ledger scripts.
- Plain tier (every file, NOT_CRAN only, 3 at a time; 07:50 to 08:14):
  510 of 510 with a RESULT line, 0 files with fail, err or warn
  (`dev/ciharden-log/r081-plain.sum`). Gated tier started 08:14.
- Gated tier (every gate; 08:14 to 09:27): 510 of 510, 36407 pass, 1 fail
  (the expected methods:326 signal), 0 err, 0 skip, 0 warn. Plain 32650
  pass, 289 skip. Counts `dev/rel081-log/counts.txt`, every `lib:` line
  rellib-r20. 0.80.0: 36118 gated, 32360 plain.
- Survey started 09:17 (frames, then 172 pairs, 1 at a time) beside the
  last gated file; key tables started beside it (1 process).
- Compat frame diff rellib-r19 -> rellib-r20, all eight extensions loaded
  (`dev/rel081-log/frame-diff-r19-r20.txt`): 12012 pairs, 7518 core, none
  only on one side; 61 status moves (43 untested -> refused, 18 works ->
  refused), 31 of them core; 9 note-only changes, all core. The lane's and
  its review's numbers exactly.
- Survey (RECORD=1, 1 pair at a time, 09:17 to 09:34; measured in 992 s;
  log `dev/rel081-log/survey.txt`): 172 of 172 `lib:` lines rellib-r20
  0.81.0, none FAILED. **Headline 89.7, Wald [85.2, 94.2]** (0.80.0:
  89.0 [84.5, 93.6]); no defect 93.0 [89.2, 96.8]; weak evidence out 89.0
  [84.4, 93.6], the one weak pair p069; agree among measured 96.4; L1+L2
  49 (57), L1 83 (99), L4 6 (9); compat-row proposals 0; "17 of 172 pairs
  differ" between stratum and live row. Translator's projection 89.7 /
  89.0 and compat's 89.0 with 0 proposals: both hold, combined.
- Verdict changes against the 0.80.0 record (`dev/rel081-equiv-delta.R`,
  `dev/rel081-log/equiv-delta.txt`): 22 pairs, 11 translator (p021, p022,
  p023, p042, p047, p048, p094, p109, p131, p136, p143) and 11 compat
  (p010, p011, p014, p059, p146, p154, p158, p166, p167, p169, p171), 0 of
  no lane. l1_const and l1_grad moved (relative 1e-6) on 7 pairs, each a
  lane pair; l2_metric and l4_z on 0.
- Key tables (`bash dev/rel081-keytables.sh`, logs `dev/rel081-log/kt/`),
  all on rellib-r20:
  - code8: at-max battery 99 fits / 32 code 8 and 220 / 77, 109 of 109
    give the message, 0 warned, gain 2.42e-11 to 9.92e-06; short stops
    "code-8 stops 99 ; short by more than 2e-5: 80 warned 80 message 0 ;
    short by 1e-6 to 2e-5: 13 message 7 warned 6 ; ratio 0.9995 1.0002";
    ridge 41 of 41 lose exactly a, k, theta_1, c0 SE 70.48 to 70.68,
    sigma 0.08333, no other condition; reviewer's seed set 52 of 52, 0
    gradient warnings; fixture 1 "h326 TRUE 21 of 21 ; grad warning 0 ;
    short warning 0 ; boundary message 21"; design-null stops (50, 51 tsv
    lines) byte-identical on rellib-r20 and rellib-r19 and to the
    reviewer's base tsv. Every number the lane's final build read.
  - support: review-support-is.R, is2.R, nq.R equal to the reviewer's base
    logs with timings removed (0 lines differ; nested seed 2 xi 0.5
    -242.9981897; 0 non-finite offsets); 33 fits (26 random-effect, 7
    fixed-effect) hex-identical to the lane's base logs, pinned to core 0;
    support-ndt-probe.R and support-gev349-probe.R identical to the lane's
    "after" logs.
  - hollow: review-hollow-probe6.R identical to the reviewer's log but
    the lib line: every changed read variable refused, the data and
    fit$data accepted for all 9 models, the cens() shortcut lines.
- 09:36 `bash dev/rel081-chain2.sh 2 3 4 5` started in the background
  (audit + floors; ledger + docs + regenerated brms-suite files; R CMD
  check core, sample, survival; plain sessions). Not yet run: OpenBLAS
  (`dev/rel081-ob.sh`), gated sessions, the commit dry run, the
  round-handoff update, the record's verification sections.
- The coordinator asked for the hand-back at this point; chain2 was
  still running at hand-back. Its results land in dev/rel081-log/
  (audit.txt, floor.txt, ledger-*.txt, docs-check.txt, suite.txt,
  check.txt, sessions-*.txt) and dev/rel081-log/chain.txt marks each
  phase.
- Phase 2 (09:36 to 09:53): translator audit on 133 pairs (the lane's 131
  plus p069, plus p154, now L1): 133 of 133 with a result, 67 pairs with
  latent blocks, 131 blocks, 2 untestable, BLIND 0, none under 100 x tol;
  133 lib lines rellib-r20. Floors: frmtmb.sample 0.28.0 over 0.80.0
  (floor lowered) loads and fails the draws readers' files (see
  dev/ciharden-log/r081-floor-frmtmb.sample.sum); 0.27.0 over 0.81.0
  loads; frmtmb.survival 0.2.3 over 0.80.0 runs 15 of 15 clean.
  frmtmb.sample's floor arm: 25 of 87 files fail or err (5 failures, 68
  errors), each `row_lpdf(..., exact = TRUE)`: "unused argument".

### Second consolidator (from 09:55)

- Took over with chain2 in phase 3 (docs build running). chain2 is not
  mine; I wait on its processes and do not touch its outputs.
- Phase 3, ledger: bin 1 recorded (14 RESULT lines), `brmsport-ledger.R`
  exit 0. **bin 2 did not run**: `dev/brmsport-bin2-all.sh` writes its
  first log to `dev/brmsport-bin2-log/survey.txt`, and that gitignored
  folder does not exist in a fresh worktree (only
  `brmsport-bin2-record.sh`, which runs later, makes it). Exit 1 before
  any R. gen and the guards ran (70 of 70). Rerun by itself as
  `dev/rel081-bin2.sh` (makes the folder, copies the bin-2 files the
  pipeline regenerates first), 1 R process, started 10:07 beside the
  last suite file. It writes only under `dev/`.
- Bin-1 ledger diff (`dev/rel081-ledger-diff.R`,
  `dev/rel081-log/ledger-diff.txt`): pass 397 -> 398, cannot transfer 48
  -> 47, divergence 49, defect 0, 494 rows. 2 rows move:
  `brmsfit-methods:326` cannot transfer / fixture -> pass (expected,
  code8) and `brmsfit-methods:394` divergence / hollow -> divergence /
  hollow, reason now "divergence by decision ... Savage-Dickey ratio
  ... no draws and flat priors" (expected, hollow). `:469`, `log_lik:580`,
  `log_lik:603`, `posterior_predict:440` are bin-2 rows: they wait on
  the bin-2 rerun.
- Docs (`dev/rel081-log/docs-check.txt`): roxygen again over all nine
  moved 0 files; installed equals the tree, 0 differ and 0 absent in
  all nine (frmtmb 2410, sample 304, latent 67, coupling 24, eam 194,
  learn 47, ode 38, spatial 44, survival 50); `check_pkgdown()` OK on
  all nine; build-docs 9 of 9 sites, DOCS ok, 10.9 min.
- loo -Inf characterization: `dev/rel081-loo-inf.R` (log
  `dev/rel081-log/loo-inf.txt`), the draws of `rel081-loo-exact.R`. 62
  cells of `log_lik()` are `-Inf`, 0 `+Inf`, 0 NaN. loo 2.10.1's
  `loo.matrix()` hands `-x` to `importance_sampling.matrix()`, whose
  `validate_ll()` refuses any `+Inf` ("All input values must be finite
  or -Inf"); the 62 cells are `+Inf` there. `validate_ll(log_lik)`
  itself passes. `loo()` and `psis()` stop; `waic()` does NOT stop
  (`waic.matrix()` validates `x`, which holds only `-Inf`): elpd_waic
  NaN, 28 pointwise rows non-finite. brms 2.23.0's
  `dgen_extreme_value()` at the same 3 cells: NaN NaN NaN ("NaNs
  produced", `log(t)` of a negative t); loo of that matrix stops with
  "NAs not allowed in input" (and waic too: `validate_ll()` tests
  `anyNA()` first). brms's Stan program maps `tmp_xi` by `scale_xi()`
  onto the region that keeps every row inside. So the backlog item and
  the sample NEWS bullet named `waic()` wrongly, and the backlog's
  "does not reach it while ..." was too broad: corrected by hand in
  `dev/test-backlog.md` and `extensions/frmtmb.sample/NEWS.md`, and in
  `dev/rel081-backlog.R` and `dev/rel081-news.R` with the same words
  (the NEWS edit after the docs build, before the checks).
- Regenerated brms-suite files (chain2 phase 3, `dev/ciharden-log/
  r081-suite.sum`, every gate): 15 of 15 clean, every `lib:` rellib-r20.
  `test-brms-suite-methods.R` pass 162 fail 0 (the focus and gated tiers
  read 161 and 1 failure at :419 before the regeneration). gen rewrote
  only `tests/testthat/test-brms-suite-methods.R` (38 lines).
- Bin 2, first rerun (`dev/rel081-log/bin2-2.txt`, 10:06 to 10:08):
  exit 0, 6 of 6 tier files clean, but `brmsfit-methods:469` stayed
  divergence (hollow). Cause: hollow edited
  `dev/brmsport-bin2-verdicts-manual.tsv` by hand (removed the :469 row,
  reworded :508 and :544), and `brmsport-bin2-all.sh` rewrites that
  file from `dev/brmsport-bin2-verdicts-2t.R`, which still had the 0.80.0
  rows (the regenerated file equalled HEAD's). The 2p files did not
  move: support and hollow put their 2p rows into the 2p generator.
  Fix (a consolidation call): `dev/rel081-bin2-moves.R`, hollow's three
  rows as its tree has them, sourced by `-2t.R` after the 0.73.0 moves,
  in the form of `dev/rel073-bin2-moves.R`.
- Bin 2, second rerun (`dev/rel081-log/bin2-3.txt`, 10:09 to 10:11):
  exit 0; 6 of 6 tier files clean (methods 58 and 58, distributions 92,
  log_lik 105, posterior_predict 65, priors 1), every `lib:` rellib-r20;
  `brmsport-bin2-verdicts-manual.tsv` now equals the merged tree's
  (hollow's) file. Bin-2 diff (`dev/rel081-log/bin2-diff.txt`): 308 rows,
  none only on one side; 3 rows move: `brmsfit-methods:469` divergence /
  hollow -> pass, `log_lik:603` and `posterior_predict:440` cannot
  transfer / not-ported -> cannot transfer / brms-internal; 1 value check
  moves: `log_lik:580` error / cannot transfer / reference -> holds /
  pass. Bins 1 and 2: 636 of 823 (77.3%), 636 of 724 without the
  decisions (87.8%); 0.80.0 634 of 823, 634 of 726. Every expected flip
  of the brief's merge step 7, and no other.
- `Rscript dev/rel081-news.R` after the NEWS edit (10:13): the three NEWS
  files reproduced byte for byte (md5 in `dev/rel081-log/news-md5-*.txt`).
- Early commit dry run (`dev/rel081-log/commit-dry-early.txt`): lanes 33,
  55, 86, 48, 51 paths (the merge's numbers), rel081 315; compat leaves
  `dev/compat-m-race1/` to `3/` and `dev/review-compat-out2/` untracked;
  fast-forward possible; index sha256 of all 7 trees unchanged.
- Record (`dev/round-20261109.md`): merge steps, review texts, conflicts,
  decisions, key tables, tiers, survey, ledger, versions, not done, for
  the user; checks, sessions, OpenBLAS and the final dry run wait.
  `dev/round-handoff.md`: a 0.81.0 section.
- R CMD check of frmtmb started 10:10 (chain2 phase 4).
- 10:24 `bash dev/rel081-after.sh 30964` started in the background: it
  waits for chain2's outer shell to exit, then runs `dev/rel081-ob.sh`
  and the gated sessions (`rel0.81.0g`, 3 at a time); its lines go to
  `dev/rel081-log/chain.txt` as AFTER.
- Phase 4 (10:06 to 10:49), `dev/rel081-log/check.txt`: frmtmb 0.81.0,
  frmtmb.sample 0.28.0, frmtmb.survival 0.2.3 each `Status: 1 NOTE` (the
  HTML manual's math rendering, V8 not installed); tests 23 min, 377 s,
  15 s. Check logs copied to `dev/rel081-log/checklogs/`; the check tree
  (41 MB) and the built site `dev/docsci-site-rel081` (47 MB) deleted.
- Phase 5 (plain sessions) started 10:49.
- Phase 5 (10:49 to 11:43), `dev/rel081-log/sessions-rel0.81.0-sum.txt`:
  27 of 27, 0 SESSION lines with fail, err or warn; 243 of 243 `lib:`
  lines rellib-r20; LEAK 46, 0 untagged. frmtmb 23487, 23652, 23755
  (forward, reverse, preload), sample 4020, survival 556 each. chain2
  END 11:43:01.
- OpenBLAS (`dev/rel081-after.sh` -> `dev/rel081-ob.sh`, 11:43 to 11:45,
  `dev/rel081-log/ob.txt`): both copies built (exit 0; Rblas 74 / Rlapack
  616 symbols to 0.3.26, 76 / 617 to 0.3.32). 0.3.26: se-tier 112,
  se-check 134, families-extreme 124, cens-trunc 99; 0.3.32: se-tier
  112; all fail 0 err 0 warn 0, `blas:` lines name the copies, `lib:`
  lines rellib-r20. `dev/ciharden-out/` deleted by the script.
- Gated sessions (`rel0.81.0g`, 3 at a time) started 11:45.
- Gated sessions (11:45 to 14:35, `dev/rel081-log/sessions-rel0.81.0g-sum.txt`):
  27 of 27, 0 SESSION lines with fail, err or warn; 243 of 243 `lib:`
  lines rellib-r20; LEAK 194 (142 a Stan model load), 0 untagged. frmtmb
  26672, 26837, 26940; sample 4387; survival 556. AFTER END 14:35:35.
- Final commit dry run (`dev/rel081-log/commit-dry-run.txt`) after the
  record's last edit; record and hand-off complete.

### Release review fixes (from 14:50)

- `dev/reviews/2026-11-10-release.md`: SHIP, five minors. The
  coordinator asked for all five; finding 2 by its call (all three stop).
- `dev/rel081-release-fix3.R` (`fix3-1.txt`, `fix3-2.txt`): the loo.R
  edits (made first by hand, then held in the script, which reads them
  "already applied"), the test appended from
  `dev/rel081-fix3-test-block.R`, the backlog item, the commit script (7
  edits), the brms entries of `dev/upstream-bugs.md` sorted (43 50 44 45
  49 46 47 48 51 -> 43 to 51), then `dev/rel081-news.R` (BREAKING on
  code8 2 to 4 and compat, rewrapped; sample NEWS: one floor line, the
  decided and measured limit bullet). The first run stopped at the
  `system2()` call (a quoted short path); the second ran clean. `bash -n`
  on the commit script clean.
- `dev/rel081-install-sample.sh` (15:01:14): roxygen wrote
  `sample-loo.Rd`, sample 0.28.0 reinstalled into rellib-r20.
- New test file on r20: pass 80 (53 before), 6 s. On r19: pass 25, fail
  7, err 6; the two new tests fail 2, err 2. `dev/rel081-exactoff.R` (r20,
  reader set to `exact = FALSE`): 8 tests, failed 2, errors 1, all in
  the two new tests.
- `dev/rel081-fix3-verify.sh` (15:02 to 15:51): sample roxygen moved 0;
  installed equals tree in all nine (sample 305); check_pkgdown needed
  pandoc on PATH (first read "Pandoc not available"; rerun OK on all
  nine; script fixed); sample plain 87/87 pass 4047 skip 18, gated 87/87
  pass 4414, 0 fail/err/warn, every lib line rellib-r20; sessions 3/3
  plain (4047 each order) and 3/3 gated (4414), LEAK 14 and 38, 0
  untagged; R CMD check sample `Status: 1 NOTE` (V8), tests 438 s.
- Commit dry run (`commit-dry-run-fix3.txt`): exit 0; lanes 33 55 86 48
  51; rel081 335; index sha256 of all 7 trees unchanged, now read before
  the first git call. `git ls-files --eol -m`: 93 `i/lf w/crlf`.

### Short check of the release fixes (from 15:55)

- The review's section "Short check of the release fixes": FIX on one
  minor. `pp_check()`'s `loo_*` types read `log_lik()` in
  `draws_ppc_psis()` and skipped `loo_check_finite()`. Also: drop the
  compat bullet's BREAKING mark.
- `dev/rel081-release-fix4.R` (`fix4-1.txt` to `fix4-4.txt`):
  `loo_check_finite(ll, "pp_check(type = '<type>')")` in
  `draws_ppc_psis()`; the second new test expects the class from the
  three `loo_*` types; `rel081-news.R` no longer marks the compat
  bullet (7 BREAKING in `# frmtmb 0.81.0`); the commit script names
  fix4. Its first run skipped the deletion edit as "already applied"
  (new text inside old); the edit check now tells insertion from
  deletion, the third run applied it, the fourth read 6 of 6 applied.
- Sample reinstalled 16:03:01; roxygen wrote no Rd.
- Reviewer's probe (`dev/review-rel081-fix3-paths.R`, log
  `dev/rel081-log/fix4-paths.txt`): every criterion path, the three
  `pp_check()` types included, raises `frmtmb_loglik_not_finite` ("62
  cells, at 28 of 40 observations and 3 of 200 draws");
  `loo_moment_match()`, `kfold()` refused as before; controls answer.
- `dev/rel081-fix4-verify.sh` (16:03 to 16:50): roxygen moved 0;
  installed equals tree in all nine (sample 305); check_pkgdown OK on
  nine; sample plain 87/87 pass 4050 skip 18, gated 87/87 pass 4417, 0
  fail/err/warn, lib lines rellib-r20; `test-loo-newdata-fitted.R` 83;
  sessions 3/3 plain (4050 each order), 3/3 gated (4417), LEAK 14 and
  38, 0 untagged; R CMD check sample `Status: 1 NOTE` (V8), tests 379 s.
- Commit dry run (`commit-dry-run-fix4.txt`): exit 0; lanes 33 55 86 48
  51; rel081 338; index sha256 of all 7 trees unchanged; 93 `i/lf
  w/crlf`. No R process left.
