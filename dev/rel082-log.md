# Running log of the 0.82.0 consolidation (round of 2026-11-11)

Release tree `C:/Users/adf44/source/r/frmtmb-wt-rel082` (branch `wt-rel082`, off
`e8abed0a`), library `C:/Users/adf44/source/r/rellib-r21`, base `rellib-r20` (0.81.0).
Machine output in `dev/rel082-log/` (gitignored).

## Lane records read (2026-11-11)

The last section of each review, read against the brief's "all twelve: SHIP":

- SHIP: quadseam, boundfam (short check), wcheck (short check 3), nlstop, silent, priortab,
  names, arspat, sampfront, ordstart (short check).
- ppdraw: the review's last section is "Short-check verdict: FIX" (one line: guard
  `sqrt(diag(H))` in `draw_edge_csd()`). The lane applied it after the review
  (R/predict-brms.R 18:13; `1 / sqrt(pmax(diag(H), 0))`); no reviewer has read it.
- gpstop: the review's last section is "Verdict: FIX, small" (S1, the span rule's lower
  end). The lane shipped `gp_relin_ls_min = 0.5` of the median gap, not the review's "a
  tenth" (a tenth does not refuse the plain point), and `gp_inner_pd()`; no reviewer has
  read it.

Edits after the last review section, unreviewed:
- boundfam: `se_null_drop_held()` in `sdr_rescue_hessian()` (R/se-check.R), its test in
  `test-se-limit.R`, a NEWS clause (the "held limit directions" fix that ppdraw's short
  check asked of boundfam).
- nlstop: S2 fix, `free_dispersion(frame)` (R/fit-restart.R), a zip test.
- names: the whole fix round (minors 1 to 8), including `pars =` order in
  frmtmb.sample `draws_accessor_args()` / `draws_select_variables()` and
  `block_display_label()` in R/multiple.R. The reviewer's names + sampfront merged run
  predates it.
- ordstart: edits at 16:39 after the 16:35 short check, before its last t5 runs and check.
- quadseam: a test comment only. arspat, sampfront, priortab: backlog filing only.
  wcheck, silent: none.

Paths the merge must not take (untracked, not ignored, machine output):
- `dev/reviews/nlstop-review-session/` (26 files: copies of test helpers and a stale
  `test-fit-restart.R`).
- `dev/reviews/names-review-core.diff`, `-ext.diff`, `-news.diff`, `-cost.prof`.
- Job lists (`dev/silent-*.jobs`, `dev/gpstop-*.jobs`, `dev/reviews/*.jobs`) are taken,
  as 117 tracked `dev/*.jobs` are.

Text fixes for the release:
- wcheck NEWS bullet 4 word order (the brief).
- gpstop's backlog diff deletes the blank line before "### MINOR and performance" (base
  line 5035): restore it.
- sampfront's CLOSED line under 4717 says "subset taken by row names"; its fix round also
  places rows by value.
- frmtmb.sample NEWS: one "Needs frmtmb 0.82.0" line for names, priortab and sampfront.

## Done (2026-11-11 to 2026-11-12)

- Merge (`dev/rel082-merge.sh`, `dev/rel082-log/merge/`): 0 conflicts; boundfam 50 paths,
  quadseam 59, ppdraw 49, wcheck 33, priortab 40, nlstop 58, gpstop 48, ordstart 44,
  silent 73, arspat 66, names 78, sampfront 62; 30 skipped (nlstop 26, names 4); 140
  modified files written CRLF (`git ls-files --eol -m`: 140 `i/lf w/crlf`).
- NEWS (`dev/rel082-news.R`): core 57 blocks in 13 sections, 31 BREAKING; sample 22
  bullets and one Needs line; spatial 0.4.0, survival 0.3.0, latent 0.7.4.
- Backlog (`dev/rel082-backlog.R`): heading renamed, gpstop's blank line restored,
  sampfront's 4717 line completed; 91 CLOSED lines of the round kept. Its first run
  dropped every blank line (strsplit of ""); `dev/rel082-remerge-one.sh` rebuilt the
  merged backlog and the fixed script reran.
- Merge steps 2 to 4 in R/fit.R and R/autoscale.R (one `cold_start`; chosen fit's
  messages; restart message after the mixture starts); `tests/testthat/test-cold-start.R`.
- Versions (`dev/rel082-release-fix.sh`); install into rellib-r21 (`dev/rel082-install.sh`,
  `dev/rel082-log/install-1.txt`).
- Tiers (`dev/rel082-chain.sh`, 2 processes, beside the key tables at 1): focus 70 of
  70, focus on rellib-r20, plain 527, gated 527. Gated: 2 errors in
  test-bcm-psychophysics.R:162 and test-bcm-retention.R:234 (helper-stan.R's
  `frm_u_term()` read `ranef(fit)[[i]]` by block index; names' layout). Helper fixed to
  index by the block label; the 4 helper-stan files rerun gated, 0 fail
  (`r082-hstan.sum`).
- test-predict-draw-edge.R: ordstart's step-start warning allowed in ppdraw's ordinal
  mixture test; the 2 warnings at :101 (fitted() NaN SEs on `zi ~ x` and `hu ~ x`) are
  boundfam's on its own lane (`dev/rel082-fitted-held.R`), left as found.
- Key tables: subagent, `dev/rel082-keytables.sh`, `dev/rel082-log/kt/`.

## Running at the time of writing

- `dev/rel082-survey.sh all` (frames, then RECORD=1 survey, 2 processes).
- `REC_JOBS=1 dev/rel082-ledger.sh bin1` then `ledger` (1 process).

## Next

bin2 (`dev/rel082-bin2.sh`), ledger diffs, audit, check (core, sample, spatial, survival,
latent), sessions plain and gated, OpenBLAS (`dev/rel082-verify.sh ob`), docs, installed
code, floors, records, commit dry run.

## Later (2026-11-12)

- Key tables done (subagent; `dev/rel082-keytables.sh`); moved: boundfam suite cmp
  (other lanes), gpstop sweep (boundfam; attribution `gpstop-attrib`).
- Survey RECORD=1 done: headline 89.7 [85.2, 94.2], 0 verdict changes, l1_grad moved on
  p131 only; 172 of 172 lib lines rellib-r21. Frame diff: 5 moves, all arspat's.
- Ledger: bin 1 0 rows move; bin 2 0 rows move; names' two reasons regenerated.
- Running: `dev/rel082-audit.sh` (2), `dev/rel082-floor.sh 1` (1), then
  `dev/rel082-after.sh` (docs check + R CMD check of five packages, then chain2 phases
  3 to 5: OpenBLAS, plain sessions, gated sessions). Progress lines in
  `dev/rel082-log/chain.txt`.
- Record draft: `dev/round-20261111.md` (PENDING markers left: audit, floor, verification,
  what is not done, for the user, commit script, state).

## Finished (2026-11-12 06:10)

Audit, floors, docs check, R CMD check of five packages, OpenBLAS, plain and gated
sessions all done; record `dev/round-20261111.md` complete; `dev/round-handoff.md` has
the 0.82.0 section; commit dry run clean (`dev/rel082-log/commit-dry-run.txt`). No
process running. The "Record draft" line above is superseded.

## After the release review (2026-11-12)

Review `dev/reviews/2026-11-12-rel082.md` (FIX). Majors 1 and 2 and minors 1 to 7
applied; lane boundfam's consolidation fixes re-merged (`dev/rel082-remerge-boundfam.sh`);
reruns by `dev/rel082-chain3.sh` (END 12:02:12). Open: `test-se-limit.R:200` fails under
OpenBLAS 0.3.26 (boundfam's new test; a rounding-decided branch of the SE check).
Record: `dev/round-20261111.md`, "After the release review". No process running.

## After boundfam's OpenBLAS fix (2026-11-12 18:30)

Re-merged (`dev/rel082-remerge-boundfam.sh`; boundfam 61 paths, MERGED updated);
`dev/rel082-chain4.sh` (13:25 to 18:25): test-se-limit 240 under reference, OpenBLAS
0.3.26 and 0.3.32; plain 33852 and gated 37615, fail 0, warn 0, only test-se-limit's
count moved; core sessions clean; core check 1 NOTE; dry run 6 clean. No process running.

## After the short check (2026-11-12 19:00)

S1: `take_list()` reads content by `git hash-object --stdin-paths` against
`git ls-files -s`, never `git diff`; proof in `dev/rel082-log/proof/` (touch in 13 trees,
two dry runs, 14 index sha256 equal outside and inside). Multi-lane drift check added
(54 pairs clean, 31 exempt). Two LOW items filed. Dry run 2's index move was the script.
