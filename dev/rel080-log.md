# Consolidation log, round of 2026-11-06 (frmtmb 0.80.0)

Running log of the consolidator. Brief: `round-briefs/rel080-brief.md`.
Release tree `frmtmb-wt-rel080` (branch `wt-rel080`, off 3cd2a0c8), library
`rellib-r19`. Machine output goes under `dev/rel080-log/` (gitignored).

## State

- Started 2026-11-06. Release tree clean at 3cd2a0c8; `rellib-r19` did not exist.

## Reading

- Reviews' last sections read for all twelve lanes. Verdicts that stand:
  ci SHIP, latperf SHIP, mfxse SHIP (second short check), loglik SHIP,
  osa SHIP, nestquad SHIP, arstep SHIP, fitedge SHIP, nlstart SHIP (third
  short check), readers SHIP (S1 confirmation), superassign SHIP.
- bigx: the review's LAST section ("Confirmation of the M2 fix") reads
  **FIX** (run `design_rank_qr()` on the centered design). The lane's
  findings section 18 says it did this after; no review confirms it. To check
  on the merged tree (test-design-offset.R, `dev/review-bigx-fx2.R`).
- fitedge: the lane edited `R/fit.R`, `R/gp-start.R` and tests after the
  SHIP short check (findings "After the short check": m1 `gp_body`, m2
  wording, m3 classed stop warnings). Not reviewed. The merge2 inputs of
  `rev-fitedge-src/merge2/resolve.R` therefore differ from the lane's current
  `R/gp-start.R` (one roxygen paragraph) and `R/fit.R`.

## Merge

(see below)

### Merge, done 2026-11-06

`bash dev/rel080-merge.sh` (log `dev/rel080-log/merge/`: `taken.txt`,
`conflicts.txt`, `resolved.txt`, `skipped.txt`). Order as the brief:
ci, latperf, mfxse, loglik, osa, nestquad, arstep, fitedge, bigx, nlstart,
readers, superassign. Each lane's whole tree (modified and untracked, not
ignored), except machine output dev/.gitignore lets through:
`dev/review-*-out*/`, `dev/*-equiv-out*/` (bigx's 172 survey tsv files),
`dev/review-*-merge/` (superassign's trial merge). NEWS files, the backlog
and upstream-bugs are merged with `--union`, then rebuilt by
`dev/rel080-news.R` and `dev/rel080-backlog.R`.

A first try (no union, no per-lane resolve) nested markers in R/fit.R and
NEWS; I put the tree back to 3cd2a0c8 by writing each taken path from
`git show` (no git write) and reran. Conflicts and resolutions
(`dev/rel080-resolve.R`, run after each lane):

| lane | file | blocks | resolution |
|---|---|---|---|
| mfxse | R/autocor.R | 1 | latperf's two paragraphs, mfxse's sentence on the wider block a carried `vcov =` inverts |
| mfxse | man/frmtmb-autocor.Rd | 1 | release side; roxygen writes it again |
| arstep | R/fit.R | 1 | brief step 2: arstep's `gp_second_descent()` lines, then nestquad's branch |
| fitedge | R/fit.R | 1 | the same two lines, already in: release side |
| fitedge | R/gp-start.R | 5 | `diff3 -m -E fitedge base arstep`, then `rev-fitedge-src/merge2/resolve.R` verbatim (in rel080-resolve.R). Result differs from the reviewed `resolved_gp-start.R` only in fitedge's later roxygen paragraph |
| bigx | R/fit.R | 1 | `gp_second_descent()`, then bigx's `nl_basis_descent()`, both in the Laplace branch, then nestquad's branch |
| nlstart | R/par-template.R | 1 | loglik's `llw` word, nlstart's `class = "frmtmb_nl_restart"` |
| superassign | R/me.R | 2 | readers' side (superassign review MINOR 3); `me_node_value()` equals readers' |

Checks: no `<<-`/`->>` in R/ or extensions' R/; 0 markers left.
Line-end fix: the union merge of the backlog left a lone CR at EOF (git
read the file as -text); merge script's `crlf()` now drops it.

### Release edits (`dev/rel080-release-fix.R`)

Step 3 (probe on `logLik_penalized()`, `llw` in the multistart message and,
beyond the brief, in the unsettled warning, which compares the same
penalized values), step 6 (basis attached to the final nested tape),
step 8 (`interop_re_part()` adds every non-field block's theta under
`allow_new_levels`), step 13, versions, floors, roxygen 8.1.1.
NEWS: `dev/rel080-news.R`. Backlog: `dev/rel080-backlog.R`.

## Build

`bash dev/rel080-install.sh C:/Users/adf44/source/r/rellib-r19` (from
rel079's, paths changed; log `dev/rel080-log/install1.txt`): roxygen on core
and the eight extensions, all nine installed 21:36 to 21:37, drmTMB 0.7.0
copied from rellib-r18. Scripts copied from rel079 with paths changed:
install, rox, env, verify, check, counts, sessions-sum, installed,
docs-check, floor, survey, equiv-delta, equiv-cmp, ledger, ledger-diff.

## Merge-step checks on the first build (rellib-r19, 21:37)

- Step 6, `dev/rel080-nestshift.R` (bernoulli, 15 schools x 3 classes x 4,
  sds 0.8, slope 0.7; `dev/rel080-log/nestshift.txt`, against
  `wt-nestquad-lib` in `nestshift-nestquadlib.txt`): at `x + 1e4` and
  `x + 1e6` the merged fit is 0.0007 to 0.039 below the centered fit, slope
  0.7182 against 0.7302 (seed 1, 0.06 SE), with `frmtmb_quad_unsettled`,
  "false convergence (8)" and SE NaN on 3 of 4. nestquad's own build: 1.2
  to 12.6 short, slope 0.0000008 to 1.003. The rounds in natural
  coordinates are what is left: loud, not silent. Added to NEWS Limits.
- Step 8, `dev/rel080-newlevel.R` (poisson `(1 | g)`, seed 7411, 3 rows at
  a new level; `newlevel.txt`): `avg_predictions()` SE 0.457649 equals the
  delta method over the coefficients and the block's log sd (ratio
  1.00000); over the coefficients alone it would be 0.406931. Estimate
  equals `fitted()` and the closed form, 2.273276. 0.18 s.
- Step 10, `dev/rel080-latnl.R` (negbinomial, `ar(cov = TRUE)` on `mu`,
  `nlf(shape ~ exp(a + b * z))`): at 150 and 400 rows, `frm_linpred(se.fit)`
  on `mu` and `shape`, in sample and on new rows, equals `frm_lp_basis()`'s
  `A V A'` to 4.4e-16; `fitted()` and `conditional_effects()` finite.

## The user's decisions, received through the coordinator (2026-11-06, evening)

1. OpenBLAS: approved. Run `dev/ciharden-openblas.sh` 0.3.26 and 0.3.32
   (LAPACK) on the merged build, on the files the lanes used, including
   readers' `test-predict-draw-mean.R`, `test-me-newdata.R` and the ce test
   files. If the permission system denies it, report the denial and stop.
2. TMB internals: keep `solveSubset()` via `get0()` and `info()` via
   `getFromNamespace()`. Filed in `dev/test-backlog.md` (0.80.0 section,
   by `dev/rel080-release-fix2.R`): ask TMB upstream to export
   `solveSubset()` and a tape-size accessor (`info()`); user will raise it
   later. Nothing posted.
3. `dev/release/checklist.md` written (how-to; first item the cran-lag
   comparison; steps taken from the 0.79.0 and 0.80.0 briefs and
   `dev/round-handoff.md`).
4. Every other decision of the brief stands.

## Focus tier, first build (`r080-focus`, 51 files, every gate, 2 at a time)

`dev/ciharden-log/r080-focus.sum`: ran 51 of 51, every `lib:` line rellib-r19.
One failure: `test-quadrature-nested.R:430`, "the re-anchoring loop does
not cycle". Semantic conflict fitedge x nestquad: fitedge's
`check_convergence()` now gives its stop warnings their own classes
(`frmtmb_conv_code`, ...) ahead of `frmtmb_warning`, and the test read
`class(w)[1]`. Fixed in the test by inheritance (`dev/rel080-release-fix2.R`).

## Running (started 2026-11-06 21:50)

- Chain `dev/rel080-chain.sh` (bash PID 897897, Git Bash): stages base
  (focus files on rellib-r18), ref (gated, all suites), plain, 2 R processes
  each. Progress: `dev/rel080-log/chain.txt`; results
  `dev/ciharden-log/r080-focus-r18.sum`, `r080-ref.sum`, `r080-plain.sum`.
- Key tables batch b1 `dev/rel080-keytables.sh b1`, 1 R process, output
  `dev/rel080-log/rr/*.txt`, progress `dev/rel080-log/rr/progress.txt`.
- Superassign: `test-no-superassign.R` pass=16 fail=0 on rellib-r19;
  `dev/review-superassign-globals.R` on the merged tree (through
  `dev/rel080-rerun.R`): 90 changed functions, 8 with new free variables,
  every one a package-level definition another lane added
  (`design_alias_tol`, `profile_run`, `prior_df_why`, `me_gh_*`,
  `me_net_tol`, `jcs_block_max`) or a base function: no reference to a
  variable the rewrite moved.

## Key tables, batch b1 (rellib-r19, 21:51 to 21:54; `dev/rel080-log/rr/`)

Run through `dev/rel080-rerun.R` (each lane library and lane-tree path in a
script and in what it sources reads rellib-r19 and this tree).

| lane | script | expected | read |
|---|---|---|---|
| loglik | `loglik-routes.R after` | max rel 3.15e-11; E REML b prior refused | 3.15e-11; refused, rel 1.92e-03 |
| latperf | `latperf-agree.R` | 10 cases sparse TRUE, max 9.99e-15 | same, 9.99e-15 (unstr) |
| latperf | `latperf-prof.R lane 8000 one` | fitted 0.42 s, se_first 0.528786047, saved 12.27 MB | 0.30 s, 0.528786047, 12.28 MB |
| latperf | `latperf-prof.R lane 8000 two` | fitted 0.50 s, saved 15.67 MB | 0.36 s, 15.67 MB |
| latperf | `review-latperf-pattern2.R 8000 zipzi / hurdle 500` | 0.63 s / 1.00 s, 0 rows by row | 0.41 s / 0.66 s, 0 rows, 50 / 500 columns |
| bigx | `review-bigx-fx2.R` | code-8 29, gains 1.95e-09 to 9.92e-06, 0 short, 0 other warning | same; 1e6 products and powers capped (`frmtmb_design_capped`) |
| nlstart | `review-nlstart-fix3-inv.R` (brief step 7) | 0 short at 1, 1e-4, 1e-6; 5 at 1e-8, 18.19 to 37.30, all warned | same (`nlstart-fix3-inv-sum.txt`) |
| readers | `review-readers-sc2.R` | 2.41468 / 2.42274 / 222685, nlf 2.81272, NaN, all warn | same |
| readers | `review-readers-sc3.R` | 10 pi row 0.509284 silent | **0.882624 with `frmtmb_mean_unsettled`**: the lane's later fix (irrational steps, findings "Confirmation minor", not reviewed); the readers lane library reads the same |
| readers | `review-readers-inf.R` | case 1 2.69294 and 9247987 warned | same |

Lanes edited after their review's last section (from their findings;
none re-reviewed): fitedge (m1 to m3), bigx (section 18, the centered rank
test that the review's last FIX asked for), readers (irrational steps),
nestquad ("Fix round 2": rule error in the warning, four memory entries,
`profile()` refusal order). The backlog items are written to match.

## New tests on rellib-r18 (`r080-focus-r18`, 21:50 to 21:59)

The 51 new or changed test files of the release tree on rellib-r18 (0.79.1 and its
extensions), every gate on, 2 at a time: 51 of 51 ran, every `lib:` line rellib-r18.
220 failures in 18 files and 77 errors in 28 files; 34 of 51 files fail or err
(`dev/ciharden-log/r080-focus-r18.sum`). The 17 that pass on 0.79.1 are tests the
lanes edited for new classes or rounding (the census files, ci's tests) or pins that
hold on both builds.

## Key tables, batch b2 (rellib-r19, 21:56 to 22:31) and summaries

Lane output folders copied into this tree (gitignored) for the before arms:
`dev/review-fitedge-out`, `review-arstep-out`, `arstep-out`, `nlstart-out`,
`review-nlstart-out`, `mfxse-rate-out`, `bigx-out`. The r19 runs write into
those copies.

| lane | script | expected | read |
|---|---|---|---|
| osa | `review-osa-time.R` | 2000 rows (1 \| g): cdf 41.9 s, old 225.0 s, 0.19; 500 rows 0.30 | 27.9 s, 167.5 s, 0.17; 500 rows 4.0 s, 12.1 s, 0.33 |
| arstep | `review-arstep-cost.R p600,p1000,a300 2` | ll 32.589404, 173.765935, 227.607924; p600 2.9 to 4.6 times, a300 1.29 to 1.34 | same ll, conv 0, w 0; p600 1.8 to 3.6 times (0.80, 0.36 s against 0.22, 0.20), p1000 1.3 to 1.6, a300 1.16 to 1.17 |
| arstep | `review-arstep-fresh.R` G, N, P + `arstep-fresh-cmp.R` | FCSUM 168 rows, base short 28, fixed short 11, worse 0, better 17 | same |
| arstep | `arstep-ng.R` 150 rows + `arstep-summary.R` | lane short 5 silent 4 worst -11; worse 0, better 8 | same |
| nestquad | `nestquad-nodes.R bern asym` | k 11 1.301e-04 ... k 31 6.468e-10 | same, every digit |
| bigx | `bigx-grid.R` + `bigx-summary.R` | 240 right, 11 warned, all at 1e8 | 240 right, 11 warned, 0 silent wrong; the code-8 warning now reads class `frmtmb_conv_code` (fitedge's class) |
| fitedge | `review-fitedge-gpfresh.R fresh` + gpsum | short 2 (silent 2), worst -36.500, worse 0, better 24 | same |
| fitedge | `... lane96` + gpsum | short 7, worst -0.635, worse 0, better 25 | short 6, worst -0.635, worse 0, better 26 (one row better than the lane's build) |
| mfxse | `mfxse-rate.R` x 4 + `mfxse-ratesum.R` | off: avg_pred 0.00, cmp_w 0.02, pred 0.01; levels 0.999998 to 1.000004; no warning; codes 149/19 | same |

Batch b3 (nlstart 450-row scan and the reviewer's fresh rows) started 22:32.

## Batch b3 (nlstart) and a semantic conflict bigx x nlstart (22:32 to 23:30)

- 450-row scan (`nlstart-scan.R` x 9 + `nlstart-scansum.R`): 0 short (scan-before
  16), higher than 0.79.1 on 16, lower 0, within 1e-4 on 434. As the lane.
- The reviewer's fresh rows (`review-nlstart-fresh.R`, `-fresh2.R` + `-freshsum.R`):
  after 14 short (5 silent) against the lane's 11 (4); probe-off differs from 0.79.1
  on 70 rows (lane 0; expected, since bigx's and fitedge's changes reach the first
  fit); 0 lower than 0.79.1. Row by row against the lane's own logs, 3 rows read
  lower: gaussian gomp 22 (-43.73), poisson two_sat_rev 33 (-2.86), bernoulli mm 11
  (-1.58), each where the lane's probe won and the merged build's did not.
- Cause (`dev/rel080-nlrows.R`, `-nlrows2.R`): on each of the three the first fit's
  tape carries bigx's `nl_basis_descent()` basis (the body's Jacobian whitened at the
  first stop), and the probe's pattern runs stepped in it from other signs. With
  `options(frmtmb.opt_basis = FALSE)` the probe wins and the fit equals the lane's.
- Fix (`dev/rel080-release-fix2.R`): `nl_sign_probe()` runs its patterns with the
  design basis alone (`opt_basis(frame, obj, ...)`, what `fit_assembled()` starts
  from) and puts the fit's basis back after. On a scratch install
  (`C:/Users/adf44/source/r/rel080-scratch-lib`, core only) the three rows equal the
  lane's: -206.714790, -462.035296, -205.167254 (`rr/nlrows-fix.txt`).
- Test: `test-nl-start.R`, "the sign probe runs without a basis whitened at the
  first stop" (the review's gomp row, seed 22004; skips if the first fit keeps no
  body basis). rellib-r19 before the fix: 2 failures (`rr/t-nlstart-r19pre.txt`);
  scratch with the fix: pass 55 (`rr/t-nlstart-fix.txt`).
- So the gated tier now running (`r080-ref`) is on the build before this fix and
  before the frm() roxygen sentence. Plan: stop the chain when `ref` ends, reinstall
  rellib-r19 (final build), rerun both tiers.

## Gated tier on the first build, chain stopped, final build (23:31 to 23:33)

- `r080-ref` (first build, 21:37): 509 of 509 files, every one fail=0 err=0 skip=0
  warn=0 (kept as `dev/rel080-log/r080-ref-prefix.sum`).
- Chain stopped at 23:31 in its first plain files (processes 897895, 897897 and
  their children, all mine; the R processes it ran were test files, which write
  nothing to a library).
- Final build: `dev/rel080-rox.R frmtmb` (wrote `frm.Rd` alone), core reinstalled
  into rellib-r19 at 23:32:17 (`dev/rel080-log/install2.txt`). The eight extensions
  stay at 21:37 (no source change since).
- Running: `dev/rel080-chain2.sh` (ref then plain, 2 processes; `chain.txt`), and
  key tables batch b4 (nl tables again, 1 process; `dev/rel080-log/rr-final/`).
- Queued (task waits for `CHAIN2 END` and `BATCH b4 END`, then runs
  `dev/rel080-chain3.sh`): obbuild (0.3.32; 0.3.26 built at 21:5x), ob (60 files:
  the 51 changed plus the 8 `test-ce-*.R` and `test-conditional-smooths.R`, under
  0.3.26 and 0.3.32 LAPACK, 4 threads), check (core, sample, latent, ode,
  survival), sessions (27 x 2 tiers, 3 at a time), docs, floor, survey (RECORD=1,
  3 jobs), ledger. `dev/rel080-record.sh` now runs 3 files at a time (rel079's ran
  every file at once).
- OpenBLAS 0.3.26 build at 21:5x briefly ran a 4th R process (its BLAS check)
  beside the chain's 2 and batch b1's 1. Noted; no other overrun.

## Key tables on the final build (batch b4, 23:33 to 00:2x; `dev/rel080-log/rr-final/`)

| table | first build (rr/) | final build | lane / review |
|---|---|---|---|
| `review-nlstart-fix3-inv.R` | 0 / 0 / 0 / 5 short (all warned) | same | same |
| nlstart 450-row scan | 0 short; 16 higher, 0 lower | same | same |
| reviewer's fresh rows (425, 335 fitted) | 14 short (5 silent); 3 rows below the lane | **11 short (5 silent), short > 1: 6 (3); 0 rows below the lane, 3 above**; 0 below 0.79.1; multistart on 87 | 11 (4 silent), short > 1: 5 (2); multistart 89 |
| fitedge 192 rows | fresh 2 short, lane96 6 | same | 2, 7 |
| readers sc2, sc3, inf | as review except the 10 pi row (lane's later fix) | same | |

One row moves against the lane on the fresh set: binomial mm seed 22 reads -643.4937
with no warning, 0.008 above the lane's -643.5014, which warned once; both 4.61
short of the best start (the review's third check lists this row as silent at the
1e-6 scale). Probe-off differs from 0.79.1 on 70 rows, max 45.4: the first fit now
carries bigx's basis and fitedge's ridge pass; nothing reads lower than 0.79.1.
(The two fitedge summaries ran as a 4th R process for a few seconds beside chain2's
2 and b4's 1.)
- ci `dev/release/cran-lag.R "<rellib-r19>;<user lib>"` (case D on this stack,
  `cranlag-r19.txt`): 285 packages, 280 current, 5 LAGS (charport, gert, igraph,
  nlme, RTMBdist 1.2.0 against 1.2.1), 0 missing, exit 1; the review's case D. Cases
  A to C and E to G (dead mirrors, proxies, scratch repositories) were not rerun.

## Gated tier, final build (`r080-ref`, 23:32 to 01:0x)

### gated tier (every gate on)

From `dev/ciharden-log/r080-ref-suite/` (509 jobs; 509 files with a RESULT line; 509 whose `lib:` line names rellib-r19 for the package and for frmtmb).

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 317 | 26430 | 0 | 0 | 0 | 0 |
| frmtmb.coupling | 11 | 557 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 33 | 1872 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 12 | 693 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 16 | 559 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 86 | 4334 | 0 | 0 | 0 | 0 |
| frmtmb.spatial | 8 | 561 | 0 | 0 | 0 | 0 |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |
| **total** | **509** | **36109** | **0** | **0** | **0** | **0** |

Files with a skip: 0.

Files with a failure, an error, an escaped warning or no RESULT line: 0 



## Plain tier, final build (`r080-plain`, 01:01 to 01:34)

### plain tier (NOT_CRAN only)

From `dev/ciharden-log/r080-plain-suite/` (509 jobs; 509 files with a RESULT line; 509 whose `lib:` line names rellib-r19 for the package and for frmtmb).

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 317 | 23245 | 0 | 0 | 239 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1848 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 12 | 687 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 86 | 3967 | 0 | 0 | 18 | 0 |
| frmtmb.spatial | 8 | 478 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |
| **total** | **509** | **32351** | **0** | **0** | **289** | **0** |

Files with a skip: 55.

Files with a failure, an error, an escaped warning or no RESULT line: 0 



## OpenBLAS (user approved; chain3 `ob`, 01:35 to 01:44)

`dev/ciharden-openblas.sh` 0.3.26 and 0.3.32 `lapack` built (no denial). The 60 files of
`dev/rel080-log/ob.jobs` (the 51 changed or new, the 8 `test-ce-*.R`,
`test-conditional-smooths.R`), 4 threads, every `blas:` line naming the copy, every
`lib:` line rellib-r19: 60 of 60 ran on each. Two files did not pass, the same on both:

- `test-nl-start.R:321` (nlstart's): b2 at the 1e-6 scale -0.06021 against -0.06022,
  1.7e-4 relative, against a tolerance of 1e-4; the log-likelihood check (1e-6) holds.
  Now compared in units of b2's standard error (0.01 SE).
- `test-se-reparam.R:453`: an escaped warning, "function evaluation limit reached
  without convergence (9). A Newton step from the stop gains 6.81e-08 ... within
  grad_tol"; the fit is at the maximum (the test's 1e-7 log-likelihood check holds).
  The stop is now an allowed outcome.

Both edits in `dev/rel080-release-fix2.R`; rerun on reference BLAS and both OpenBLAS
copies (`r080-obfix-*`). The core `R CMD check` (chain3 `check`) built its tarball at
01:44, before these two test edits.
- Rerun after the edits: both files pass on ref, ob26 and ob32 (nl-start 55, se-reparam 79; 0 fail, 0 warn).

## R CMD check --as-cran (chain3 `check`, 01:44 to 02:31; `dev/rel080-check/<pkg>/check.log`)

`Status: 1 NOTE` on each of frmtmb 0.80.0 (tests 22 min), frmtmb.sample 0.27.0,
frmtmb.latent 0.7.3, frmtmb.ode 0.7.3, frmtmb.survival 0.2.2; the NOTE is the HTML
manual's math rendering (V8 not installed). The core tarball was built at 01:44, before
the two OpenBLAS test edits (test files only; they pass on all three BLAS setups).

## Sessions (chain3 `sessions`, 02:31 to 06:11)

Plain (`dev/release/rel0.80.0.log`):

tier rel0.80.0: REL0.80.0 ran 27 of 27
SESSION lines: 27; with fail, err or warn not 0: 0
lib: lines: 243; naming rellib-r19: 243
LEAK lines: 46; without a tag: 0 ; tagged as a Stan model load (rstan): 0 ; naming package:emmeans: 0
| suite | files | pass forward, reverse, preload | skip forward, reverse, preload |
|---|---|---|---|
| frmtmb | 304/317 | 23245, 23395, 23488 | 239, 239, 241 |
| frmtmb.coupling | 11/11 | 551, 551, 551 | 5, 5, 5 |
| frmtmb.eam | 33/33 | 1848, 1848, 1848 | 5, 5, 5 |
| frmtmb.latent | 12/12 | 687, 687, 687 | 2, 2, 2 |
| frmtmb.learn | 16/16 | 474, 474, 474 | 13, 13, 13 |
| frmtmb.ode | 11/11 | 549, 549, 549 | 1, 1, 1 |
| frmtmb.sample | 81/86 | 3967, 3967, 3967 | 18, 18, 18 |
| frmtmb.survival | 15/15 | 552, 552, 552 | 0, 0, 0 |
| frmtmb.spatial | 8/8 | 478, 478, 478 | 6, 6, 6 |

Gated (`rel0.80.0g`):

tier rel0.80.0g: REL0.80.0G ran 27 of 27
SESSION lines: 27; with fail, err or warn not 0: 0
lib: lines: 243; naming rellib-r19: 243
LEAK lines: 194; without a tag: 0 ; tagged as a Stan model load (rstan): 142 ; naming package:emmeans: 0
| suite | files | pass forward, reverse, preload | skip forward, reverse, preload |
|---|---|---|---|
| frmtmb | 317/317 | 26430, 26580, 26673 | 0, 0, 2 |
| frmtmb.coupling | 11/11 | 557, 557, 557 | 0, 0, 0 |
| frmtmb.eam | 33/33 | 1872, 1872, 1872 | 0, 0, 0 |
| frmtmb.latent | 12/12 | 693, 693, 693 | 0, 0, 0 |
| frmtmb.learn | 16/16 | 559, 559, 559 | 0, 0, 0 |
| frmtmb.ode | 11/11 | 551, 551, 551 | 0, 0, 0 |
| frmtmb.sample | 86/86 | 4334, 4334, 4334 | 0, 0, 0 |
| frmtmb.survival | 15/15 | 552, 552, 552 | 0, 0, 0 |
| frmtmb.spatial | 8/8 | 561, 561, 561 | 0, 0, 0 |

## Docs (chain3 `docs`, 06:11 to 06:23; `dev/rel080-log/docs-check.txt`)

Roxygen again over all nine: 0 files moved, so the `devtools::document()` diff is empty.
Installed code against the tree: 0 functions differ and 0 absent (frmtmb 2395, sample 302, latent 67, coupling 24, eam 194, learn 47, ode 38, spatial 44, survival 49).
`pkgdown::check_pkgdown()`: OK on all nine. `dev/release/build-docs.R`: 9 of 9 sites, DOCS ok, 10.7 min (site folder deleted after).

## Floors (chain3 `floor`, 06:23 to 0x; `dev/rel080-log/floor.txt`)

- `logLik_penalized` in NAMESPACE: 0 at 3cd2a0c8, 3 in this tree.
- frmtmb.sample 0.27.0 built here, floor lowered, over 0.79.1 (rellib-r18): loads; 86 of
  86 files, 1 fails: `test-me-gr-new-level-draws.R:37` (posterior_epred() refuses a new
  `me(gr = )` level on 0.79.1). On 0.80.0: passes (gated tier). Converse: sample 0.26.0
  loads over 0.80.0.
- frmtmb.latent 0.7.3, same: 12 files, `test-hmm-starts.R` errs 16 times
  ("'logLik_penalized' is not an exported object"). Converse: 0.7.2 loads over 0.80.0.
- The six others (coupling, eam, learn, ode, spatial, survival), the release tree's
  tests on rellib-r18 alone: 94 of 94 files, 0 failures, 0 errors; so their floors
  need no change (ode and survival change only internally).

## The survey (chain3 `survey`, 07:16 to 07:23; `dev/rel080-survey.sh`)

- Frames: `dev/equiv-dump-compat.R` on rellib-r17 (0.79.0) into `dev/equiv-frame-out`
  (the 0.77.0 frame of rellib-r15 is gone) and on rellib-r19 into
  `dev/rel080-log/frame-r19`: 7518 core pairs each; 33 core pairs move status between
  them (25 conditional -> works, 7 refused -> works, 1 conditional -> refused).
- `FRMTMB_EQUIV_LIB="<rellib-r19>;<rellib-r18>" RECORD=1 bash dev/equiv-run.sh 3
  dev/rel080-equiv-out`: 172 of 172 logs read `lib: .../rellib-r19/frmtmb 0.80.0`,
  0 FAILED, `SURVEY END 0`.
- **Headline: 89.0, Wald [84.5, 93.6]**, the same as 0.79.0's; no defect 91.3,
  weak evidence out 85.3. `dev/rel080-equiv-delta.R` (`equiv-delta.txt`): 0 of 172
  pairs change verdict, level or agreement; `l1_const`, `l1_grad`, `l2_metric` and
  `l4_z` move on 0 pairs (relative 1e-6). The lanes expected none: fitedge p136
  (moves 4e-7), bigx all 172, nlstart 14 nl pairs, loglik 21 prior pairs.
- The summary's corrected-weights line now reads **NA**: 0.79.0's frame holds
  features `weights.tsv` does not weigh. `dev/rel080-frame-line.R` (na.rm, as
  `rel079-frame-line.R`): 88.9 on the r17 frame and 88.9 on the r19 frame (0.79.0:
  89.0 on 0.77.0's frame, 88.9 on 0.79.0's). Not changed in `equiv-summary.R`.

## Ledger (chain3 `ledger`, then rerun)

The first `ledger` step stopped: `dev/brms-suite/` (brms 2.23.0's tree) was absent from
the release tree and every lane tree. Copied from the main checkout (read only;
gitignored), then `dev/rel080-ledger.sh ledger` (exit 0) and `bin2` (running).
- Bin 1 (record mode, 14 of 14 files, 3 at a time), the ledger, the generated files,
  70 of 70 guards, bin 2: `dev/rel080-log/ledger-diff.txt`, `bin2-diff.txt`. Bin 1:
  pass 397, cannot transfer 48, divergence 49, defect 0, as at 0.79.1; 0 rows move.
  Bin 2: 308 rows, 0 move, 0 value checks move. Bins 1 and 2: 634 of 823 (77.0%),
  634 of 726 without the decisions (87.3%). The regenerated ledger files equal the
  tracked ones (git shows no change). No row flips.

## Close (2026-11-07 08:0x)

- Core `R CMD check --as-cran` again on the final tree (07:27 to 08:01): `Status: 1
  NOTE` (V8). A first try at 07:26 failed to build because I was deleting
  `dev/ciharden-out/` while `R CMD build` copied the tree; rerun clean. Check logs kept as
  `dev/rel080-log/check-<pkg>.log`, `00check-<pkg>.log`; `dev/rel080-check/`,
  `dev/ciharden-out/`, `dev/brms-suite/` (copy), the scratch and floor libraries deleted.
- `bash dev/rel080-commit.sh --dry-run`: DRY RUN DONE, index sha256 of 14 trees unchanged.
- Record: `dev/round-20261106.md`; `dev/round-handoff.md` has the 0.80.0 section.
- No process of the consolidation is running.

## Release review fixes (2026-11-07; `dev/reviews/2026-11-07-release.md`, FIX)

The coordinator's instruction: apply both majors and every minor as the review words
them, ship item 2 in 0.80.0, add a test for the multistart wording, rerun everything.

`dev/rel080-release-fix3.R`:
- item 1: `fit_assembled()` skips `nl_basis_descent()` on the refit of a probe's winner
  (`attr(start, "frmtmb_nl_probe")`); the reviewer's test in `test-nl-start.R`.
- item 2: `design_center()` as the review gives it (a `0 + f` design is centered outside
  the span of the constant), its roxygen; the reviewer's test in `test-design-offset.R`.
- item 6: the backlog's three missing items. The review keeps the bernoulli mm 55 item
  as it is once item 1 is applied (it is true again), so nothing is removed; the
  coordinator's message said "remove the stale item", the review's text wins.
- item 7: the checklist's commands (steps 1, 4, 10, 13) and step 6's wording.
- item 8: `na.rm` in `dev/equiv-summary.R`.
- item 11: the commit message names `rel080-release-fix2.R` and `-fix3.R`; a test of the
  multistart wording under a prior ("penalized log-likelihood", last value equals
  `logLik_penalized()`), skipping only if neither a restart nor the probe speaks.
- NEWS (items 2, 3, 4) through `dev/rel080-news.R`; the record (items 5, 8) by hand.
- Not applied: item 9 (NIT, `nl_local_message()`; filed) and item 10 (NIT).

Core reinstalled 2026-11-07 08:59:26 (`install3.txt`; roxygen moved no file).
`test-nl-start.R` 59 pass, `test-design-offset.R` 80 pass, 0 fail, 0 skip, 0 warn
(`r080-fix3`).

Running (09:0x): `dev/rel080-chain4.sh` (ref, plain, sessions; 2 R processes) and
`dev/rel080-chain5.sh` (nlstart 450 and fresh rows into `rr-fix3/`, bigx grid and
`bigx-prod.R`, OpenBLAS on `ob2.jobs`, core R CMD check after chain4, ledger with one
bin-1 file at a time, RECORD=1 survey with one job). Progress: `dev/rel080-log/chain.txt`.

### nlstart on the review-fixed build (`rr-fix3/`)

- 450-row scan: 0 short, 16 higher than 0.79.1, 0 lower (unchanged).
- Fresh rows: 11 short (5 silent), short > 1: 6 (3), multistart message on 87 rows, 0
  below 0.79.1. Against the lane's own logs: 0 rows lower, 2 higher; binomial mm 22
  warns again (-643.5014, the lane's value); bernoulli mm 55 is silent at 1.84 short
  again. One row's warning state still differs from the lane: bernoulli sig4 55 reads
  -189.696351, 2.9e-4 above the lane's warned -189.696640, with no warning, 2.8e-4
  below the best start (the review counted it "at the best value"; freshsum's best is
  -189.6961, so it lists the row as short by 2.8e-4). NEWS: "now 11 are, 5 with no
  warning" and "on 1 of the 87 fits", from this run.

### bigx on the review-fixed build

- Grid: 240 of 240 right, 11 warned, 0 silent wrong (unchanged).
- `dev/bigx-prod.R` (`rr-fix3/bigx-prod.txt`): `w1 * w2` and `w1 + I(w1^2)` kept and
  within 2.5e-7 at 1e4, within 1.9e-5 (code 8 at the maximum) at 1e5, capped with
  `frmtmb_design_capped` at 1e6 (no "rank deficient" message); raw polynomials on
  [1000, 1010] keep 5 columns with finite SEs, the capped warning and the message for
  the aliased ones. As the lane's `prod-after2.txt`.

### OpenBLAS on the review-fixed build (09:27 to 09:34)

`dev/rel080-log/ob2.jobs`: every `test-nl*.R` (9), `test-design-offset.R`,
`test-se-reparam.R` and the reviewer's 5 (`test-formula-cmc.R`,
`test-aliased-grouping.R`, `test-frame.R`, `test-edgecases.R`, `test-ordinal.R`), 16 files,
on the reference BLAS, 0.3.26 and 0.3.32 with LAPACK, 4 threads (copies rebuilt, no
denial): 16 of 16 each, every file 0 fail, 0 err, 0 skip, 0 warn; `blas:` lines name the
copies, `lib:` rellib-r19.

### Tiers on the review-fixed build (09:0x to 11:01)

Gated:

From `dev/ciharden-log/r080-ref-suite/` (509 jobs; 509 files with a RESULT line; 509 whose `lib:` line names rellib-r19 for the package and for frmtmb).

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 317 | 26439 | 0 | 0 | 0 | 0 |
| frmtmb.coupling | 11 | 557 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 33 | 1872 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 12 | 693 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 16 | 559 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 86 | 4334 | 0 | 0 | 0 | 0 |
| frmtmb.spatial | 8 | 561 | 0 | 0 | 0 | 0 |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |
| **total** | **509** | **36118** | **0** | **0** | **0** | **0** |

Files with a skip: 0.

Files with a failure, an error, an escaped warning or no RESULT line: 0 

Plain:

From `dev/ciharden-log/r080-plain-suite/` (509 jobs; 509 files with a RESULT line; 509 whose `lib:` line names rellib-r19 for the package and for frmtmb).

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 317 | 23254 | 0 | 0 | 239 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1848 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 12 | 687 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 86 | 3967 | 0 | 0 | 18 | 0 |
| frmtmb.spatial | 8 | 478 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |
| **total** | **509** | **32360** | **0** | **0** | **289** | **0** |

Files with a skip: 55.

Files with a failure, an error, an escaped warning or no RESULT line: 0 

### Sessions on the review-fixed build (11:01 to 15:54)

tier rel0.80.0: REL0.80.0 ran 27 of 27
SESSION lines: 27; with fail, err or warn not 0: 0
lib: lines: 243; naming rellib-r19: 243
LEAK lines: 46; without a tag: 0 ; tagged as a Stan model load (rstan): 0 ; naming package:emmeans: 0
| suite | files | pass forward, reverse, preload | skip forward, reverse, preload |
|---|---|---|---|
| frmtmb | 304/317 | 23254, 23404, 23497 | 239, 239, 241 |
| frmtmb.coupling | 11/11 | 551, 551, 551 | 5, 5, 5 |
| frmtmb.eam | 33/33 | 1848, 1848, 1848 | 5, 5, 5 |
| frmtmb.latent | 12/12 | 687, 687, 687 | 2, 2, 2 |
| frmtmb.learn | 16/16 | 474, 474, 474 | 13, 13, 13 |
| frmtmb.ode | 11/11 | 549, 549, 549 | 1, 1, 1 |
| frmtmb.sample | 81/86 | 3967, 3967, 3967 | 18, 18, 18 |
| frmtmb.survival | 15/15 | 552, 552, 552 | 0, 0, 0 |
| frmtmb.spatial | 8/8 | 478, 478, 478 | 6, 6, 6 |

tier rel0.80.0g: REL0.80.0G ran 27 of 27
SESSION lines: 27; with fail, err or warn not 0: 0
lib: lines: 243; naming rellib-r19: 243
LEAK lines: 194; without a tag: 0 ; tagged as a Stan model load (rstan): 142 ; naming package:emmeans: 0
| suite | files | pass forward, reverse, preload | skip forward, reverse, preload |
|---|---|---|---|
| frmtmb | 317/317 | 26439, 26589, 26682 | 0, 0, 2 |
| frmtmb.coupling | 11/11 | 557, 557, 557 | 0, 0, 0 |
| frmtmb.eam | 33/33 | 1872, 1872, 1872 | 0, 0, 0 |
| frmtmb.latent | 12/12 | 693, 693, 693 | 0, 0, 0 |
| frmtmb.learn | 16/16 | 559, 559, 559 | 0, 0, 0 |
| frmtmb.ode | 11/11 | 551, 551, 551 | 0, 0, 0 |
| frmtmb.sample | 86/86 | 4334, 4334, 4334 | 0, 0, 0 |
| frmtmb.survival | 15/15 | 552, 552, 552 | 0, 0, 0 |
| frmtmb.spatial | 8/8 | 561, 561, 561 | 0, 0, 0 |

### Check, ledger and survey on the review-fixed build (15:55 to 16:5x)

- Core `R CMD check --as-cran` (15:55 to 16:27): `Status: 1 NOTE` (V8).
- Ledger: chain5's run stopped at `dir.exists(suite)`, because I had deleted the copy of
  `dev/brms-suite/` after the first ledger. Copied again from main (read only) and
  kept; `REC_JOBS=3 bash dev/rel080-ledger.sh` (nothing else running): bin 1 14 files,
  ledger exit 0, 70 of 70 guards, bin 2 exit 0. Bin 1: 494 rows, 0 move; bin 2: 308
  rows, 0 move, 0 value checks move; 634 of 823 (77.0%), 634 of 726 (87.3%). The
  regenerated ledger files equal the tracked ones.
- Survey, `RECORD=1`, 1 job (frames dumped again first): 172 of 172 logs on
  `rellib-r19/frmtmb 0.80.0`, 0 FAILED, `SURVEY END 0`. Headline 89.0, Wald
  [84.5, 93.6]; 0 of 172 pairs change verdict, level or agreement, and no metric
  moves (`equiv-delta.txt`). The corrected-weights line now reads 88.9 (the `na.rm`
  fix). No pair moves.

### Close of the review fixes (2026-11-07 17:xx)

- Docs on the review-fixed build: roxygen moved 0 files; installed code 0 functions
  differ (all nine); `check_pkgdown()` OK on all nine; 9 of 9 sites (deleted). The
  reviewer's stronger check (`dev/review-rel-installed.R`): 0 functions differ; 22
  load-time values differ, as in the review (`installed-strong.txt`).
- `bash dev/rel080-commit.sh --dry-run`: DRY RUN DONE, rel080 779 paths, index sha256 of
  14 trees unchanged.
- Record regenerated (`dev/round-20261106.md`: counts blocks, totals, sessions, survey,
  section "After the release review"). `dev/ciharden-out/` and `dev/rel080-check/`
  deleted; `dev/brms-suite/` (copy from main, gitignored) kept, since the ledger needs it.
- No process of the consolidation is running.

## Short check of the release fixes: SHIP (2026-11-07)

Its two text edits and its NIT, text only: NEWS names the fifth silent fresh row (a
four-parameter logistic under `bernoulli()` 2.8e-4 short; `dev/rel080-news.R`); the
backlog's ridge item adds the saddle-restart refit (sig4 55) and a NIT on
`design_center()` and `0 + w1 * w2` (`dev/rel080-release-fix4.R`, named in the commit
message). NEWS parses with `tools:::.build_news_db_from_package_NEWS_md` for core and the
four changed extensions (first versions 0.80.0, 0.27.0, 0.7.3, 0.7.3, 0.2.2).
