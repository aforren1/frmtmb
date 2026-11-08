# Lane ci: the open CI items (round of 2026-11-06, base 3cd2a0c8)

Worktree `frmtmb-wt-ci`, library `wt-ci-lib`. Logs are local, under
`dev/ci-gh-out/` (gitignored); every number below is pasted from them.
GitHub was read with `gh run list`, `gh run view` and the read-only
`gh api .../actions/jobs/<id>` (step timings). Nothing was re-run,
dispatched, cancelled or pushed.

## Last items after the short check (verdict SHIP)

Applied by `dev/ci-cache-edit.py` (each old text must occur once; CRLF
kept, 16 of 16 workflows CRLF after).

- stan-next (`ubuntu-next.yaml`): job limit 150 to 180, the push job's.
  At 150, a slow mirror's 41 minutes of setup plus the identity tier's
  120-minute step limit (161) would cancel the job first, and
  `!cancelled()` would then skip the gp and methods tiers. The stale
  "Wider than the push job's 90 minutes" comment is rewritten, and the
  "Differences from it" comment no longer lists the time limit, which
  is now the same.
- The record step of ubuntu-next's core job has
  `continue-on-error: true`, so its failure cannot skip the check and
  the gate.
- The Stan cache in `brms-likelihood.yaml` and in stan-next is split
  into `actions/cache/restore@v4` (`id: stancache`) and a last step
  `actions/cache/save@v4` with
  `if: always() && steps.stankey.outcome == 'success' &&
  steps.stancache.outputs.cache-hit != 'true' &&
  hashFiles('dev/stan-cache/*.rds') != ''`. So a job that fails or hits
  a limit keeps what it compiled; an exact hit saves nothing, as
  before; a job that compiled nothing does not try to save.

  Can a partial cache overwrite a good one? Not by key: GitHub never
  overwrites a cache key, and a save is attempted only when the primary
  key missed. But saved under the primary key, a partial cache would be
  an exact hit on every later run, and the programs it lacks would be
  compiled every time and never saved. So a save from a job whose
  `job.status` is not `success` goes under the primary key plus
  `-partial-<run id>`. The restore-keys prefix matches that key, and
  `actions/cache` restores the newest prefix match, so the next run
  starts from the partial set (a superset of what it restored), and a
  successful run saves the full set under the primary key.

  Residual risks. (1) A truncated RDS from a killed write: CLOSED, see
  the next section. (2) Each failed or timed-out run adds one partial
  entry to the repository's 10 GB cache budget; GitHub evicts the least
  recently used. Not changed.
- Checks: every workflow parses with pyyaml; `dev/review-ci-yaml.py .`
  reports `gates: 12` (`ci-gh-out/review-yaml-final.log`); actionlint
  1.7.12 gives the same 2 `runner-label` findings and nothing else, so
  it accepts `job.status` in `with:` and `hashFiles()` in `if:`
  (`ci-gh-out/actionlint-final.log`). Whether the save step runs after a
  step timeout, and whether the next run restores the partial key, is
  CI only.

## Atomic writes to the Stan cache (coordinator's last item)

Every test code path that writes a compiled program to the Stan cache
now calls `stan_cache_write(mod, path)`: `saveRDS()` to
`<name>-XXXX.part` in the same directory, then `file.rename()` to the
final path, with a "why" comment. Applied by `dev/ci-atomic-edit.py`
(each old text once; CRLF kept). Ten sites in eight files of four
packages; every `saveRDS()` of a compiled program was one of them:

    tests/testthat/helper-brms.R                         defines it; brms_stan_model()
    extensions/frmtmb.learn/tests/testthat/helper-stan.R    defines it; ln_stan_model()
    extensions/frmtmb.spatial/tests/testthat/helper-brms.R  defines it; brms_stan_model()
    extensions/frmtmb.sample/tests/testthat/helper-brms.R   defines it
    extensions/frmtmb.sample/tests/testthat/test-default-priors-fit.R  2 sites
    extensions/frmtmb.sample/tests/testthat/test-mo-coding.R           2 sites
    extensions/frmtmb.sample/tests/testthat/test-loglik-autocor.R      1 site
    extensions/frmtmb.sample/tests/testthat/test-mi-draws-brms.R       1 site

Windows. On R 4.6.1 here, `file.rename()` over an existing file
succeeds and replaces it (measured: TRUE, target holds the new content,
no warning). It can still fail while another process holds the target
open. The helper then drops its copy when the target exists (a whole
program another run wrote) and stops only when the rename failed and
there is no target. The rename's warning is suppressed only on that
call, and its result is then checked.

By construction (`dev/ci-atomicprobe.R`, log `ci-gh-out/atomicprobe.log`,
a fresh compile into an empty cache, lane core):

    write: files 6a164bc061373577609fc4aa539ad7b8.rds
    crash: error 'killed mid-write'; final name exists: FALSE; .part files: 0
    crash, old write: final name exists: TRUE; readRDS: TRUE (TRUE = error)
    exists: warnings 0; replaced: TRUE; .part files: 0
    read: class stanmodel ; draws 10

"crash" makes `saveRDS()` write 100 bytes and stop: under the new
helper nothing appears under the final name; under the old plain write
a truncated `.rds` appears and fails `readRDS()`. "read" loads the
program in a second process and samples from it, the cross-session
load that the cache exists for.

The test files that can reach `stan_cache_write()`, found by
`dev/ci-cachecallers.py` (helper functions that call it, transitively,
then the test files that call any of them): 36 in core, 1 in
frmtmb.learn, 2 in frmtmb.spatial, 4 in frmtmb.sample. Each ran one per
process, 3 at a time (`dev/ci-atomic-run.sh`), with NOT_CRAN=true,
FRMTMB_BRMS_FIT_TESTS=true and FRMTMB_STAN_CACHE at the worktree's
`dev/stan-cache`. The core is the lane build (36 `lib:` lines read
wt-ci-lib/frmtmb); the extensions are rellib-r18's builds, which run
the worktree's changed test files. Every log has a RESULT line and
`R.home: C:/PROGRA~1/R/R-46~1.1`. Totals over 43 files: pass=3832
fail=0 error=0 skip=0 warning=0 (`ci-gh-out/atomic-results.txt`, one
line per file, every one fail=0 error=0 skip=0 warning=0). Largest:
test-brms-methods.R 1013, test-brms-likelihood.R 540,
test-latres-brms.R 281, test-nl-stan-funs.R 162, test-loglik-autocor.R
138 (sample), test-sar.R 126 (spatial), test-stan-identity.R 71
(learn), test-brms-likelihood-car.R 67 (spatial), test-mo-coding.R 50
(sample), test-default-priors-fit.R 28 (sample), test-mi-draws-brms.R
12 (sample); the 13 BCM chapter files 12 to 52 each.

The write path was exercised at scale, not only read: the seeded cache
held 238 programs and holds 434 after the run, so 196 programs were
compiled and written through `stan_cache_write()` and then read back by
the tests that use them. No `.part` file is left in it.

Left: a process killed outright (not an R error) skips `on.exit()`, so
a `<name>-XXXX.part` can stay in the cache directory. No reader opens a
`.part` file, and CI's save step uploads the directory as it is, so at
worst it carries a dead file. `helper-brms.R` is one of the five files
hashed into brms-tiers' cache key, so the first run after this change
misses the exact key, restores by prefix, and saves under a new key.

## Fix round (review `dev/reviews/2026-11-06-ci.md`, verdict FIX)

- MAJOR 1, `cran-lag.R` failed open: fixed; old exit 0 and new exit 1
  on the bad-mirror and offline cases, plus empty-index,
  unknown-Suggests, all-current (exit 0) and missing-library cases
  (item 4).
- Minor 4, one library: it reads a `;`-separated stack, default
  `.libPaths()` (item 4).
- Minor 2, amd64v3: recorded (item 1a); ubuntu-next prints the OpenBLAS
  package and its index; no local reproduction.
- Minor 3, brms-tiers 44 minutes: findings corrected; limits raised to
  180 (job) and 120 (identity tier step) (item 1b).
- Minor 1, ordinal nl headroom: filed in `dev/test-backlog.md`.
- NIT, gate's second red line: kept, with the reason (item 3). The
  other NITs: the 1.5 factor (item 1b), NEWS says "tries to show",
  `ci-runfile.R` prints `R.home()`, the tarball NEWS (R CMD check).

## Summary

| Item | Verdict |
|---|---|
| 1 ubuntu-next red | Two causes. (a) 26.04 core check: `test-mo-coding.R:236`, a bitwise identity between a 3-row and a 4-row prediction, 1 ulp apart under OpenBLAS 0.3.32. Fixed in the test. (b) "Stan tiers, ubuntu-24.04" cancelled at its 90-minute job limit: every weekly run starts from an evicted cache. Limit raised to 150. Do NOT unpin yet: see "Unpinning". |
| 2 hurdle_negbinomial gradient warning | It already escaped `allow_warnings()` at 0.79.1 (WARN 1 in the ubuntu release log); an escaped warning does not fail a test. Cause: Poisson data, so the NB shape's maximum is at infinity and the stopping point is rounding. Fixed by giving the five NB families and `student()` data with an interior maximum. |
| 3 skip gate after a failed check | Done in all 11 gate steps (R-CMD-check, 8 extension checks, 2 coverage jobs), plus a new gate in ubuntu-next's core job. |
| 4 CRAN lag script | `dev/release/cran-lag.R`. Fixed in the fix round: it failed open when CRAN was unreachable (review, MAJOR 1); it now fails closed and reads a library stack. No release checklist exists under `dev/release/` (or anywhere in the repository): needs the user. |
| 5 pkgcheck | 4 failing checks at 0.79.1; 2 closed by 3cd2a0c8, 1 (`<<-`) is superassign's, 1 (CI) should close if test-coverage is green on the release push. The 6-hour default limit is safe today with about 50 minutes to spare in the worst case, but the margin shrinks each round. |
| 6 other red or flaky | See section 6. brms-tiers can exceed 90 minutes on a cold cache and a slow mirror: limits raised in the fix round. |

## 1. ubuntu-next

Run 41454170880 (schedule, 2026-11-02, d5d0056d), four jobs:

    ubuntu-26.04 (release)     failure    125306525470
    ubuntu-24.04 (release)     success    125306525624
    Stan tiers, ubuntu-26.04   success    125306525845
    Stan tiers, ubuntu-24.04   cancelled  125306525887

### 1a. The 26.04 core check (log `dev/ci-gh-out/un-26-core.log`)

    [ FAIL 1 | WARN 0 | SKIP 234 | PASS 20197 ]
    ── Failure ('test-mo-coding.R:236:3'): new data are coded from the fitted minimum ──
    Expected `fitted(f1, newdata = data.frame(k = 2:4))` to be identical to `p1[2:4, ]`.
     [3] 6.63854071086290976 - 6.63854071086290887 [3]
     [9] 6.40285698705515749 - 6.40285698705515660 [9]
    [12] 6.87422443467066202 - 6.87422443467066113 [12]

Only the k = 4 row's Estimate differs, by 8.9e-16, which is 1 ulp at
6.6; Q2.5 and Q97.5 follow it and Est.Error is identical. The image:
R 4.6.2, Ubuntu 26.04.1, BLAS and LAPACK
`openblas-pthread/libopenblasp-r0.3.32.so`, gcc 15.2.0. The 24.04 column
(OpenBLAS 0.3.26) passed the same commit. The line came in with
bb221eb3 (2026-10-27), so the 10-19 and 10-26 runs never held it.

It did NOT reproduce here. `dev/ci-moprobe.R` (seed 5, the test's data)
under `dev/ciharden-openblas.sh 0.3.32 lapack`, with
`OPENBLAS_CORETYPE` unset, Haswell, Zen, Sandybridge and Nehalem, and
under the reference BLAS: identical in every arm, matprod "default" and
"internal", and a direct 3-row against 4-row `X %*% b` probe differs in
0 of 2000 draws (log `moprobe-ob332.log`, `moprobe-ref.log`). Whether the
Windows OpenBLAS build honors `OPENBLAS_CORETYPE` was not checked. So
the mechanism (which kernel, which operation) is not established; the
runner's CPU was not recorded. ubuntu-next now records it (below).

Lead, from the review (minor 2), checked in the log: 26.04 installs
Ubuntu's x86-64-v3 build, not a generic one. `un-26-core.log`:

    Preparing to unpack .../3-libopenblas0-pthread_0.3.32+ds-5_amd64v3.deb ...
    Preparing to unpack .../4-libopenblas0_0.3.32+ds-5_amd64v3.deb ...

(334 lines of that log name `amd64v3`; the apt indexes are
`resolute-updates/main amd64v3 Packages` and so on.) Ubuntu compiles an
amd64v3 package for x86-64-v3, which includes AVX2 and FMA, so the
compiler may contract a multiply and an add into one FMA in OpenBLAS's
C code. `dev/ciharden-openblas.sh` uses upstream's generic Windows
build, `OpenBLAS-0.3.32-x64.zip`, which is a different binary. That can
explain why it did not reproduce. Upstream publishes no Windows build
for x86-64-v3; a local reproduction would mean building OpenBLAS from
source with `-march=x86-64-v3`, which is not cheap, so it was not
tried. 24.04 has no amd64v3 packages. The record step now also prints
`dpkg -l 'libopenblas*'` and `apt-cache policy libopenblas0-pthread`,
whose index line names the variant.

The fix (`tests/testthat/test-mo-coding.R`): the bitwise identity keeps
p1's shape, `fitted(f1, newdata = data.frame(k = c(2L, 2:4)))[2:4, ]`
against `p1[2:4, ]` (4 rows each, the same rows at the same positions;
the new data's minimum is still 2, so a code from the new data's own
minimum still changes these rows). The 3-row prediction is held to
`max(|p3 - p1[2:4, ]| / |p1[2:4, ]|) < 64 * .Machine$double.eps`, the
file family's ulp bound (`test-arma-na-newdata.R`).

Seen to fail on the defect it guards (`dev/ci-modefect.R`, base build
rellib-r18, `mo_codes()` mocked to code from the new data's minimum, the
behavior bb221eb3 fixed; log `modefect.log`):

    real coder:
      identity on 4 rows: pass; 3-row bound: pass (ratio 0)
    coder from the new data's own minimum:
      identity on 4 rows: FAIL; 3-row bound: FAIL (ratio 0.551)

The old line's failure on 26.04 cannot be shown here (it did not
reproduce); the 26.04 log above is the record.

### 1b. Stan tiers, ubuntu-24.04: cancelled at 90 minutes

Log `un-24-stan.log`: the identity tier passed, `[ FAIL 0 | WARN 0 |
SKIP 0 | PASS 494 ]`, `Duration: 5296.4 s`, then the job was cancelled
during the methods tier at 15:05:51 (job start 13:34:40). The
compiled-model cache missed on every weekly run (GitHub evicts a cache
not used for 7 days, and the schedule is 7 days):

    run          image  cache                       identity tier
    39240311179  24.04  Cache not found (key+prefix) 2936.0 s
    39240311179  26.04  Cache not found              4692.0 s
    40313144422  24.04  Cache not found              3826.3 s
    40313144422  26.04  Cache not found              4687.8 s
    41454170880  24.04  Cache not found              5296.4 s
    41454170880  26.04  Cache not found              3956.1 s

(26.04 job totals: 85 min on 10-19 and 10-26, 72 min on 11-02, against
a 90-minute limit.) `stan-next`'s `timeout-minutes` is now 150. The
worst job was cut off at 90; from its identity tier of 88.3 minutes plus
about 7 to 10 minutes of setup and methods tier (the 26.04 jobs), its
full length would have been about 98, so 150 is about 1.5 times that.
The identity tier also has a step limit of 120 minutes now, as in
brms-likelihood.yaml (below).

Correction (review, minor 3). The first report said the push job
`brms-likelihood.yaml` keeps 90 minutes because "its cache hits, and
its last five runs took 5 to 10 minutes". That was wrong for the run on
3cd2a0c8: brms-tiers 42214681196 took 44 minutes, and its steps
(`gh api .../runs/42214681196/jobs`) were setup-r 25 minutes and
setup-r-dependencies 15 minutes (the slow mirror), then the identity and
methods tiers about 1 minute each from a warm cache. Its cache prefix
carries the rstan version, so a new rstan release gives a cold cache,
and ubuntu-next measured a cold identity tier on 24.04 at 49 to 88
minutes. A cold cache and a slow mirror together: about 41 + 88 + 3 =
132 minutes, which exceeds 90. So brms-likelihood.yaml now has a job
limit of 180 and a step limit of 120 on the identity tier (88 at the
slowest measured, a factor of 1.36), the test-coverage pattern: the
step limit catches a wedged run, the job limit leaves room for setup.
The job limit is 180 rather than 150 because 41 + 120 + the gp and
methods tiers exceeds 150. stan-next keeps 150 (it runs weekly and
reports only) with the same 120 on its identity tier.

### Also changed in ubuntu-next.yaml

- "Record the numerical libraries" prints the CPU model
  (`/proc/cpuinfo`) and runs with `OPENBLAS_VERBOSE=2`, which makes a
  DYNAMIC_ARCH OpenBLAS print its kernel ("Core: ...") at load. Not
  verified that Ubuntu's build prints it: CI only. In the fix round it
  also prints `dpkg -l 'libopenblas*'` and `apt-cache policy
  libopenblas0-pthread` (each followed by `; true`, so a missing tool
  cannot fail the step).
- The core job gets `id: check` and the skip gate, as R-CMD-check.yaml,
  whose steps its header says it copies.

### Unpinning

Do not unpin yet. On the fix's reasoning, the only 26.04 core failure
of 11-02 is fixed and the 26.04 Stan column passed (10-19, 10-26,
11-02), but the tree has changed since d5d0056d (0.79.0, 0.79.1, and
this round's merge) and nothing has run on 26.04 since; my OpenBLAS
0.3.32 shim did not reproduce the one failure that did happen, so a
local green on it is weak evidence about 26.04. Dispatch ubuntu-next
after the release push (or wait for Monday's schedule); when both
26.04 columns are green, change the pins ("Pinned 2026-10-06" in
R-CMD-check, test-coverage, pkgcheck, pkgdown and the others).

## 2. test-nl-family.R:93, hurdle_negbinomial

0.79.1 run 42116483453, job 127557370985 (ubuntu-24.04 release), log
`rcc-0791-ubrel.log`:

    [ FAIL 0 | WARN 1 | SKIP 270 | PASS 22261 ]
    ══ Warnings ══
    ── Warning ('test-nl-family.R:93:5'): the nl body reaches the likelihood of hurdle_negbinomial ──
    Large maximum absolute gradient at the optimum (0.00201 at shape_(Intercept)): the curvature there is unusable, ...
      2.   ├─frmtmb:::allow_warnings(...) at test-nl-family.R:65:3

It recurred on the base of this round: R-CMD-check run 42214681186
(3cd2a0c8), on both Ubuntu jobs and on neither of the others:

    ubuntu-24.04 (devel):     [ FAIL 0 | WARN 1 | SKIP 270 | PASS 22261 ] Warning ('test-nl-family.R:93:5') ... hurdle_negbinomial
    ubuntu-24.04 (release):   [ FAIL 0 | WARN 1 | SKIP 270 | PASS 22261 ] Warning ('test-nl-family.R:93:5') ... hurdle_negbinomial
    windows-latest (release): [ FAIL 0 | WARN 0 | SKIP 270 | PASS 22264 ]
    macos-latest (release):   [ FAIL 0 | WARN 0 | SKIP 270 | PASS 22264 ]

(The 3 passes fewer on Ubuntu were not traced: the check reporter does
not print per-file counts.)

So it already escaped `allow_warnings()` (its allowed list is
"Standard errors", "did not report convergence", "standard error"). An
escaped warning is counted in WARN; it does not fail the test or R CMD
check. The assertion itself had a wide margin: the file's ratio
`|nl - lin| / (lin - int)` was at most 5.06e-8 against the bound 1e-3
on every arm below.

Cause: the five negative binomial families were fitted to Poisson
counts (and hurdle_negbinomial to `cnt + 1`), so the shape's maximum is
at infinity: log-shape stopped between 16.1 and 22.1. Where it stops,
and which warning it gives, is the platform's rounding.
`dev/ci-nlfamprobe.R` fits the file's 41 families x 3 models from the
file's own builders, on the base build, under three arms (reference
BLAS, OpenBLAS 0.3.26 lapack, 0.3.32 lapack):

    before: 123 fits, logLik differs across arms in 11, warnings differ in 7: negbinomial/nl negbinomial2/lin negbinomial2/nl zero_inflated_negbinomial/lin zero_inflated_negbinomial/nl hurdle_negbinomial/lin hurdle_negbinomial/nl
    after3: 123 fits, logLik differs across arms in 1, warnings differ in 0:

Before, hurdle_negbinomial gave "Standard errors are not available"
on the reference arm and "Optimizer did not report convergence: false
convergence (8)" on both OpenBLAS arms (`hnb-*.log`); ubuntu gave the
gradient warning. After: `student()` was the only other family with a
dpar beyond 8 on the log scale (intercept-only nu at 19.7, the same on
all arms); `xbeta` stays at its boundary by design (no exact 0 or 1),
with its own allowed warning, the same on every arm. The one logLik
that still differs is zero_one_inflated_beta/lin, in the 10th decimal.

The fix (`tests/testthat/test-nl-family.R`, `nl_fam_response()`):
negbinomial, negbinomial2, nbinom1: `rnbinom(n, size = 2, mu =
exp(eta))`; hurdle_negbinomial: the same plus 1 where not zero;
zero_inflated_negbinomial: size 8, because at size 2 the count's own
zeros absorbed the inflation and the intercept-only fit sent zi to
-16.8 (`dev/ci-zinbprobe.R`, `zinbprobe.log`; size 8 gives zi -1.816
and -1.695, truth logit(0.2) = -1.39); student: `eta + rt(n, df = 4)`.
The draws are made inside the family's branch of the switch, after
`nl_fam_data()`, so no other family's data change (their logLik values
are bitwise equal before and after in the probe). Allowed list
unchanged, so a real gradient warning still escapes and shows.

After, with the NB families (generated from `nlfam-after3-*.tsv`):
every NB and student dpar between 0.175 and 2.41, no warning, logLik
bitwise equal across the three arms; largest ratio in the file 2.19e-09
(bernoulli) on all three arms.

## 3. Skip gate after a failed check

Every call of `.github/scripts/skip-gate.sh` (11 steps in 10 files):
the step that writes testthat.Rout got an id (`check` on
`r-lib/actions/check-r-package@v2`, `measure` on "Measure coverage"),
and the gate got

    if: ${{ !cancelled() && steps.<id>.outcome != 'skipped' }}

so it runs after a pass or a failure of that step, and not after a
cancel or when that step never ran. (Corrected in the fix round: the
first text said a failed setup would add a second "no testthat.Rout"
error. It cannot: a failed setup skips the check, so the gate does not
run either.) A failed check leaves `testthat.Rout.fail`,
which the gate's `testthat.Rout*` glob reads. A step that hits its own
`timeout-minutes` fails rather than cancels the job, so the gate runs
and, on a cut-off output, reports "no testthat summary line", which is
true. Applied by `dev/ci-gate-edit.py` (refuses a file already edited,
and a count of gates or anchors other than expected); the coverage
comment was reworded after it ran. Line endings kept CRLF.

Every workflow parses (`uv run --with pyyaml`), and the parse lists 12
gates, each with its `if` and a step of the matching id. actionlint
1.7.12 (`uvx --from actionlint-py`, shellcheck and pyflakes off) reports
2 findings, both `label "ubuntu-26.04" is unknown` in ubuntu-next.yaml:
its label list predates the image; the label exists (the 11-02 run used
it). The base workflows (3cd2a0c8) give the same 2. Logs
`actionlint-after.log`, `actionlint-base.log`.

Only CI can show the gate running after a red check.

The second red line (review, NIT). A check that dies before its tests
(an install error, an examples error) leaves no testthat.Rout, so the
gate now adds "no testthat.Rout under check, so nothing was checked" as
a second red step under the red check. I left it. The line is true,
and the job is red either way. The alternative, passing the check's
outcome to the gate and passing when the check failed and no output
exists, adds to the gate a branch that exits 0 with nothing checked,
which is the fail-open shape the gate exists to prevent; one mistake in
that condition hides a green check with no output. That is not a
trivial change, so I did not make it.

Fix-round recheck: the reviewer's `dev/review-ci-yaml.py .` reports
`gates: 12`, each with its `if`, a target of the matching id, no
`continue-on-error`, and every workflow CRLF (log
`ci-gh-out/review-yaml-rerun.log`). Step and job limits after the fix
round, from a pyyaml parse:

    brms-likelihood.yaml  brms-tiers          identity tier step 120, job 180
    test-coverage.yaml    test-coverage       Measure coverage step 105, job 150
    test-coverage.yaml    extension-coverage  Measure coverage step 60, job 150
    ubuntu-next.yaml      stan-next           identity tier step 120, job 150
    pkgdown.yaml          build               job 90 (18 minutes on 3cd2a0c8; not changed)

actionlint after the fix round: the same 2 `runner-label` findings and
no others (`ci-gh-out/actionlint-fix.log`).

## 4. CRAN lag script

### Fix round: it failed open (review, MAJOR 1), now fixed

`available.packages()` only warns when it cannot read an index, and
returns an empty table. Every package then read "not in any
repository", which the script did not count, and it ended "nothing
lags" with exit 0. Seen on the first version (a copy kept at
`ci-gh-out/oldlag/`, run against copies of the DESCRIPTIONs; logs
`lag-old-*.log`):

    old, mirror http://127.0.0.1:9:          exit 0, closure 78, "nothing lags"
    old, proxy to 127.0.0.1:9 (both repos):  exit 0, closure 73, "nothing lags"
    old, a fake Suggests zzzNoSuchPkg:       "5 package(s) lag or are missing"
                                             = 4 LAGS + drmTMB; zzzNoSuchPkg
                                             listed but not counted

The fix:
- each `available.packages()` call runs under
  `withCallingHandlers(warning = stop)`, so an unreadable index stops
  the run with exit 1 (the Windows binary index, information only,
  prints "binary index not read" and goes on);
- an empty CRAN index stops the run;
- every status other than "current" counts toward exit 1, "not in any
  repository" included, and a version string `package_version()` cannot
  read is its own status, "unreadable version", not a silent NA;
- the first argument is a `;`-separated library stack, default
  `.libPaths()` (review, minor 4); the header names the default; the
  `where` column prints the library a package came from;
- the repository root falls back to the working directory when
  `--file` names another script (the reviewer's harness sources it).

The cases, on the final script (`ci-gh-out/lag-final-summary.txt`):

    badmirror exit=1 | Error: URL 'http://127.0.0.1:9/src/contrib/PACKAGES.rds': status was 'Could not connect to server'
    offline exit=1 | Error: URL 'https://cloud.r-project.org/src/contrib/PACKAGES.rds': status was 'Could not connect to server'
    emptyindex exit=1 | Error: the CRAN index at file:///C:/Users/adf44/source/r/frmtmb-wt-ci/dev/ci-gh-out/emptyrepo is empty
    fakesuggests exit=1 | closure: 286 packages from 74 direct ones of 9 DESCRIPTION files 6 package(s) lag, are missing or are in no repository
    default exit=1 | closure: 285 packages from 73 direct ones of 9 DESCRIPTION files 6 package(s) lag, are missing or are in no repository
    stack exit=1 | closure: 285 packages from 73 direct ones of 9 DESCRIPTION files 5 package(s) lag, are missing or are in no repository
    allcurrent exit=0 | closure: 285 packages from 73 direct ones of 9 DESCRIPTION files nothing lags
    nosuchlib exit=1 | Error: no such library: C:/nope

"offline" sets `http_proxy` and `https_proxy` to 127.0.0.1:9.
"emptyindex" is a `file://` repository with an empty PACKAGES.
"fakesuggests" adds `zzzNoSuchPkg` to core's Suggests in a copy of the
DESCRIPTIONs; it reads "not in any repository" and is counted (6 = 5
LAGS + 1). "allcurrent" uses a `file://` index written from the
installed stack itself (`dev/ci-lagrepo.R`), so the exit-0 path is now
exercised. "stack" is rellib-r18 then the user library: drmTMB now
reads current (from rellib-r18), and RTMBdist reads LAGS: CRAN
published RTMBdist 1.2.1 (source) between this lane's first run and
the fix round; the Windows binary is still 1.2.0. RTMBdist is the
package whose 1.2.0 release turned 0.79.0's CI red, so read its NEWS
before the final suites.

The reviewer's `dev/review-ci-cranlag.R` (which sources the lane's
script) first stopped at `file.exists(descs[[1L]])`: under `sys.source`
the root came from the reviewer's `--file`. The root fallback fixed it;
rerun (`ci-gh-out/review-cranlag-rerun.log`): lane exit status 1,
independent closure 285 from 73 direct against the lane's 285, nothing
in either difference, the version comparisons as before, and "library
dirs before/after: 407 407 newest mtime unchanged: TRUE".

### First version (kept for the record)

`dev/release/cran-lag.R [library] [mirror]`: the closure is Depends,
Imports, LinkingTo and Suggests of core and the eight extensions, plus
covr, rcmdcheck, roxygen2, testthat, srr, pkgstats, then recursive
Depends/Imports/LinkingTo; base packages and the frmtmb packages out;
non-CRAN packages from the DESCRIPTIONs' Additional_repositories
(RTMBode, kaskr.r-universe.dev). A package resolves to the library, then
`.Library`. The reference is the source version (what the Ubuntu
runners install); the Windows binary version is printed beside it.
Exit 1 when anything lags or is missing. It reads only.

Run on the user library, 2026-11-06 (log `cran-lag-user.log`):

    closure: 285 packages from 73 direct ones of 9 DESCRIPTION files
      current                280
      LAGS                   4
      MISSING                1
      not in any repository  0
      package  status installed   where repository  binary from
     charport    LAGS     0.1.1 library      0.2.0   0.1.1 CRAN
       drmTMB MISSING      <NA>    <NA>      0.7.0   0.7.0 CRAN
         gert    LAGS     2.4.1 library      2.5.0   2.5.0 CRAN
       igraph    LAGS     2.3.3 library      2.3.4   2.3.4 CRAN
         nlme    LAGS   3.1-171 library    3.1-172 3.1-171 CRAN

charport and nlme have no newer Windows binary yet. drmTMB is in
rellib-r18, not the user library. Spot check: RTMBdist 1.2.0, Matrix
1.7-6, mgcv 1.9-5 and RTMBode 1.0 read current. The first version's
"nothing lags" path was reachable and fired falsely (above); the claim
that it "was not exercised" was true but missed that.

No release checklist exists under `dev/release/` (it holds scripts and
.jobs files) or anywhere in the repository or its history; the one
mention, `dev/organizer-rules.md:119`, refers to one the organizing
session keeps. The nearest thing is the "Evidence required at
hand-back" list of each round's release brief in `round-briefs/`,
outside this worktree. The line to add, before the final suites:

    - Compare the library with CRAN: `Rscript dev/release/cran-lag.R`
      (exit 0, or each lagging package updated or explained).

## 5. pkgcheck, run 42116483431, job 127557370511 (03004e4d, 0.79.1)

Failing checks (`##[error]`, log `pkgcheck-0791.log`):

| Finding | Closed this round? |
|---|---|
| no documented return value: `[frmtmb-gp]` | Closed by 3cd2a0c8 (`@return` in R/gp-by.R; man/frmtmb-gp.Rd), the base of this round |
| no examples: `[frmtmb-gp]` | Closed by 3cd2a0c8 (`@examples`) |
| global assignment operator `<<-` | Lane superassign's; not touched here |
| fails continuous integration checks | Not closed by code. At 03004e4d R-CMD-check was green and test-coverage was CANCELLED at its 90-minute job limit (setup 44 min, then "Measure coverage" cancelled after 46 min). The workflow's own comment says pkgcheck scores the latest push runs of R-CMD-check and test-coverage; the pkgcheck log does not print the table, so this attribution is inferred. Closes when test-coverage is green on the release push; the 150/105 limits are this round's fix |

Passing or informational: package coverage 91.3%; srr: 120 standards
complied with, 10 N/A; R CMD check no errors, no warnings; "Some
goodpractice linters failed" (45 linters, informational); 123 function
names duplicated in other packages (mostly brms's own vocabulary, by
design); `\dontrun` in examples of `frm`, `frmtmb_register_aterm`,
`frmtmb_register_covstruct`, `frmtmb_register_frame_check`;
goodpractice "Add examples to all exported functions" names
`man/frmtmb-gev-xi-chart.Rd` (it has no `\examples`; not in pkgcheck's
failing list); "Chunk start/end indices failed sanity checks in
vignette" 5 times.

Timing. The job has no `timeout-minutes`, so GitHub's 6-hour maximum
for a hosted runner applies, and it cannot be raised. Step timings
(`gh api .../actions/jobs/<id>`), newest first:

    commit    wait (min)  pkgcheck step (min)  job
    03004e4d      89          154               4h06
    01a1f3da      93          143               3h59
    e548cd52      54          120               2h57
    b9d9d8f0      76          123               3h23
    d5d0056d      65          101               2h49
    ac4906bd      61          125               3h09
    30d5f354      48          120               2h51

Inside the pkgcheck step at 03004e4d: covr 87.6 min, rcmdcheck 48.7 min,
cyclocomp 7.0, lintr 5.8. The wait ends when test-coverage ends, and
that job is now capped at 150 minutes, so the worst case is about
4 (pull) + 150 + 154 = 308 minutes against 360. Safe today, with about
50 minutes to spare, but the pkgcheck step grew 101 to 154 minutes
from 0.78.0 to 0.79.1. If it grows by the same amount again, a slow
coverage day will hit the 6-hour cap. Not changed here; options for the
user: cap the wait at a total budget instead of 100 + 100 polls, or
accept that pkgcheck may then score a running coverage job as a failure.

## 6. Last run of every workflow on main, read 2026-11-06 12:41 UTC

Update at 13:32 UTC: every 3cd2a0c8 run finished green except pkgcheck
(still running): R-CMD-check 42214681186 (4 of 4), brms-tiers,
check-frmtmb-coupling, test-coverage 42214681230. The test-coverage job
had a fast mirror this time: setup-r 0.5 min, dependencies 4.3 min,
"Measure coverage" 82.0 min (12:09:48 to 13:31:48), and the skip gate
passed. 82 minutes against the new 105-minute step limit is a margin of
1.28; it was 74 minutes at 0.79.0. Watch it; raise it with the suite.

On 3cd2a0c8 (pushed 12:02 UTC) still running then: R-CMD-check,
brms-likelihood, check-frmtmb-coupling, pkgcheck, test-coverage. Done
and green on 3cd2a0c8: the eam, latent, learn, ode, sample, spatial and
survival checks, and pkgdown. Last completed of the rest:

- R-CMD-check 03004e4d: green on all four. The three before it
  (b9d9d8f0, e548cd52, 01a1f3da) were red on every platform: the 0.79.0
  defects fixed by 0.79.1 (dev/test-backlog.md, "Filed at the 0.79.1
  release"). 03004e4d had one escaped warning (item 2).
- test-coverage: cancelled on 03004e4d and 01a1f3da (slow apt mirror,
  44-minute setup, 90-minute job limit), red on e548cd52 and b9d9d8f0.
  The 3cd2a0c8 run still has the old 90-minute limit, because the new
  limits are uncommitted: it can be cancelled again for the same reason.
  Flaky on the mirror, not on the code.
- pkgcheck: red on every run since at least 8413abf9 (item 5).
- ubuntu-next: red on 11-02 (item 1). 10-12 was red on both core
  columns at 48c1e6c9, fixed by 47cf659d (green 10-19 and 10-26).
- brms-likelihood: green on its last five runs. Corrected in the fix
  round: they took 5 to 10 minutes before 3cd2a0c8, and 44 minutes on
  3cd2a0c8, 40 of them setup on the slow mirror (item 1b). The review
  found the same slow setup on two more jobs of that push:
  check-frmtmb-coupling (setup-r 26.2, dependencies 19.4, job 48 min)
  and the frmtmb.coupling coverage job (25.2, 18.6, 47 min). The core
  coverage job of 3cd2a0c8 ran 87.1 minutes in all against the old
  90-minute job limit.
- check-frmtmb-coupling: last completed 01a1f3da, green.
- absent-helper-identical (70476dc3) and scale-helper-identical
  (c84587f0): green; they run only on their paths.

## Tests run (one file per process, `dev/ci-runfile.R`)

    mo-before-ref    rellib-r18/frmtmb | RESULT test-mo-coding.R pass=164 fail=0 error=0 skip=0 warning=0
    mo-before-ob332  rellib-r18/frmtmb | RESULT test-mo-coding.R pass=164 fail=0 error=0 skip=0 warning=0
    mo-after-ref     rellib-r18/frmtmb | RESULT test-mo-coding.R pass=165 fail=0 error=0 skip=0 warning=0
    mo-after-ob326   rellib-r18/frmtmb | RESULT test-mo-coding.R pass=165 fail=0 error=0 skip=0 warning=0
    mo-after-ob332   rellib-r18/frmtmb | RESULT test-mo-coding.R pass=165 fail=0 error=0 skip=0 warning=0
    mo-lanelib       wt-ci-lib/frmtmb  | RESULT test-mo-coding.R pass=165 fail=0 error=0 skip=0 warning=0
    nl-after-ref     rellib-r18/frmtmb | RESULT test-nl-family.R pass=320 fail=0 error=0 skip=0 warning=0
    nl-after-ob326   rellib-r18/frmtmb | RESULT test-nl-family.R pass=320 fail=0 error=0 skip=0 warning=0
    nl-after-ob332   rellib-r18/frmtmb | RESULT test-nl-family.R pass=320 fail=0 error=0 skip=0 warning=0
    nl-lanelib       wt-ci-lib/frmtmb  | RESULT test-nl-family.R pass=320 fail=0 error=0 skip=0 warning=0

The "ob" arms are `dev/ciharden-openblas.sh <ver> lapack` copies of R
4.6.1 (Rblas 76 symbols and Rlapack 617 to OpenBLAS). Only test files
changed and no R/ code, so no extension suite can read the change; the
extensions were not run.

## R CMD check

Core, `dev/ci-check.sh` (R_LIBS = wt-ci-lib, rellib-r18, user library;
NOT_CRAN=true; pandoc and TinyTeX on PATH;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`), log `ci-gh-out/check.log`:

    * checking tests ...
      Running 'testthat.R' [64m]
     [64m] OK
    * checking HTML version of manual ... [105s] NOTE
    Skipping checking math rendering: package 'V8' unavailable
    Status: 1 NOTE

    [ FAIL 0 | WARN 0 | SKIP 270 | PASS 22265 ]   (tests/testthat.Rout)

22265 is the Windows CI count at 3cd2a0c8 (22264) plus the one
assertion the mo() test gained. The NOTE is this machine's known V8
note. The tarball holds the old NEWS: R CMD build copied the tree at
"preparing", before the NEWS.md edit, and `tar -xOzf` on the tarball
printed "# frmtmb 0.79.1" as its first line. (The review read the file
times and thought the build probably came after the edit; the tarball
itself was read in this session.) The new NEWS.md parses with
`tools:::.build_news_db_from_package_NEWS_md()` (246 entries, no "bad"
attribute; the development-version heading is skipped by the parser,
as in earlier rounds). No extension changed, so no extension check.
`dev/ciharden-out/` and `dev/ci-check/` were deleted after use.

## Not done, and why

- Unpinning the ubuntu-24.04 jobs: needs a green 26.04 run (above).
- The mechanism of the 26.04 ulp: not reproduced locally.
- The release checklist line: no checklist file exists (item 4).
- pkgcheck's timeout: measured and judged, not changed (item 5).
- `<<-`: superassign's.
- `frmtmb-gev-xi-chart.Rd` has no examples: informational in pkgcheck,
  outside this lane.
- The ordinal nl fits with NA headroom (review, minor 1): filed in
  `dev/test-backlog.md`, "Filed by lane ci, round of 2026-11-06", not
  fixed.
- A local reproduction on Ubuntu's amd64v3 OpenBLAS: no Windows build of
  that variant exists; building one is not cheap (item 1a).
- The gate's second red line after a check that died before its tests:
  kept on purpose (item 3).

## Other files

- `DESCRIPTION`: roxygenise with roxygen2 8.1.1 changed
  `Config/roxygen2/version` from 8.1.0 to 8.1.1; nothing else in R/ or
  man/ changed. Every lane that roxygenises makes the same change.
- `dev/test-backlog.md`: one filed item (fix round).
- `.github/workflows/brms-likelihood.yaml`: job limit 180, identity
  tier step limit 120 (fix round).
- `dev/ci-runfile.R`: prints `R.home()`, so a log shows which R ran
  (fix round; the logs above predate it).
- `dev/ci-lagrepo.R`, `dev/ci-findings-fix.py`: fix-round helpers.
- `NEWS.md`: a development-version section, three bullets. No user-
  visible change, so the work needs no version bump of its own (patch
  if released alone) and no new floor.

## Only CI can verify, after a push

- ubuntu-next on 26.04: `test-mo-coding.R` passes, the record step
  prints the CPU, "Core: ..." and the OpenBLAS package and its variant,
  the core job's gate runs, and the 24.04 Stan column finishes inside
  150 minutes from an empty cache.
- brms-tiers' new 120-minute step limit and 180-minute job limit: only
  a cold cache on a slow mirror would test them.
- A step-level timeout on "Measure coverage" fails the step without
  cancelling the job, so the gate runs.
- The hurdle_negbinomial warning is gone from both Ubuntu R-CMD-check
  jobs (it reproduced on 3 of 3 Ubuntu jobs read: 03004e4d release,
  3cd2a0c8 release and devel).
- The gate runs after a red check (needs a red check to show).
- pkgcheck's "fails continuous integration checks" clears: at 3cd2a0c8
  R-CMD-check and test-coverage are both green, so pkgcheck run
  42214681159 (still running at 13:32 UTC) should report it passing,
  and the two frmtmb-gp findings closed.
