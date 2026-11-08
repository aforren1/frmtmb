# How to check a release before you commit it

Use this list after you merge the lanes of a round into the release tree and
install the release library. Do each step on the final build. Read the `lib:`
line of each log, and make sure that it names the release library.

## Before you start

- Use absolute library paths. Install only into the new release library.
- Run at most three R processes at a time. Run one test file per process.
- Do not stop a process that you did not start.
- Put Rtools at the end of `PATH` under Git Bash.

## Steps

1. Compare the library with CRAN:
   `Rscript dev/release/cran-lag.R "<release library>;<user library>"` (exit
   0, or each lagging package updated or explained).
2. Point `dev/release/run-tests.R` and `dev/release/run-session.R` at the new
   release library (`LIB`).
3. Set `Config/roxygen2/version` in every `DESCRIPTION` to the roxygen2
   version that the user library holds.
4. Run the test suites of core and of every extension, one file per process,
   with `ROUND_BASE=<previous release library> bash dev/ciharden-scan.sh <run>
   ref <release library> 3`:
   - the plain tier (`NOT_CRAN=true` only, `CIHARDEN_UNGATED=true`);
   - the gated tier (every gate on: brms, drmTMB, fuzz, scale);
   - frmtmb.sample with its gates on.
   Each count must show 0 failures, 0 errors and 0 escaped warnings.
5. Run the new and changed test files on the previous release library. Record
   how many fail there.
6. Run each suite in one session, in three orders (forward, reverse,
   preload), with `dev/release/run-sessions.sh`. Each `SESSION` line must show
   `fail=0 err=0 warn=0`, and a `LEAK` line must have a tag.
7. Run the numeric test files under OpenBLAS 0.3.26 and 0.3.32 with LAPACK
   (`dev/ciharden-openblas.sh <version> lapack`, then `dev/ciharden-scan.sh`
   with the `ob<version>-lapack-t4` configuration). Delete
   `dev/ciharden-out/` when you are done.
8. Run `R CMD check --as-cran` on core and on each package whose code
   changed, with `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on
   `PATH`.
9. Regenerate the documents with `dev/release/build-docs.R`. Make sure that
   `devtools::document()` changes no file and that `pkgdown::check_pkgdown()`
   finds no problem.
10. Make sure that the installed code equals the tree, with
    `dev/rel080-installed.R` (copy it and change the library): count the
    functions that differ (0 expected). Compare the top-level constants too,
    not only the functions.
11. Rerun the key table of each lane on the merged build with its own script,
    and compare each number with the lane's and the review's numbers.
12. Run `dev/equiv-dump-compat.R`, then the equivalence survey on all pairs
    with `RECORD=1`:
    `FRMTMB_EQUIV_LIB="<release lib>;<previous release lib>" RECORD=1 bash dev/equiv-run.sh <N> <outdir>`.
    Give the headline beside the previous one, and explain each verdict that
    changed.
13. Regenerate the ported brms-suite ledger once, with `dev/rel080-ledger.sh`.
    Report each row that flips.
14. Check NEWS:
    - one new section per package, with its version number;
    - no lane names and no `dev/` paths;
    - do not edit an earlier section;
    - each behavior change is marked BREAKING;
    - no em dashes, no emojis, US English.
15. Check the versions and floors. Raise an extension's `frmtmb (>= )` floor
    when it calls a function that the release adds or changes. Test the floor:
    the extension's tests on the previous core must fail where the floor says
    they do.
16. Merge the lanes' items into one section of `dev/test-backlog.md`, and
    file upstream bugs only in `dev/upstream-bugs.md`. Post nothing outside.
17. Write the release record (`dev/round-<date>.md`), update
    `dev/round-handoff.md`, and write the commit script. Dry-run the script.
    Do not run it without the user's approval.
