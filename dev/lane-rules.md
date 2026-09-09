# Standing rules for a lane worker or reviewer

Give this file's path to every agent that writes or reviews code in a
worktree. It was reconstructed into the repository after the session
scratchpad holding it was cleaned mid-round, which is why it lives here
now: a process document that agents depend on cannot live in a
session-local temporary directory.

## Where you work

Your worktree is the ONLY tree you touch. Never edit
`C:\Users\adf44\source\r\frmtmb`, which is main, and never edit another
lane's worktree. Your branch is checked out for you.

## Git

DO NOT COMMIT. Do not stage, branch, merge, push or stash. Leave every
change in the working tree. The user has not authorized git operations
from lane agents; consolidation is done by the session that owns the
round.

## The R library rule, and why it is absolute

This machine has had its shared R library destroyed twice by concurrent
installs. R installs by deleting a package directory's contents and then
unpacking, so an install into a library another R process is reading
leaves hollow directories.

- Create ONE private library and install ONLY into it:
  `install.packages(..., lib = LIB)` and `R CMD INSTALL --library=LIB`.
  Always pass `lib=` explicitly. NEVER `.libPaths()[1]`, never the
  default.
- The user library `C:/Users/adf44/AppData/Local/R/win-library/4.6` is
  READ-ONLY to you. Put it LAST in `.libPaths()` as a fallback. Never
  install into it, for any reason, including a missing dependency.
- Never use `dependencies = TRUE`.
- Other lanes may be running. Your private library is yours alone, so
  your installs are safe; installing anywhere else is not.

## Toolchain

- R 4.6.1 at `C:\Program Files\R\R-4.6.1\bin`.
- `.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))`
  at the top of every script.
- Tests: `NOT_CRAN=true`, ONE TEST FILE PER R PROCESS. A whole-suite run
  in one process has repeatedly hidden state leakage.
- Gated tiers need `FRMTMB_BRMS_FIT_TESTS=true` and `NOT_CRAN=true`.
  `FRMTMB_STAN_CACHE` defaults to `dev/stan-cache` in the tree; leave it
  there or every model recompiles. That directory is GITIGNORED, so a
  fresh worktree has none: copy it from the main checkout before running
  anything gated.
- `R CMD check --as-cran` needs pandoc and TinyTeX on PATH or it reports
  a bogus `pdflatex` ERROR:

      $env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
      $env:PATH = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools;C:\Users\adf44\AppData\Roaming\TinyTeX\bin\windows;" + $env:PATH

  Only the V8 math-rendering NOTE on the HTML manual is expected. A
  CRAN-incoming WARNING about frmtmb not being on CRAN is also
  pre-existing; `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE` removes it.

## Shell traps that have cost this project real time

- Bash heredocs here EAT ONE BACKSLASH even when quoted, which breaks R
  regexes. Write R scripts with the Write tool, or use bracket classes
  (`[(]` not an escaped paren). This bit a lane again in the 0.55.1
  round, caught by the run rather than by review.
- PowerShell `Set-Content` writes a UTF-8 BOM that R chokes on. Use
  `[System.IO.File]::WriteAllText(path, text, (New-Object System.Text.UTF8Encoding($false)))`.
- PowerShell flattens `@(@("a","b"))` into `@("a","b")`. Build job lists
  with `New-Object System.Collections.ArrayList` and `[void]$x.Add(...)`.
- `$_.Name` inside `ForEach-Object` piped to an external command does
  not reach the child. Use `foreach ($f in $files) { $nm = $f.Name ... }`.
- Prefix every scratch file and log with your lane name.

## House style, which the reviewer will check

- US English. No em dashes, no emojis, and no spaced hyphens standing in
  for an em dash: restructure instead.
- 80 columns in R sources and in markdown you write.
- Documentation follows ASD-STE100 and the Google developer style guide;
  structure follows diataxis.
- Comments say WHY, not how.
- A `@noRd` block must not be followed by any text; pkgcheck fails on
  it. Never let two roxygen blocks run together.
- NEVER write an absolute numeric tolerance in a test. Express it as a
  ratio to something the run itself measures. Four have broken on CI or
  in review: `test-deriv.R`, `test-gratia.R`, the fuzz harness's
  permutation check, and a spectral assertion that passed only on a
  pinned seed.

## The evidence standard

- "No test reaches it" is NOT evidence that code is unreachable. Prove
  reachability or unreachability by construction.
- A test that pins a bug is worthless unless you have SEEN IT FAIL.
  Run it against the unfixed code and record the failure. An error
  because a symbol does not exist is a weak form of this; construct the
  behavioural failure.
- A check that fires on a correct model is worse than no check. Measure
  the false-alarm rate, on designs the field actually produces.
- Report numbers, not adjectives: what you measured, on what, with how
  many replicates.

## Added after the 0.55.0 and 0.55.1 rounds, where each was found only in review

- **Check the instrument before believing a measurement.** Two timing
  claims evaporated. `proc.time()` ticks at 10.0 ms here, so 200 calls
  at 200 us is four ticks and a "0.75x speedup" was four ticks against
  three. Another measured a call that runs ONCE PER FIT. Interleave the
  arms in one process, grow blocks past 1.2 s, take a minimum of several
  rounds, and carry a CONTROL built from the same code that must report
  1.0. Better still, count something load-independent: AD nodes settled
  one of these where the clock could not.
- **A count is not a count if it was capped, stale or skipped.** Three
  lanes reported one. testthat's summary reporter stops at ten failures
  and says so. A runner that sums `failed` and not `error` prints a
  clean line for a file that aborted halfway, and one did: PASS=51
  FAIL=0 on a file that asserts 59.
- **A number that will not reproduce means you have not found the
  construction.** Two lanes called a recorded figure false and both were
  wrong. Search `dev/` and `dev/reviews/` before concluding a record is
  wrong, and say which construction you used when you record a number.
- **A printed zero is not a measured zero.** `%.6f` renders 1.8e-11 as
  0.000000. Prefer `identical()` or a ulp count when you mean bitwise.
- **Construct the case where the guarded thing is ABSENT.** Every guard
  built in the 0.55.1 round failed open on its first try: a loop that
  hit `next` on every name and asserted nothing, a hash check that
  passed on a deleted file, an argument whose affirmative value silently
  disabled the guard it gated.
- **New R files need the hazard guard.** A new package or `R/` directory
  needs `frm_hazard_reads()` in its own `test-bracket-access.R`. A `$`
  read on a hazard container partial-matches to a neighbouring slot.
- **A sibling dependency needs its CI job to install it.** An extension
  that depends on another extension must have that sibling installed
  from the checkout in its workflow and listed in `paths:`.
  `tests/testthat/test-ci-siblings.R` asserts both.
- **The robust dpar accessors are public.** A density needing a
  probability near 0 or 1 uses `dpar_log()`, `dpar_log1m()`,
  `dpar_log_complement()` or `dpar_complement()`. Do not reimplement the
  log-scale gate and do not read `.eta_<dpar>` directly: it is reserved,
  not API.

## Cost, which is a real constraint

- **Use the round's shared reference library. Do not build one.** The
  organizing session installs ONE copy of the base commit per round and
  gives you its path. It is READ-ONLY: never install into it, and put it
  after your own private library in `.libPaths()`. Your private library
  holds your change; the shared one holds the base you measure against.
- **`R CMD check --as-cran` runs ONCE, on your final pass.** It is ten
  to twenty-five minutes, and the consolidating session runs the
  authoritative one anyway. Running it after every edit is the single
  largest avoidable cost in a lane.
- **While you iterate, run the test files your change can reach.** The
  whole package suite runs once, on that same final pass.
- **Record the seed and the script path beside every number.** A number
  without its construction gets re-derived from scratch by whoever reads
  it next, and that has twice produced a confident and wrong claim that
  a shipped figure was false.

## What you deliver

1. The change in your worktree, roxygenised, installed into your private
   library, with the affected test files rerun one per process.
2. `dev/<lane>-findings.md`: what you changed, what you measured, the
   numbers, what you decided NOT to do and why, and any defect you found
   but did not fix. An honest "I could not settle this" with the
   measurement attached is worth more than a guess.
3. A final report naming every file you touched and every test file you
   ran with its counts.

If part of your task turns out to be wrong, say so with the measurement
rather than doing it anyway. If part is blocked, finish the rest in full
and say exactly what you left and why.
