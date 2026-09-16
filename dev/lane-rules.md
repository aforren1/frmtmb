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

## The pinned-package library, and why `.libPaths()` order matters

`C:/Users/adf44/source/r/pinlib` holds packages this project pins to a
version the shared user library does not have. It is READ-ONLY to you:
never install into it.

Today it holds one, and it is not optional if you touch Stan. `rstan`
2.32.7 declares `StanHeaders (>= 2.32.0)`, an open bound, so the user
library carries 2.39.1, which rstan cannot compile against. Any test
that compiles a FRESH Stan program dies in `compileCode()` at
`make: *** Error 1`, while anything served from `FRMTMB_STAN_CACHE`
passes, so a green suite is not evidence. `dev/machine-library.md` has
the measurement.

Order your paths so the pin wins over the user library and your own
library wins over both:

    .libPaths(c(LIB,
                "C:/Users/adf44/source/r/pinlib",
                "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

Check it rather than assume it, because the failure is quiet:

    packageVersion("StanHeaders")   # must be 2.32.10, not 2.39.1
    packageVersion("rstan")         # 2.32.7

Rebuilding the pin, if it is ever lost. The tarball is kept beside it
in `pinlib/.src`, and its source is the CRAN archive:

    https://cran.r-project.org/src/contrib/Archive/StanHeaders/StanHeaders_2.32.10.tar.gz

`R CMD INSTALL --library=<pinlib>` with Rtools 4.5 on PATH, which takes
about a minute. A dated Posit Package Manager snapshot would serve the
same purpose on a fresh machine; the archive URL is used here because
it names the exact version rather than a date that has to be looked up.

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
- **You cannot `local_mocked_bindings()` a generic frmtmb or
  frmtmb.sample binds to its owner.** 53 exported names, 25 in frmtmb
  and 28 in frmtmb.sample, are ACTIVE bindings installed by
  `frm_install_generics()`. Assigning to an active binding CALLS its
  function with the value, and that function takes no argument, so the
  mock dies with "unused argument" naming a quoted function. Mock the
  METHOD instead, `log_lik.frmtmb_draws` rather than `log_lik`, and
  know what that covers: the mock reaches the namespace binding and the
  package's own method table, never the owner's, so with the owner
  loaded only a caller inside the package's namespace dispatches to it.
  This is a rule and not a guard on purpose: a guard would have to grep
  test sources for 53 names and would fire on a correct test that mocks
  a same-named local, while this failure is loud and one run away.

- **The robust dpar accessors are public.** A density needing a
  probability near 0 or 1 uses `dpar_log()`, `dpar_log1m()`,
  `dpar_log_complement()` or `dpar_complement()`. Do not reimplement the
  log-scale gate and do not read `.eta_<dpar>` directly: it is reserved,
  not API.

## Added after the Phase 2 round, 2026-09-10 to 09-14

Four lanes ran at once. Every rule below is something a lane or a
reviewer paid for in that round, and most were found by a lane in its
own work rather than in review.

### Provenance: the number must come from the thing it describes

- **Generate counts into the document; do not type them.** A lane wrote
  "24 of 24 settled" from its target count while 9 files existed. Its
  own verifier did not catch it, and the reason generalizes: **the
  verifier checks files and this was a sentence.** The fix is
  structural, not vigilance. Emit counts from the summariser into a
  marked block and paste that block verbatim. A generated block also
  makes a PENDING figure resolve itself instead of waiting to be
  remembered. When that lane did this, the first run found its
  convergence rate was 2 of 54 by code and 5 of 54 by Hessian, not the
  "about 1 in 9" it had been carrying.
- **Count from results, never from launches.** A launcher tallied a
  replicate that ran 44 minutes and wrote nothing. Print distinct seeds
  and gaps against the grid.
- **Results live in the worktree.** Two lanes wrote replicate output to
  the session scratchpad. One was moved in time and kept 129 files
  through a session restart AND a library loss; the other was not and
  lost four arms of 60. `dev/release/` is in the repository for the
  same reason.
- **Do not edit a runner while `Rscript` is reading it.** R reads a
  script incrementally, so an edit mid-run produces a half-old,
  half-new execution. It cost one 44-minute fit and left 90 files with
  a drifted schema. Write a new file and switch the launcher instead.

### Guards

- **The guard's positive condition must be one you have OBSERVED to be
  true.** A stall guard grepped `ps -W` for a script name, but Git
  Bash's `ps -W` prints the executable path and never the arguments, so
  the "alive" branch was unreachable and the guard fired immediately.
  This is the complement of the standing rule about constructing the
  absent case, and it would have caught this in one step.
- Two guards in this round failed CLOSED on their first spelling, after
  three rounds in which every guard failed open. Both were written with
  their inverse case at the same time as the assertion, rather than
  after it. One was deliberately built to FAIL when a filed defect gets
  fixed, with a comment saying to flip it rather than delete it.

### Reporting a number

- **Observation first, fit second.** A lane put "about 1 in 2,200" in
  the plan from a lognormal fit to 40 seeds with nothing in the tail.
  Shapiro-Wilk rejected neither lognormal nor raw normal; four shapes
  that fit the same body gave 1 in 2,202, 1 in 568, 1 in 99 and 1 in
  90; and the model-free bound from 0 of 40 was about 1 in 4. Report
  what happened, then the fit, and name a disagreeing shape.
- **Say whether a relation is an IDENTITY or a measurement.** Three
  claims this round were arithmetic presented as confirmation: a
  ratio of 1.0000 with sd 0.0000 whose two sides were the same
  expression; a residual said to equal another rung's value, which it
  does exactly because paired differences subtract; and a variance
  ratio quoted over 30 replicates that was algebra plus rounding. State
  the identity, then give the residual as a numerical check at full
  precision. `%.6f` renders 1.8e-11 as 0.000000.
- **Power for a coverage study is ONE-sample.** A count is compared
  against a nominal rate known exactly. Using `power.prop.test` gave
  435 replicates where 202 were needed, a 2.1x overbuy. Power is also
  NOT monotone in n, because the rejection region advances in whole
  counts: 0.80 is first reached at n = 180 and only sustained from n =
  202. And the sawtooth is paid for in SIZE, 0.0374 at 180 against
  0.0259 at 202, so report the size beside the count.
- **An early signal at a small count is not a finding.** Two dissolved
  this round: a coverage of 3 of 6 that became 0.85 at 20, and an
  `se/sd` of 0.694 at 14 replicates that became 0.984 at 60. The second
  is the instructive one: 0.694 sat at the 0.3rd percentile of 20,000
  random 14-subsets, so it was an unlucky prefix rather than a noisy
  statistic, and the coverage count agreed with it at 14. Both
  instruments were too few. Record a dissolved signal rather than
  deleting it; it is a result about the count.
- **A shortfall in a variance component is usually ordinary.** ML
  shrinkage at q groups is about `sqrt(1 - 1/q)`. Two lanes reported an
  apparent bias that this explains: 0.4889 against 0.5 at 40 centres,
  and two eam components at 30 subjects that a single common factor
  covers. Check it before calling a component biased.
- **A withdrawn refutation must be withdrawn ON THE PAGE.** A lane
  refuted a candidate, later found its own argument wrong, and had to
  restore the candidate. A retraction left only in conversation is
  worse than never having made it, because the next reader will not
  re-open it.

### House style and Rd

- **`%` starts a comment in Rd even inside `\preformatted{}` and other
  verbatim macros.** Five unescaped signs in a new coverage table
  silently dropped the minus signs from the rendered column. Escape as
  `\%`, and **verify an Rd by RENDERING it**, with `Rd2txt` and a grep
  on the output, not by reading the source.

### The machine and the release harness

- **Do not pass `--no-manual` to `R CMD check`.** It skips the HTML and
  PDF manual sections, which is where the `%` defect above would have
  surfaced, and it makes the result incomparable to other lanes.
- **An examples-timing NOTE on this box measures LOAD.** Established
  with a control: the fixed arithmetic control swings a factor of 4.7
  on identical work, and the BASE build crossed the 5 second threshold
  at 6.75 s while the lane's own minimum was faster. The release
  harness's own timing file has the same example at 1.05 s when the
  machine was quiet and 5.28 s when it was not.
- **A stale CI expectation count is a record defect, not a breakage**,
  unless a step actually compares it. In this repo those counts sit in
  a comment block that no workflow parses. Fix them, but do not file
  them as failures.
- **Verify a restored library with a FIT, not a version string.** A
  hollow directory still answers `packageVersion()`.
  `dev/machine-library.md` carries two ten-digit reference values for
  this.
- **Remove an exemption only against a run.** A count of zero written
  call sites is not evidence that nothing reaches the exempted path,
  because dots, `do.call()` and other packages' generics dispatch
  without any call site in this repository. Run every suite that can
  reach the path, including the other extensions' suites. The 0.58.0
  dots refusal broke frmtmb.learn, which no grep of core would have
  found.


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
