# Lane wt-records: the documentation and repository records of 0.64.0

Worktree `C:/Users/adf44/source/r/frmtmb-wt-records`, branch
`wt-records`, base `ec0e6d5f` (frmtmb 0.64.0). Private library
`C:/Users/adf44/source/r/wt-records-lib`, which holds `frmtmb` and
`frmtmb.eam` built from this worktree and nothing else. Every install
passed `--library=` explicitly. No git operation was run except
read-only `git show`, `git diff` and `git status`.

Nothing here changes behavior. Two changes are code-adjacent: a
roxygen `@return` and a moved roxygen tag. Neither reaches compiled or
interpreted logic, and both were checked by rendering the Rd.

## 1. The two pkgcheck failures

### `R/me.R`, topic `frmtmb-me`: no `@return`

pkgcheck requires a `\value` section on every exported topic.
`frmtmb-me` had none.

`me()` is a formula special, not a function, so there is no return
value in the ordinary sense. Two topics of exactly the same kind
already solve this and their wording is the model followed here:
`frmtmb-multimembership` (`R/parse.R`, `mm()` and `mmc()`) and
`frmtmb-student-re` (`R/covstruct.R`, `gr(g, dist = "student")`). Both
say the term is read by `bf()` at parse time, name what the term
contributes to the model, name the accessor that reaches it, and then
say the page itself returns nothing.

The new block, at `R/me.R` just above `@name frmtmb-me`:

    @return `me()` is a formula term, not a free-standing function:
      `bf()` reads it at parse time, and the value it contributes is
      the noise-free predictor's coefficient together with the latent
      values and the measurement-model hyperparameters, reachable
      through [fixef()] and the Noise-free Terms section of
      [summary()]. This page itself documents the term grammar and
      returns nothing.

The accessors named are the ones the page's own "Parameters and
methods" section already names, so the two agree.

### `R/ad-env.R` line 85: text after `@noRd`

The `nl_rtmb_shadow` block held `@noRd`, then six lines of prose about
`atan2`, then `@noRd` again. pkgcheck refuses any text after `@noRd`.

Fixed by deleting the FIRST `@noRd` (line 85), which puts the `atan2`
paragraph back into the block it belongs to and leaves one `@noRd` at
the end. The why is kept word for word; only the tag moved. A plain `#`
comment would also have worked, but the paragraph is part of the same
argument as the `qchisq` paragraph above it, and splitting them would
have separated a reason from its subject.

### Verification

- `roxygen2::roxygenise()` with roxygen2 8.1.0 (`dev/records-roxy.R`).
  It rewrote exactly one file, `man/frmtmb-me.Rd`, and nothing else in
  the tree, confirmed by `git status`. No Rd appeared or disappeared
  for `nl_rtmb_shadow`, so the second `@noRd` still hides it.
- `tools::Rd2txt("man/frmtmb-me.Rd")` (`dev/records-rd.R`, rendered
  output kept as `dev/records-me-rd.txt`). The Value section renders at
  line 24 of the output, with the full paragraph under it. Read the
  rendered text rather than the Rd source, because `Rd2txt` writes the
  heading as `_V_a_l_u_e:` with backspace overstrike, which a grep for
  `^Value:` misses.
- All 105 Rd files in `man/` parse with `tools::parse_Rd()`, 0
  failures.
- One Rd has no `\value`: `man/frmtmb-shared-generics.Rd`. It carries
  `\keyword{internal}`, which is the exemption pkgcheck honors, so it
  is not a second instance of this defect and was left alone.

### The scan for other `@noRd` blocks followed by text

`dev/records-nord-scan.R` scans `R/` and every `extensions/*/R`. It
reports a roxygen line with non-blank, non-tag content anywhere inside
a block that has already reached `@noRd`, which is a wider net than
"the line immediately after", because pkgcheck's rule is about the
block and not the line.

- On this worktree: 120 files scanned, 0 hits.
- On `R/ad-env.R` as `HEAD` has it: 6 hits, lines 86 to 91. The
  instrument was seen firing on the defect it exists to find before it
  was believed on the clean tree.

Nothing in the test suite asserts this rule. It is stated in
`dev/lane-rules.md` and enforced only by pkgcheck on CI, which is why
the defect reached a release. The scanner is left in `dev/` so that a
later lane can promote it to a test; that was not done here because the
promotion needs its own failing construction and a decision about
running it inside `R CMD check`, where an installed package has no
sources to scan.

### Test files rerun

Installed into the private library with
`R CMD INSTALL --library=C:/Users/adf44/source/r/wt-records-lib`, then
one test file per R process (`dev/records-run1.R`, `NOT_CRAN=true`):

| file | pass | fail | err | skip |
|---|---|---|---|---|
| `tests/testthat/test-me.R` | 70 | 0 | 0 | 0 |
| `tests/testthat/test-methods-audit.R` | 58 | 0 | 0 | 0 |
| `tests/testthat/test-message-uniqueness.R` | 6 | 0 | 0 | 0 |
| `tests/testthat/test-bracket-access.R` | 33 | 0 | 0 | 0 |

The last three were run because they read `R/` sources or the
documented method surface, which is what these two edits touch.

## 2. `dev/test-backlog.md`

Three edits.

**The five wt-mvprior punch-round-1 items are closed.** The heading
`### Open - medium` under "Filed by wt-mvprior after punch round 1,
2026-09-24" became `### Closed at 0.64.0`. Each of the five entries is
kept verbatim and gains a `SHIPPED:` line naming what closed it and the
lane that did it, from NEWS.md's 0.64.0 section: `bf + bf + bf` (one
`+` method, lane mv), `cumulative()` in a multivariate model (lanes mv
and thres), `me()` (lane me, with a pointer to `?frmtmb-me`),
`0 + Intercept` (lane icpt0) and `student()` with `set_rescor(TRUE)`
(lane mv). The section's lead sentence "Not fixed; outside that lane"
became "Not fixed by that lane", which is what it now means.

**A new section, "Filed at the 0.64.0 release (2026-09-25)."** The ten
"Left open" bullets of `dev/parity-round-20260925.md`, split into
`Open - high priority` (the four that are wrong answers or lost
information), `Open - medium`, `Upstream`, and a closing subsection for
the two older record inconsistencies this lane closed. Each bullet
names its lane and its findings file. Four carry "In progress:"
with the lane now working on them: wt-gradcheck (the convergence
check), wt-thresrefit (refits recounting thresholds, and
`draw_prior_entry()`), wt-csfactor (`y ~ x + cs(x)` unidentified),
wt-arcovsample (frmtmb.sample `log_lik()`/`loo()` for `cov = FALSE`
ARMA).

**Two existing items marked in progress.** The wt-predfix `cs()`
integer-codes entry now opens "IN PROGRESS: wt-csfactor" instead of
"OPEN", keeping "A SILENT WRONG ANSWER, pre-existing". The
`re_formula = NA` smooth entry now opens "IN PROGRESS: wt-resmooth"
instead of "DECIDED, waiting for a lane", keeping the user's decision
and its date. Neither entry's measurements were touched.

## 3. `dev/round-handoff.md`

Rewritten above "Decisions the user made on 2026-09-24" and below
"What the user has settled", with the middle kept.

- Header dated "rewritten 2026-09-28 at the 0.64.0 brms-parity
  release", with a pointer to `dev/parity-round-20260925.md`.
- The nine lanes, each with its findings file, from that record's
  table.
- Versions read from the DESCRIPTIONs rather than retyped: frmtmb
  0.64.0; frmtmb.sample 0.12.0 floored on frmtmb 0.64.0; frmtmb.eam
  0.11.0, frmtmb.latent 0.6.0 and frmtmb.learn 0.7.0 still floored on
  0.63.0; frmtmb.coupling 0.6.0, frmtmb.spline 0.8.0 and frmtmb.ode
  0.7.0 still floored on 0.61.0.
- The verification section says plainly that it ran on a Linux
  container, that `R CMD check` was `--no-manual` there because the
  container has no LaTeX, and that a Windows run for the manual
  sections is still owed. It carries the suite, gated, ported-brms,
  build and check numbers from the round's record, and the two files
  that fail identically on the base build there.
- A paragraph on the fresh Windows R install of 2026-09-28 and the
  restore from `dev/restore-cran-2026-09-28.txt`, naming what this lane
  corrected because of it.
- The merge defects are summarized in one paragraph, with the table
  left in the round's record instead of duplicated.
- "What is next" now leads with the six lanes running, then me()
  hyperparameter priors, `gr(g, by = f, cov = A)`, latent-residual AR
  for non-gaussian families, the multivariate family-fill decision, and
  Phases 4 and 5 of `dev/extension-gaps-plan.md`. The backlog and
  upstream paragraphs were updated for what shipped: three of the
  drmTMB items are gone because 0.64.0 shipped them, and the brms
  `posterior_predict_hurdle_negbinomial()` defect is added to the
  unfiled upstream list.
- "Worktrees" now names the six live worktrees off `ec0e6d5f`.
- "What this round is evidence for" was retitled "What the 0.63.0
  round was evidence for", with its four lessons kept verbatim in
  substance, and a short 0.64.0 section added above it: a lane cannot
  test the seam it shares with another lane, and verification on one
  platform leaves a hole.

Verbatim checks, by diff against `HEAD`:

- "Decisions the user made on 2026-09-24", including "Settled earlier,
  do not reopen", is byte-identical.
- "How a round runs here" and "What the user has settled" differ in
  exactly three places, all of them deliberate: "this release" became
  "the 0.63.0 release" in the `R CMD check` paragraph, and the memory
  and StanHeaders items were rewritten (below).

## 4. The StanHeaders pin, which no longer exists

`dev/release/run-tests.R` lines 16 to 32 and
`dev/tmbstan121-findings.md` both record the retirement on 2026-09-17.
Measured on this machine today: `C:/Users/adf44/source/r/pinlib` does
not exist at all, and `C:/Users/adf44/Documents/.R/Makevars.win` holds

    CXX17FLAGS += -mtune=native -O3 -mmmx ... -msse4.2
    CXX17FLAGS += -std=gnu++17

The brief said `pinlib` was an empty directory; on this machine the
directory is gone. The rewritten text says it holds nothing, which is
true either way.

The two assertions `dev/release/run-tests.R` makes before it runs a
test, now named in the rules:

    packageVersion("tmbstan") >= "1.2.1"
    any(grepl("-std=gnu++17", readLines(tools::makevars_user())))

**`dev/lane-rules.md`.** The section "The pinned-package library, and
why `.libPaths()` order matters" became "Stan on this machine, and what
replaced the pinned library". It names the Makevars file and the line,
says why the line is needed (StanHeaders 2.39 needs C++17 and R does
not select it), gives the two assertions, says that a suite served from
`FRMTMB_STAN_CACHE` is weaker evidence than one fresh compile, keeps
the history in one sentence with a pointer to
`dev/tmbstan121-findings.md`, and replaces the three-entry
`.libPaths()` example with the two-entry one. The rebuild recipe and
the CRAN archive URL are gone, because there is nothing to rebuild.

**`dev/machine-library.md`.** Three places:

- "Restoring it": the StanHeaders step is removed from the restore
  order, and a paragraph says what replaced it. It makes the point that
  the Makevars file is outside the library, so a library restore does
  not write it and nothing reports its absence until a fresh Stan
  compile fails.
- "The StanHeaders trap, which a restore walks straight into": the
  2026-09-09 measurement is untouched. Its closing instruction, "after
  any restore, pin StanHeaders", is replaced by a paragraph headed
  SUPERSEDED ON 2026-09-17 with the current remedy and what to check
  instead. The trap itself is still described, because the open version
  bound that causes it has not changed.
- "What was decided, 2026-09-09": the pin paragraph is in the past
  tense and ends with the date the arrangement ended.
- The lesson "The pin lives outside `%LOCALAPPDATA%`", in the
  five-losses section, keeps its measurement and gains one sentence:
  the pin went, and the Makevars file that replaced it has the same
  property.

The remaining `pinlib` mentions in `dev/machine-library.md` are inside
dated tables of what survived each library loss. They are history and
were left alone.

**`dev/organizer-rules.md` does not mention the pin, StanHeaders or
`.libPaths()` at all**, so there was nothing to change there. Checked
by grep for `pin`, `StanHeaders`, `2.32` and `libPaths`; the only hits
are "A test that pins a bug" and "false-alarm".

Twenty other files under `dev/` name `pinlib`: lane findings files,
review files and lane scratch scripts, all dated records of runs that
really did use the pin. Rewriting them would falsify the record, so
they were left.

## 5. Memory, after the machine change

Instructed mid-lane by the user: this machine has 63 GB of RAM and a
fast disk, so the "at most 3 fitting processes, 5 GB free" cap is
lifted. One test file per R process stays.

- `dev/lane-rules.md`, the memory bullet under shell traps: the
  standing rule is now "on the machine of 2026-09-28 test runs need no
  process cap; keep one test file per R process". The 2026-09-24 crash,
  the 11 processes at up to 2.4 GB, the 34 `std::bad_alloc` failures
  and the clustered NaN gradients are kept as history, as is the rule
  that a `bad_alloc` or NaN gradient under memory pressure is not
  evidence about the model until it reproduces alone.
- `dev/round-handoff.md`, "How a round runs here": the same
  replacement, with the crash kept in one sentence and a pointer to
  `dev/lane-rules.md`.

This lane ran three test files at once to use it, with no trouble.

## 6. `codemeta.json`

`codemetar` 0.3.7 is in the user library, so the file was regenerated
with `codemetar::write_codemeta()` (`dev/records-codemeta.R`) rather
than hand-edited. The previous file was kept for the diff.

The diff is 97 insertions and 22 deletions, and every part of it
matches the root DESCRIPTION:

- `version` 0.50.0 to 0.64.0, and the citation's
  "R package version 0.47.0" to 0.64.0. The citation string was two
  releases staler than the version field.
- Added to `softwareSuggestions`: `ordinal` (the one the brief asked
  for), and `brokenstick`, `drmTMB (>= 0.7.0)` and `ggplot2`, which
  are all in Suggests today and were added after the file was last
  written.
- Added to `softwareRequirements`: `generics` and `nlme`. `nlme` moved
  from Suggests to Imports in the meantime, and the regeneration
  follows the move, which renumbered the requirement keys 2 to 14.
- Version bounds that DESCRIPTION gained: `posterior (>= 1.0.0)` and
  `tmbstan (>= 1.2.1)`.
- `fileSize` 2779.599KB to 6511.271KB, which is the tree, not the
  package.

Reviewed and judged wrong: nothing. Two things are worth knowing
rather than fixing. `codemetar` reaches the network (CRAN metadata,
Bioconductor metadata and the GitHub API) and writes a
`"sameAs": "https://CRAN.R-project.org/package=<pkg>"` for every
non-base dependency it resolved, including `drmTMB` and `RTMBdist`;
that field was already there for `RTMBdist` before this round, so the
shape is unchanged, but it is a claim about CRAN that this lane did not
independently check. And the file carries no `dateModified`, so the
regeneration does not record when it ran; `codemetar` did not add one.

`codemetar` printed "Added codemeta.json to .Rbuildignore".
`.Rbuildignore` already had `^codemeta\.json$` at line 16 and is
unchanged, confirmed by `git diff`.

## 7. `SPEC.md`

Two minimal edits, no restructuring.

- The Status paragraph's first sentence, "design specification; v0.19
  implemented ... The original roadmap and its extensions are
  complete", became "the roadmap of this specification is complete, and
  the package has gone past it. NEWS.md is the changelog, through
  0.64.0. What follows is the record of how the roadmap was built,
  milestone by milestone". The forty lines of milestone history after
  it are untouched, and now read as the record they are.
- The Deferred line named six things, five of which had shipped. It now
  defers `CAR/SAR` alone, and a following sentence names what shipped
  and when: `me()` at 0.64.0 with a pointer to `?frmtmb-me`, `mo()` at
  v0.15, `cs()` and `mixture()` at v0.17, `gp()` at v0.18, and `mi()`,
  which the line called Excluded, at v0.16 and v0.17. Every date comes
  from the Status paragraph in the same file.

Only `me()` was in the brief. The other five were changed in the same
sentence because leaving them would have left the paragraph false in
five more places, and separating them would have cost more words than
the fix.

## 8. frmtmb.eam's Suggests and its CI job

**`test-sampling.R` already skips without frmtmb.sample.** It calls
`skip_if_not_installed("frmtmb.sample")` at line 13, before the
`frmtmb.sample:::tmbstan_build_broken()` call that would otherwise
error, with a comment explaining why. So the missing Suggests never
produced a failure; it produced a silent skip and an undeclared
dependency.

Measured, not assumed. `frmtmb.sample` is in NEITHER the user library
nor this lane's private library, so the skip fires here:

    frmtmb.eam tests/testthat/test-sampling.R
    pass=94 fail=0 err=0 skip=1

**What was added.** `frmtmb.sample` to
`extensions/frmtmb.eam/DESCRIPTION` Suggests, in alphabetical position
after `brms`.

Nothing else was added. Every `::` and `:::` call in
`extensions/frmtmb.eam/tests/` and `vignettes/` was harvested:
`EMC2`, `RTMB`, `RWiener`, `Rmpfr`, `WienR`, `base`, `brms`, `frmtmb`,
`frmtmb.eam`, `frmtmb.learn`, `frmtmb.sample`, `knitr`, `rmarkdown`,
`rtdists`, `statmod`, `stats`, `testthat`, `tinyplot`, `tools`,
`utils`. Of the three not declared, `EMC2` and `frmtmb.learn` appear
only in prose and comments: `1 - EMC2:::pWald()` is inside a sentence
of `vignettes/ddm.Rmd`, not a chunk, and `frmtmb.learn::rlddm()` is
named in three test-file comments explaining why a seam exists. By
this repository's Suggests-follow-use rule neither belongs in
Suggests.

**The workflow.** `.github/workflows/check-frmtmb-eam.yaml` now:

- lists `'extensions/frmtmb.sample/**'` in both `paths:` filters;
- installs frmtmb.sample from the checkout, after core and before the
  step that reads this package's Suggests, so that step finds it in
  `installed.packages()` and does not ask a repository for a package
  that is in no repository;
- says in its header that two packages come from the checkout, which
  is the shape `check-frmtmb-learn.yaml` established, and drops the
  claim that every suggested dependency is on CRAN.

frmtmb.sample imports only frmtmb and base packages, and core's own
dependency resolution already brings RTMB, so the install needs no
extra step of its own. Its Suggests are not installed, because nothing
in this job runs its tests. That is the same choice
`check-frmtmb-learn.yaml` makes for frmtmb.eam.

**One stale record fixed while there.** Both `paths:` filters listed
`.github/workflows/check-frmtmb-ddm.yaml`, a file that does not exist:
the workflow was renamed when the package was renamed from frmtmb.ddm
and the self-reference was not. So an edit to this workflow did not
re-run its own job. Corrected to `check-frmtmb-eam.yaml`. The job id
`check-frmtmb-ddm` is left alone, because
`tests/testthat/test-ci-siblings.R` finds the file by its `name:` line
and the id is only a label.

**Evidence that the test asserts what it is supposed to.**
`tests/testthat/test-ci-siblings.R`, one file per R process, at three
points:

| tree | pass | fail |
|---|---|---|
| base, before the DESCRIPTION change | 11 | 0 |
| DESCRIPTION changed, workflow not yet | 12 | 0 with 2 failures |
| DESCRIPTION and workflow both changed | 14 | 0 |

The middle row is the one that matters: the test was SEEN FAILING on
exactly the state the brief describes, with both of its assertions
firing, the install one and the `paths:` one. The count rises by three
because the file's `expect_gt(checked, 0)` and both sibling assertions
are new work: the loop had one sibling pair before (frmtmb.learn on
frmtmb.eam) and has two now.

## What was found and NOT changed

- **`man/frmtmb-shared-generics.Rd` has no `\value`.** It is
  `\keyword{internal}`, which pkgcheck exempts, and the pkgcheck run
  named only `frmtmb-me`. Left alone.
- **No test asserts the `@noRd` rule.** `dev/records-nord-scan.R` is
  the instrument; promoting it to a test needs its own decision about
  running it under `R CMD check`, where an installed package has no
  sources. Filed here rather than done.
- **`codemeta.json` claims CRAN for `drmTMB`.** `codemetar` resolved
  it against the CRAN metadata it downloaded. Not independently
  checked by this lane.
- **Twenty findings and review files under `dev/` still describe the
  pin as live.** They are dated records of runs that used it.
- **`dev/suite-baseline.tsv` was not touched**, as instructed. The
  consolidating session regenerates it from the suite now running.
- **No version number was changed anywhere**, per
  `dev/organizer-rules.md`. See the report for which bump the eam
  change needs.

## Scripts and outputs in this lane

| file | what it does |
|---|---|
| `dev/records-roxy.R` | roxygenises the worktree from the private lib |
| `dev/records-rd.R` | renders the me Rd, greps Value, parses all Rd |
| `dev/records-me-rd.txt` | the rendered Rd, kept as the evidence |
| `dev/records-nord-scan.R` | the noRd prose scanner, exit 1 on a hit |
| `dev/records-run1.R` | one test file, one R process, from the private lib |
| `dev/records-codemeta.R` | regenerates `codemeta.json` |
