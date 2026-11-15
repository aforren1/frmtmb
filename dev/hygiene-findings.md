# Lane hygiene: findings (round of 2026-11-13)

Base ec8178d3 (frmtmb 0.82.0), worktree `frmtmb-wt-hygiene`, branch
`wt-r83-hygiene`, library `wt-hygiene-lib`. Brief:
`round-briefs/r083-hygiene.md`. Line numbers are 0.82.0 lines of
`dev/test-backlog.md`.

Tools this lane wrote:

- `dev/hygiene-run-tests.R`: one test file per process. `HYG_ARM=mine`
  puts the lane library first; `HYG_ARM=base` is rellib-r21 alone.
  `HYG_BARE=true` attaches the package and calls `test_file(f)` with no
  `package =` and no `env =`, which is how a reader runs one file by hand.
- `dev/hygiene-par.sh`: `dev/release/run-par.sh` over that runner, logs
  in `dev/hygiene-out/<tier>.log`.
- `dev/hygiene-barenames.R`: lists the names a test file reads that an
  attached session cannot see but the package namespace (or its imports)
  can.

## Item state

| line | item | state |
|---|---|---|
| 131 | constant-weight gaussian | closed: test-weights-constant.R, vignette sentence |
| 764 | ode mocks without `.package` | closed |
| 769 | three files need load_all() | closed: qualified, pass bare |
| 1166 | wiener coverage | open: about 17 process-hours |
| 1222 | REML x mi() vs me() | closed: two compat rows conditional, measured |
| 1252 | codemeta.json | closed: regenerated, test-codemeta.R, checklist step 3 |
| 1519 | rate() x cens()/trunc() | closed: agrees with brms to 5.1e-13; rows works |
| 2733 | fixture 1 theta under OpenBLAS | open: selimit2's area |
| 3018 | gr(g, by = f) fixture | open: selimit2's area |
| 3042 | R7 staleness via setup | closed: helper fix, test-port-harness.R |
| 3045 | R8 `lcdf` name | documented; rename is a user question |
| 3154 | reload duplicates registrations | closed: register_once() merges only the provably same (fix round 2) |
| 3164 | extension session state | closed: guarded in run-session.R |
| 3171 | registries helper hooks, proof by count | closed |
| 3177 | coverage path filter | closed: plan job, carryforward |
| 3309 | bin-2 lncor, category order | open: M, ledger machinery |
| 3631 | translator, per-category groups | closed: helper fix, new gated block |
| 3649 | posterior_epred, stancode files | open: bin 3 by the audit |
| 4042 | guardrail1 first home only | closed |
| 4259 | bin-2 survey gap note | closed |
| 4448 | bin-2 adapters `../../` | closed |
| 4628 | survey sample redraw | open: consolidation's, equiv-data |
| 5516 | conditional refused half count | closed |
| 5533 | translator gaps | open: M |

## 764 and 769: test files that ran only under the namespace

Before (rellib-r21, `HYG_ARM=base HYG_BARE=true`, log
`dev/hygiene-out/bare-base.log`, original test files):

    frmtmb RESULT test-compat-register.R pass=77 fail=0 err=7 skip=0 warn=0
    frmtmb RESULT test-structure.R pass=63 fail=5 err=16 skip=0 warn=0
    frmtmb.sample RESULT test-loo.R pass=119 fail=0 err=9 skip=2 warn=0
    frmtmb.ode RESULT test-ode.R pass=70 fail=0 err=1 skip=0 warn=0

The ode error is the item-764 one: `local_mocked_bindings()` without
`.package` says "No packages loaded with pkgload". Two mocks in
test-structure.R had the same defect (not named in the item).
`dev/hygiene-barenames.R` listed 7 hidden names in
test-compat-register.R, 10 in test-structure.R and 20 in test-loo.R (of
which `arma`, `unstr` appear only inside formulas and `brmsfamily` only
as `brms::brmsfamily`, so they need nothing). Each other name is now
qualified (`frmtmb:::`, `frmtmb::` for exports, `frmtmb.sample:::`), and
each mock names its package. Assignments into an internal environment go
through a local alias (`contrib <- frmtmb:::frmtmb_compat_contrib`),
because `frmtmb:::x$y <- v` is not valid R.

After (same arm, edited files; `bare-after.log` and the normal runner,
`norm-after.log`):

    frmtmb RESULT test-compat-register.R pass=103 fail=0 err=0 skip=0 warn=0
    frmtmb RESULT test-structure.R pass=128 fail=0 err=0 skip=0 warn=0
    frmtmb.sample RESULT test-loo.R pass=138 fail=0 err=0 skip=2 warn=0
    frmtmb.ode RESULT test-ode.R pass=72 fail=0 err=0 skip=0 warn=0

The bare and the normal runner give the same counts. The 2 skips are the
brms-gated tests (`FRMTMB_BRMS_FIT_TESTS` unset).

## 1252: codemeta.json

Before: `codemeta.json` said `"version": "0.79.1"` with DESCRIPTION at
0.82.0 (its last regeneration was by hand at 01a1f3da, which edited the
version line only; the citation entry still said "R package version
0.65.0").

Change: regenerated with `codemetar::write_codemeta(pkg = ".")`
(codemetar 0.3.7). Three lines move: `version` 0.82.0, the citation's
"R package version 0.82.0", and `fileSize` 6719.856KB to 12207.244KB.
`fileSize` is codemetar's sum over the files `.Rbuildignore` keeps
(`dev/hygiene-out/cm-size.R`: R/ 5108.9 KB, tests/ 4902.4 KB, man/
912.6 KB, NEWS.md 765.4 KB; 593 files). The package grew; it is not
worktree junk.

Where the release flow keeps it current: step 3 of
`dev/release/checklist.md` now says to regenerate it after the version
bump, and the new `tests/testthat/test-codemeta.R` fails in step 4 (the
one-file-per-process run from the source tree) until the file carries the
DESCRIPTION version and the same dependency set. The test reads the source
tree and skips under R CMD check and covr, which hold no codemeta.json
(it is in `.Rbuildignore`).

Seen to fail (rellib-r21, the 0.82.0 codemeta.json restored with `git
show ec8178d3:codemeta.json`; `dev/hygiene-out/cm-before.log`):

    RESULT test-codemeta.R pass=3 fail=1 err=0 skip=0 warn=0
    DETAIL FAIL [codemeta.json carries the version of DESCRIPTION] ...
      `actual`: "0.79.1" | `expected`: "0.82.0"

The dependency arm, with `ordinal` renamed to `ordinalX` in a scratch
copy (`cm-dep.log`): pass=3 fail=1, at `[42] "ordinalX"`. With the
regenerated file (`cm-after.log`): pass=4 fail=0.

## 4042: test-guardrail1.R checked the first extension special only

Before (rellib-r21, `g1-orig-plain.log` and `g1-orig-pre.log`, the
latter with frmtmb.spatial, frmtmb.eam, frmtmb.survival loaded first):

    plain    RESULT test-guardrail1.R pass=37 fail=0 err=0 skip=0 warn=0
    preload  RESULT test-guardrail1.R pass=36 fail=0 err=0 skip=0 warn=0

That is the 37 against 36 of the item. The block also read `cs_names()`
by bare name and missed the autocor specials (`sar` is registered through
`frmtmb_register_autocor()`, not as a covstruct). It now checks every
home in `extension_special_homes` (car, sar, spde): a column named like
the special is data, and a call is refused by name when nothing
registered it, or left to its extension when something did. One
assertion per branch, so the count is the same whatever loaded first.
After (`g1-new-plain.log`, `g1-new-pre.log`): pass=42 on both.

## 3171: local_frmtmb_registries() dropped hooks added in scope

`frmtmb_local_registries()` (R/registries.R) put back the saved hook
lists at the restore, which threw away a load hook another party set
inside the scope. It now removes only its own hooks, by identity, and
schedules the restore before it sets any hook (the review's two notes,
`dev/reviews/2026-10-14-testiso.md`). The proof in
test-registry-helper.R compared `length()` of the compat rules; it now
compares the rule list itself, so a change of order (precedence) fails.

New test "a load hook another party sets inside the scope survives it"
(child process through `child_with_libs()`; frmtmb.latent not loaded at
the start, which the test asserts). Before (rellib-r21,
`reg-base.log`):

    RESULT test-registry-helper.R pass=8 fail=2 err=0 skip=0 warn=0
    DETAIL FAIL [...] `actual`: "HOOKS 0 " | `expected`: "HOOKS 1 "
    DETAIL FAIL [...] "THIRD PARTY HOOK RAN" %in% out ... FALSE

After (lane library, `reg-mine.log`): pass=10 fail=0.

## 3154: a reload duplicated frame checks and prior-default providers

`frmtmb_register_frame_check()` and `frmtmb_register_prior_defaults()`
now go through `register_once()` (R/registries.R): a function registered
again under the same key replaces its entry in place. The key is
superseded by fix round 1 (F1, below): it is now the owner's name, the
deparsed code and the serialized values of the function's own
environment. The first key (namespace plus the name found by scanning
every binding, or code alone) collapsed two closures one factory made
and forced each namespace's lazy-load database.

Measured with `dev/hygiene-reload.R` (`dev/hygiene-out/reload.log`):
frmtmb.eam, frmtmb.ode, frmtmb.sample loaded, then
`pkgload::load_all("extensions/frmtmb.ode")` twice, then frmtmb.sample's
`.onLoad()` again.

    arm base
    eam, ode, sample loaded      checks=3 providers=1
    load_all(frmtmb.ode) 1       checks=4 providers=1
    load_all(frmtmb.ode) 2       checks=5 providers=1
    frmtmb.sample .onLoad() again checks=5 providers=2
    arm mine
    eam, ode, sample loaded      checks=3 providers=1
    load_all(frmtmb.ode) 1       checks=3 providers=1
    load_all(frmtmb.ode) 2       checks=3 providers=1
    frmtmb.sample .onLoad() again checks=3 providers=1

New `tests/testthat/test-register-once.R` (core only, inside
`local_frmtmb_registries()`), with the absent case (two checks of
different code stay two). Before (`reg-base.log`): pass=4 fail=4 err=1
(the error is `registration_key` not found, the weak form; the four
failures are the behavioral ones, lengths 3 and 5 for 2 and 3, providers
2 and 3 for 1 and 2). After (`reg-mine.log`): pass=9 fail=0.

## 131: constant weights on a gaussian fit

The behavior was right and unpinned. New
`tests/testthat/test-weights-constant.R` (no extension, no gate): a
constant `weights()` of 2 on `y ~ x + z` and of 3 on `bf(y ~ x + z,
sigma ~ x)` leaves every estimate (sigma's included) within 1e-3 of its
own SE of the unweighted fit, and the logLik is 2 (3) times the
unweighted one within a relative 1e-8. Passes on rellib-r21
(`wc-base.log`: pass=8 fail=0): the item asked for the pin, not a fix.

Margins over 20 perturbations of y by `1 + 2e-15 * k`, k = 0..19, seed
31 (`dev/hygiene-wconst.R`, `dev/hygiene-out/wconst.log`):

    MAX est/se 1.65e-05 (bound 1e-3)  MAX rel logLik 3.15e-12 (bound 1e-8)
    logLik at k=0: unweighted -82.0586, w=2 -164.1172

The sentence is in `vignettes/inputs.Rmd` (the `weights()` bullet), not
in `?bf`: the `?bf` paragraph on weights is the one lane wcheck2's
unbounded-weights rule lives in, and a hunk there would collide. It says
that in a maximum-likelihood fit with no random effect and no prior a
constant weight multiplies the log-likelihood and leaves every estimate,
sigma included; that a prior, a random effect or `REML = TRUE` changes
that; and that lme4's `weights =` are precision weights, under which
sigma moves. (The first wording said "with no random effect" only; fix
round 1, F3, scoped it.)

## 1222: REML x mi() registered "works", REML x me() "conditional"

The registrations disagreed and the `mi()` one was wrong in part. In
`frm_compat()`, `mi()` is the response addition term and `mi_pred()` the
predictor; REML x `mi_pred()` had no rule (it resolved to "untested,
See the specific rules").

Measured with `dev/hygiene-reml-mi.R` (seed 7, n = 120, 12 groups, 20
missing x; `dev/hygiene-out/reml-mi.log`, rellib-r21). Laplace is exact
exactly when the integrand is gaussian in what it integrates, so the
instrument is how far the inner Hessian moves between the inner mode and
a point 0.3 away in every random coordinate (relative to its largest
entry). `spHess()` writes into one sparse matrix in place, so the first
spelling compared the matrix with itself and read 0 on the positive
control too; it copies now.

    REML y | mi() with no reader minus complete case logLik: 0
    inner Hessian relative move, 20 latent x: 0.0342
    inner Hessian relative move, nothing missing: 0
    control, poisson GLMM: 0.341
    random: b 12 beta 4 miss 20

So with a predictor `mi(x)` reading latent values, REML integrates the
slope with the values it multiplies (`beta` and `miss` are both in the
random vector) and is a Laplace approximation, as for `me()`. With
nothing missing, and for a missing response no density reads, it is
exact. Compat rows changed (R/compat.R), each on its own line:

- `REML` x `mi()`: works -> conditional, with that note.
- `REML` x `mi_pred()`: new rule, conditional (was the default
  "untested").

No test reads either row (grep of all nine suites).

## 1519: rate() with cens() or trunc() against brms

`dev/hygiene-ratecens.R` (seed 1519, n = 300) runs
`helper-brms.R`'s `brms_lp_check()` (check A: `rstan::log_prob()` at
frmtmb's estimates against `logLik()`; check B: brms's gradient there)
on brms 2.23.0's compiled programs, rellib-r21
(`dev/hygiene-out/ratecens.log`):

    RESULT poisson rate+cens const 5.12e-13 grad 6e-05
    RESULT poisson(sqrt) rate+cens const -5.68e-14 grad 0.000114
    RESULT negbinomial rate+cens const 0 grad 1.18e-07
    RESULT geometric rate+cens const 0 grad 0.000332
    RESULT poisson rate+trunc(lb=0,ub) const 0 grad 4.42e-06
    RESULT poisson rate+trunc(lb=1) const 1.14e-13 grad 0.000586
    RESULT negbinomial rate+trunc(lb=0,ub) const 1.14e-13 grad 1.46e-06

No defect: frmtmb agrees with brms on every shape. Compat rows changed:
`rate()` x `cens()` untested -> works, `rate()` x `trunc()` untested ->
works, each naming what was verified. The gated test is a new block
"row 24b" in `tests/testthat/test-brms-likelihood.R` (six of the seven
shapes; geometric is covered by the negbinomial program at shape 1 and
by row 24). It cannot fail on rellib-r21: it pins agreement that was
already there.

## 3045: the `lcdf` slot

Docs only. The `@param lcdf` of `frmtmb_family()` now says that despite
the `l`, the slot returns `F(q)` on the probability scale, and names the
two log-scale slots. The internal `fam_lcdf()` comment said "Call a
family's log-CDF", which is the misreading itself; it now says what the
slot returns. Not done: renaming the slot (see "Questions for the user").

## 4259: the bin-2 survey's sar() and fcor() rows

The note is generated, not hand-typed: `dev/brmsport-bin2-survey.R`
gave the 13 rows `bp_absent("sar")` / `bp_absent("fcor")` ("frmtmb has no
sar term ..."). 0.76.0 ported all 13 (the ledger reads them as pass);
the survey's class is read by the ledger only for a row that is not
ported, so the text only misled. The five list entries and `bp_absent()`
are gone and the comment says why. Rerunning the survey before the edit
reproduced both tsv files byte for byte; after it, exactly 13 rows of
`dev/brmsport-bin2-survey.tsv` change (class and reason now empty: ids
log_lik:150, :152, :157, :159, :174, :176, :230, posterior_predict:111,
:113, :118, :120, :133, :135) and `dev/brmsport-bin2-groups.tsv` does not
change. No ledger row reads them.

## 4448: the bin-2 adapters' `../../brms-suite`

`brms_adapt_suite_dir` (dev/brmsport-bin2-adapters.R) now walks up from
the working directory to the nearest directory holding
`dev/brms-suite/brms/tests/testthat` when `FRMTMB_PORT_ROOT` is unset,
and otherwise yields a path that names the variable, so the `source()`
that reads it fails with the cure. Measured with `dev/hygiene-bin2-run.R`
(a copy of `dev/brmsport-bin2-run.R` whose stage directory and adapters
file are variables), `FRMTMB_PORT_ROOT` unset, rellib-r21, the
distributions file staged four levels down:

    before (old adapters)  pass=88 fail=4 err=0  (dsratio/dcratio,
      inv_link_sratio/cratio, dcategorical, inv_link_categorical)
    after                  pass=92 fail=0 err=0

## 3042: staleness lapsed through a setup line (R7)

`brms_setup()` (tests/testthat/helper-brms-suite.R, and its generated
copy in frmtmb.sample with the generator's three-line header) now reads
the stale names before the line runs; a line that reads one and succeeds
marks what it assigns stale, with the STALE message, instead of fresh.
The setup row is still recorded as having run, so no recorded setup
verdict moves. New `tests/testthat/test-port-harness.R` is the
construction of `dev/remeasure-bin2-stale-setup.R` plus the absent case
(a fresh reassignment clears it). Before (the 0.82.0 helper restored,
`ph-before.log`): pass=9 fail=2, x:5 and x:9 held TRUE. After
(`ph-after.log`): pass=11 fail=0.

Reach: `dev/remeasure-stale-reach.R` on the bin-1 records in this tree
(865 rows, 177 blocks): "held assertions that read a name tainted
through a setup line: 0" (`dev/hygiene-out/stale-reach.log`). The bin-2
records are not in the tree; the bin-2 files were run with the new
helper instead (see the suite runs).

## 5516: refusal correct on a conditional row, counted

`dev/equiv-summary.R` prints one more generated line under the
definition of "agree": the pairs scored "refusal correct" whose row is
"conditional", with their ids. The record is not touched: the summary
was rerun into a scratch outdir from the recorded results, split back
into per-pair files by `dev/hygiene-equiv-split.R`
(`dev/hygiene-out/equiv-sum/summary.md`), and `git status` shows
`dev/equiv-data/` and `dev/equiv-findings.md` unchanged.

    0.82.0 record:  Refusal correct on a conditional row (...): 0 unique (0 draws).
    control (p010's status set to conditional in a scratch copy):
                    ... 1 unique (1 draws): p010 cens() x disc.

## 3164: extensions' once-per-session state

No test seam was added to the extensions (package code in three
packages, for state no test depends on today). The guard is in the
one-session runner, where a dependency between files would show:
`dev/release/run-session.R` snapshots every top-level environment of
each loaded `frmtmb.*` namespace (found, not listed: today
`frmtmb.ode::ode_warned`, `ode_failure_log`, `frmtmb.eam::gd_cache`,
`frmtmb.sample::draws_call_state`) and prints each change as
`LEAK <file> <pkg>::<env>: <names> [extension session state]`. The tag
keeps the checklist's rule (a LEAK line has a tag) and names the file
that filled a once-per-session cache. The RNG is not compared: most
files seed it, so a line per file would carry no signal. Exercised by
the session runs below.

## 3177: the coverage matrix runs only what changed

`.github/workflows/test-coverage.yaml` gains a job `plan` that runs
`.github/scripts/coverage-plan.py` and gives the core job an `if:` and
the extension job its matrix (`include: ${{ fromJSON(...) }}`). The
script lists the changed files (`git diff --name-only` between
`github.event.before` and `github.sha` on a push, the base and the merge
commit on a pull request) and keeps a job when a file it reads changed:
core reads everything outside `extensions/` except `dev/`, `docs/`,
`pkgdown/`, `*.md`, `codemeta.json`, `_pkgdown.yml` and the other
workflows; an extension reads that plus its own directory and its
siblings'. When the list cannot be made (a new branch, a force push),
every job runs. The matrix rows moved into the script, with their
comments.

Codecov's join, read with its public API (read only, nothing posted):
at 0.81.0 (e8abed0a) the core flag holds 79 files at 91.3 percent, the
frmtmb.eam flag 15 files at 84.86 percent, all named
`extensions/frmtmb.eam/R/...`, and the commit 143 files over 9 sessions.
So the path rewrite in the extension job joins as intended. Skipping a
job leaves its flag without an upload for that commit, so
`.github/codecov.yml` turns on carryforward for every flag (Codecov
reads that location; the root would need an `.Rbuildignore` line).

Dry run of the script on constructed file lists; the 60-commit line is
the push diff the workflow computes (first parent to commit, no rename
detection), read with `git show --diff-merges=first-parent --no-renames
--name-only` (`dev/hygiene-out/plan-dryrun.log`). The first version of
this line used `git show --name-only`, whose merge output is the combined
diff, and said 350 and 38 (review F5); the review's bare-clone run of the
script's own `git diff` gives the same 398 and 43:

    dev-only (dev/a.R, dev/b.md, NEWS.md)     core false, matrix []
    another workflow (pkgdown.yaml)           core false, matrix []
    this workflow                             all nine
    extensions/frmtmb.spatial/R/car.R         spatial only
    extensions/frmtmb.eam/R/x.R               eam, learn, sample
    ac4906bd (one frmtmb.sample test file)    eam, sample
    before = 000... (new branch), or no env   all nine
    last 60 first-parent commits of main      398 of 540 jobs; 43 need all nine

The YAML parses (`yaml::read_yaml()`), but the plan job has never run on
GitHub: the first push is its test. `fromJSON()` of an empty matrix is
guarded by the job's `if:`.

## Items not done, and why

- **1166** (wiener coverage of log sd(log bs | s)): a coverage study.
  One replicate at 30 subjects by 400 trials takes a median of 261 to
  424 s (dev/phase3b-findings.md:657, :741, :765); a one-sample design
  that can include or exclude 95 percent needs about 200 replicates
  (lane-rules.md, "Power for a coverage study"), about 17 process-hours.
  Out of an S lane's budget with three processes on a shared machine.
  Left open.
- **2733** (fixture 1's theta_1 and theta_2 under OpenBLAS LAPACK) and
  **3018** (the gr(g, by = f) fixture pinned 1.44e-4 below the maximum):
  both are SE-verdict tests in test-se-flags.R, test-se-check.R and
  test-se-tier.R, the files of lane selimit2's area (R/se-check.R),
  which is changing the verdicts this round. A test written against
  0.82.0's verdicts would be merged against selimit2's. In 3018 the
  stop below the maximum is the subject of the test ("a small loading
  on a downward direction"), given by an explicit start; replacing it
  needs a fixture that reaches such a point unforced. Left open for
  selimit2 or the round after it.
- **3045** (the `lcdf` slot): documented, not renamed. See the question.
- **3309** (bin-2 tier: logistic-normal correlations, category order):
  brms's logistic-normal blocks set `prep$lncor`, which neither brms
  nor the adapter reads (review 2026-10-14-bin2port.md, finding 6), so
  no brms row can carry the path; a tier check needs a value check of
  frmtmb's own with a prep brms does read, and the simulator rows are
  shape-only by design. M, in the ledger machinery that is regenerated
  once at consolidation. Left open.
- **3649** (brms tests.posterior_epred.R, tests.stancode.R): the suite
  audit put both in bin 3 by design (dev/brms-suite-audit.md:80, :89):
  stancode is 74 blocks and 756 assertions, 708 of them matches against
  generated Stan, and posterior_epred's 37 read brmsprep internals. The
  audit's own costing of harvesting stancode's 23 refusal assertions is
  0.25 of a lane (:470). Not S. Left open.
- **4628** (redraw the survey sample): the sample lives in
  `dev/equiv-data/`, which lanes must not edit; the item says to redraw
  "at the next frame change", a consolidation decision. The other parts
  (`inv_logit()` in an nl body, translator gaps) are features or 5533.
  Left open.
- **5533** (translator gaps, M): six independent gaps in the survey's
  translator; not reached.

## Questions for the user

1. **Rename the `lcdf` slot of `frmtmb_family()` to `cdf`?** It returns
   probabilities while its siblings `lccdf` and `cdf_log` are on the log
   scale; the docs now say so. A rename touches every family in core
   and the extension families (frmtmb.eam, frmtmb.survival, the
   families-extreme file of lane famend), so it belongs in a round
   with no families lane. Recommendation: yes, with `lcdf` kept as an
   error that names `cdf` for one release (no existing users, so no
   alias).

## 3631: the translator read uncorrelated per-category blocks wrong

Diagnosed and fixed in the test helper (no package code). `(1 | g)`
under `multinomial()` gives one VarCorr block per category, "muy2: 1 |
g" and "muy3: 1 | g", each with one column named "(Intercept)"; the
`(1 | ID | g)` block names its columns "Yc.muy2:(Intercept)". Pooled by
`brms_group_pars()`, the bare name matched the first block for every
category (`brms_frm_coef()`'s `want %in% cols` branch), so each
category's z was computed with muy2's sd. New `brms_block_qualify()`
(tests/testthat/helper-brms.R) gives a column name the pooled blocks
repeat its block's dpar ("muy3:(Intercept)", the form `brms_frm_coef()`
already reads); unique names are left alone, so no other shape moves.

`dev/hygiene-3631.R` (seed 65, n = 300, three categories, 15 groups,
multinomial; rellib-r21; logs `dev/hygiene-out/3631-before.log`,
`3631-after.log`):

    before (0.82.0 helper)  (1 | g):       check C off by 4.199 (tol 0.001)
                            (1 | ID | g):  const 9.09e-13, grad 2.88e-12
    after                   (1 | g):       const 0, grad 9.2e-12
                            (1 | ID | g):  const 9.09e-13, grad 2.88e-12

New gated block in tests/testthat/test-nb2dm-brms.R, "uncorrelated
per-category (1 | g), on the joint density" (multinomial and
dirichlet_multinomial on `nbdm_dm_data(65)`), and the comment of the
`(1 | ID | g)` block that said check C cannot speak for `(1 | g)` is
gone. frmtmb.spatial keeps its own, already diverged copy of
helper-brms.R; its car()/sar() tests have no per-category blocks, so it
is not changed.

The test file, gated, lane library (`nb-before.log` runs it beside the
0.82.0 helper-brms.R in a scratch directory; `nb-after.log` in the
tree):

    before  RESULT test-nb2dm-brms.R pass=22 fail=4
      multinomial: check C 8.043 >= 0.002, gradient 2.455
      dirichlet_multinomial: check C 126.675 >= 0.001, gradient 1.454
    after   RESULT test-nb2dm-brms.R pass=26 fail=0

These are the item's own 8.04 and 126.7.

Positive control of the snapshot (`dev/hygiene-extstate.R`, which takes
the two functions out of the runner's source): it reads
`frmtmb.ode::ode_failure_log`, `frmtmb.ode::ode_warned`,
`frmtmb.eam::gd_cache` and `frmtmb.sample::draws_call_state`; with
nothing changed it prints nothing, and after a key is set in
`ode_warned` it prints
`LEAK probe.R frmtmb.ode::ode_warned: probe-key [extension session state]`.

## Test runs (lane library; every `lib:` line names wt-hygiene-lib for
## frmtmb)

One file per process, `NOT_CRAN=true`, no other gate
(`dev/hygiene-par.sh plain-mine`, `dev/hygiene-out/plain-mine.log`,
533 of 533 files ran):

    frmtmb           files=334 pass=24428 fail=0 err=0 skip=241 warn=0
    frmtmb.coupling  files=11  pass=551   fail=0 err=0 skip=5   warn=0
    frmtmb.eam       files=33  pass=1851  fail=0 err=0 skip=5   warn=0
    frmtmb.latent    files=13  pass=697   fail=0 err=0 skip=2   warn=0
    frmtmb.learn     files=16  pass=474   fail=0 err=0 skip=13  warn=0
    frmtmb.ode       files=11  pass=549   fail=0 err=0 skip=1   warn=0
    frmtmb.sample    files=91  pass=4238  fail=0 err=0 skip=18  warn=0
    frmtmb.spatial   files=8   pass=529   fail=0 err=0 skip=6   warn=0
    frmtmb.survival  files=16  pass=575   fail=0 err=0 skip=0   warn=0

Every skip is a gate (201 brms fit, 52 Stan tier, 15 scale, 13 drmTMB,
2 slow sampler, 1 fuzz, and 6 brms rows of the same gates). The ported
bin-1 suite is NOT in the plain tier: `skip_unless_brms_suite()` needs
`FRMTMB_BRMS_FIT_TESTS=true`, and its 15 files report pass=0 skip=0 there
(a file-level skip). The first version of this record said otherwise
(review F4). They ran gated on the final helper in the review
(`dev/reviews/hygiene-review-suite.sh mine`,
`dev/reviews/hygiene-review-suite-mine.log`), every file fail=0 err=0
skip=0: frmtmb brm 23, brmsfit-helpers 3, brmsformula 16, brmsterms 5,
data-helpers 6, emmeans 11, families 84, methods 162, priors 36,
standata 87; frmtmb.sample brmsformula 16, brmsterms 4, families 84,
helper-copy 36, methods 61. A pass row turned STALE would fail there, so
no bin-1 pass row moves.

Gated (`FRMTMB_BRMS_FIT_TESTS=true`), every file that reads the
translator or that this lane changed (`gated-mine.log`,
`gated2-mine.log`, `nb-after.log`; 0 skips in each):

    test-brms-likelihood.R 556, test-brms-lik-cdfs.R 70,
    test-brms-lik-gp-exact.R 38, test-fcor.R 57, test-gp-by.R 86,
    test-latres-brms.R 281, test-mi-mm-brms-l1.R 41,
    test-mo-factor-re.R 44, test-mo-terms.R 74, test-nl-brms.R 91,
    test-refcat-brms.R 17, test-simplex-brms.R 15, test-nb2dm-brms.R 26,
    frmtmb.sample test-loo.R 146; all fail=0 err=0 warn=0

Bin 2, every generated file in tier mode with this tree's helpers and
adapters, `FRMTMB_PORT_ROOT` unset, staged at another depth
(`dev/hygiene-bin2-all.sh`, `dev/hygiene-out/bin2-all.log`):

    frmtmb distributions 92, log_lik 105, methods 58,
    posterior_predict 65, priors 1; frmtmb.sample methods 58;
    all fail=0 err=0 skip=0

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
lapack`, deleted after): test-weights-constant.R pass=8 fail=0, and
`dev/hygiene-wconst.R`'s 20 perturbed draws give the same margins as
the reference BLAS (MAX est/se 1.65e-05, MAX rel logLik 3.15e-12;
`wconst-ob.log`). No other new assertion reads an SE or a warning.

R CMD check --as-cran (`dev/hygiene-check.sh`, built and checked in
`dev/hygiene-check/<pkg>/`, R_LIBS lane, rellib-r21, user):

    frmtmb.ode     Status: 1 NOTE (V8 math rendering); tests 45 s OK
    frmtmb.sample  Status: 1 NOTE (V8 math rendering); tests 276 s OK
    frmtmb         Status: 2 NOTEs (examples over 5 s: frm_hazard_reads
                   6.27 s, dharma_residuals 6.48 s elapsed, which is load
                   on a shared machine (lane-rules.md); V8); tests 67 min OK

One session per suite, three orders (`FRMTMB_RELLIB=<lane lib> NOT_CRAN=true
bash dev/release/run-sessions.sh hygsess 2`; `dev/release/hygsess.log`,
"HYGSESS ran 27 of 27"):

    frmtmb   forward files=321/334 pass=24428 / reverse 24593 /
             preload 24697, fail=0 err=0 warn=0 in each
    coupling 551, eam 1851, latent 697, learn 474, ode 549, spatial 529,
             survival 575 in all three orders, fail=0 err=0 warn=0
    frmtmb.sample reverse and preload files=86/91 pass=4238 fail=0 err=0;
             forward pass=4236 err=1

The one error is test-parallel-chains.R, `socketConnection(... port =
port ...)`: "cannot open the connection". My two sessions ran at once
without distinct `R_PARALLEL_PORT` values (lane-env's rule, missed by
this lane); the file passes one file per process (plain tier) and in the
other two orders. Rerun of that session alone with `R_PARALLEL_PORT` set:
see below. The 13 core files and 5 frmtmb.sample files the sessions do
not count (321/334, 86/91) are the ported-suite files and their
like, which take 0 s in a session; the 0.82.0 round's records show
the same shortfall of 13 (305/318 and 304/317 in
`dev/arspat-findings.md`, `dev/quadseam-findings.md` and others).

No LEAK line without a tag. The new tagged lines, all expected:

    3  frmtmb.eam::gd_cache: tri                                (test-gddm-gradients.R)
    12 frmtmb.ode::ode_failure_log: when groups penalty n_groups
    6  frmtmb.ode::ode_failure_log: when n_groups
    3  frmtmb.ode::ode_warned: states8                          (test-ode.R)

Order of the runs against the edits: the plain tier and the core
forward and reverse sessions ran before the 3631 helper change and the
new nb2dm block; the core preload session, the core R CMD check (its
whole suite in one session, OK), the gated files above and `nb-after`
ran after it.

Rerun of the frmtmb.sample forward session alone with
`R_PARALLEL_PORT=11873` (`dev/hygiene-out/sample-forward-rerun.txt`):
`SESSION frmtmb.sample forward files=86/91 pass=4238 fail=0 err=0
skip=18 warn=0`, no untagged LEAK line.

## Version and floor

Patch is enough for this lane alone: no function is added to the
exported surface, two compat rows change status, and two extension-seam
functions stop duplicating on a re-registration. No extension calls
anything new, so no `frmtmb (>= )` floor moves. frmtmb.ode and
frmtmb.sample change only in tests (no NEWS entry, no version bump
needed by this lane).

## Files touched outside the obvious area

R/compat.R (four rule lines: REML x mi(), REML x mi_pred(), rate() x
cens(), rate() x trunc()), R/families.R (two comment blocks, no code),
R/structure.R and R/sampling-api.R (one registration line each plus
comments), tests/testthat/test-brms-likelihood.R (one new block after
row 24), tests/testthat/test-nb2dm-brms.R (one comment removed, one
block added), tests/testthat/helper-brms.R (one new function, one line
in brms_group_pars()), vignettes/inputs.Rmd (one bullet).

## Fix round 1 (review `dev/reviews/2026-11-13-hygiene.md`, verdict FIX)

### F1: register_once() collapsed two closures from one factory

The first key was the namespace plus the name found by `get()`-ing every
binding of the closure's namespace, else the code alone. Two closures
one factory made for different values have the same code, so the second
replaced the first and its refusal was lost; the scan also forced the
namespace's lazy-load database (+50.4 MB in the review).

New `registration_key()` (R/registries.R): a list of the owner's name
(`topenv()`), the deparsed code (no srcref), and, unless the function's
own environment has a name (a namespace, the global environment), the
serialized values of that environment (`as.list(..., sorted = TRUE)`).
No binding of any namespace is read. A value that cannot be serialized
makes the key unique, so such a function is never collapsed. Only the
function's own environment is read: a value one level further up (a
factory inside a factory) is not part of the key, which the roxygen
says. Keys compare with `identical()`. Docs (`?frmtmb_register_frame_check`,
the comment in R/structure.R, R/sampling-api.R) and NEWS now say what
counts as the same: the same code with the same captured values.

The review's reproducer, `dev/reviews/hygiene-review-regonce.R`, on
rellib-r21 and the fixed lane build (`dev/hygiene-out/regonce-base.log`,
`regonce-mine.log`):

    base  FACTORY checks added: 2; fit refused by check refuse_negative_x
          FACTORY providers added: 2
          LOAD eam+ode+sample 0.27 s, +4.6 MB; checks=3 providers=1
          RELOAD ode: checks 4; RELOAD sample: providers 2
          COST one registration from a test env: 0.00 s, +0.0 MB
    mine  FACTORY checks added: 2; fit refused by check refuse_negative_x
          FACTORY providers added: 2
          LOAD eam+ode+sample 0.37 s, +4.5 MB; checks=3 providers=1
          RELOAD ode: checks 3; RELOAD sample: providers 1
          COST one registration from a test env: 0.00 s, +0.0 MB
          (first key, from the review: +11.6 MB on the load, 0.17 s and
          +50.4 MB on the registration)

`proc.time()` ticks at 10 ms, so a single registration reads 0.00 s;
1000 re-registrations of one closure from a testthat environment, with
eam, ode and sample loaded (`dev/hygiene-regcost.R`,
`dev/hygiene-out/regcost.log`):

    base  0.00 s, +0.0 MB, checks 3 -> 1003 (appended every time)
    mine  0.85 s (0.85 ms each), +0.1 MB, checks 3 -> 4

`pkgload::load_all("extensions/frmtmb.ode")` twice and frmtmb.sample's
`.onLoad()` again (`dev/hygiene-reload.R mine`, `reload-f1.log`): checks
3, 3, 3 and providers 1, 1, 1, as before the fix.

test-register-once.R gains the factory pair as the case that must
register twice (a refusing and a logging check; the fit is refused; a
third registration of the refusing value is a reload and replaces it),
and two providers from one factory. Against the first key put back into
the lane build (`dev/hygiene-oldkey.R`, `oldkey.log`) that block fails
3 assertions. On rellib-r21 (`ro-base-f1.log`) the file gives pass=6
fail=6 err=1: the reload halves fail, the factory pair's first part
passes, as base appends. Lane build: pass=14 fail=0.

### F2: the coverage plan and renames

`.github/scripts/coverage-plan.py` calls `git diff --name-only
--no-renames`, and `neutral()` no longer treats a `.md` file under
`tests/` or `extensions/*/tests/` as unread (testthat snapshots). In a
scratch repository in the session scratchpad (outside this repository,
deleted after), `git mv R/y.R extensions/frmtmb.eam/R/y.R`: git's
default diff prints only `extensions/frmtmb.eam/R/y.R`; the fixed script
reads 2 files and runs core and all eight extension jobs (the review's
log had core false, eam/learn/sample only). A commit adding only
`tests/testthat/_snaps/print.md` now runs core and all eight. The
60-commit dry run (first parent, no renames, `plan-dryrun.log`):
398 of 540 jobs, 43 commits need all nine, the review's numbers.

### F3: the constant-weight sentence

`vignettes/inputs.Rmd` and NEWS now say "a maximum-likelihood fit with
no random effect and no prior", and that a prior, a random effect's
density and the REML term of the integrated fixed effects are not
multiplied, with sigma^2 going from RSS / (n - p) to c RSS / (c n - p)
under REML (the review's `hygiene-review-wreml.log`: sigma ratio 0.9904
at w = 2, equal to sqrt(2(n-p)/(2n-p))). The test file's header says the
same; its assertions are ML only and unchanged.

### F4, F5, nits

The record now quotes the review's gated bin-1 run for 3042 and the
398/43 dry run for 3177 (above, in place). The 3164 CLOSED line says the
RNG half was declined and why. Checklist step 3 no longer refers to a
version-bump step; it says to regenerate once core's DESCRIPTION
carries the release version.

### Reruns on the final tree (lane library, every `lib:` line names it)

One file per process (`fix1-mine.log`, 75 of 75): core
test-register-once 14, test-registry-helper 10, test-structure 128,
test-compat-register 103, test-conditions 150, test-core-boundary 10,
test-get-prior-route 82, test-prior-table-brms 96,
test-weights-unbounded 160, test-weights-constant 8, test-codemeta 4;
every frmtmb.eam file (33 files, 1851 pass, 5 gate skips), every
frmtmb.ode file (11, 549, 1), and the 20 frmtmb.sample files that touch
frame checks or prior defaults (1000 pass, 1 gate skip); all fail=0
err=0 warn=0.

Forward sessions, each with its own `R_PARALLEL_PORT`
(`dev/hygiene-sessions-fwd.sh fix1`, `dev/hygiene-out/sess-fix1.log`):

    frmtmb   files=321/334 pass=24433 fail=0 err=0 skip=242 warn=0
             (5 more than before: the new test-register-once.R block)
    sample 4238, eam 1851, coupling 551, learn 474, ode 549, spatial 529,
             survival 575; all fail=0 err=0 warn=0
    frmtmb.latent: NO SESSION LINE in that run

The latent session's log stops after `FILE test-mixture-mvn.R` with no
error text and no SESSION line, so the process ended from outside the
test (no exit code is kept by the driver). In that run test-hmm.R took
11.9 s against 80.9 s, loaded fewer packages, and left
`LEAK test-hmm.R search path: +AD-overloads` (RTMB's taping
environment), the only untagged LEAK line of the round. Not reproduced:
the latent session rerun alone (`R_PARALLEL_PORT=11977`,
`sess-fix1/frmtmb.latent-rerun.txt`) gives `files=13/13 pass=697 fail=0
err=0 skip=2 warn=0` and no untagged LEAK, as did all three latent
orders of the first round. Recorded as an unexplained process death on a
loaded machine; worth a look if it recurs at consolidation.

R CMD check --as-cran of core on the final tree (R/ changed: the key in
registries.R and the roxygen of `frmtmb_register_frame_check()`):
`Status: 2 NOTEs`: examples over 5 s (`frm_hazard_reads` 5.51 s
elapsed, load) and V8; tests 39 min OK. frmtmb.ode and frmtmb.sample
changed nothing in the fix round, so their earlier checks stand.

## Fix round 2 (short check `dev/reviews/2026-11-13-hygiene-short.md`, F1 only)

The fix-round-1 key read only the closure's own environment and
serialized environments by their contents, so it still merged: a
factory inside a factory with different outer values, `local()` inside
one factory, two closures over two distinct (empty) environments, and
two closures over two distinct external pointers. Each merge dropped the
refusing check (short check, `hygiene-short-regkey-mine.log`).

Rule now (decided by the coordinator): never merge unless provably the
same; when in doubt keep both. Implementation (a), the short check's
prototype, chosen over (b) owner + id because it needs no new argument
on the exported `frmtmb_register_frame_check()`, no change in any
extension's `.onLoad()`, and it passed every case as it stood:

- `registration_key()` walks from `environment(fn)` up to the first
  named environment (namespace, global, base, empty, or any environment
  with a name), serializing each frame's bindings
  (`as.list(..., sorted = TRUE)`). A `refhook` that is called at all (an
  unnamed environment or an external pointer among the values) or a
  serialization error makes the key a fresh environment, which is
  `identical()` to nothing. The key is the owner name, the deparsed code
  and those frames.
- `register_once()` merges when the new function is `identical()` to an
  entry (the same object) or its key is `identical()` to the entry's
  stored key, and appends otherwise. Each entry carries its key as the
  attribute `frmtmb_registration_key`, computed once at registration, so
  a later registration does not serialize every other entry again.
- No `<<-` (the refhook writes into a small environment).
- `?frmtmb_register_frame_check` (exported) now says exactly what is
  merged and lists what is kept, the bounded extra provider of
  `load_all()` after `library()` included. NEWS, the comments in
  R/structure.R and R/sampling-api.R, and the 3154 CLOSED line say the
  same.

test-register-once.R (rewritten; its factories live in the stats
namespace, as a package's own would, because a factory made inside a
test reaches the test environment, which holds environments, and would
never merge): new blocks "a value one frame further up keeps two checks
apart" (nested factory, nested same values merge, `local()` pair, the
refusal survives) and "an environment or a pointer among the values is
never merged" (two `new.env()`, then the first turned to refuse; two
`stats` native-symbol pointers). On the fix-round-1 build
(`ro-f2-before.log`): pass=14 fail=5, the nested, `local()` and
pointer pairs merged. Final build (`ro-f2-after.log`): pass=19 fail=0.
rellib-r21 (`ro-f2-base.log`): pass=11 fail=8, the reload merges.

The short check's script on the final build
(`dev/reviews/hygiene-short-regkey.R`, logs
`dev/hygiene-out/regkey-{mine,base}.log`; the extension copies went to
the session scratchpad, deleted after):

    case                                     base           final build
    nested factory, outer value differs      2, refused     2, refused
      same, reverse order                    2, refused     2, refused
    local() in one factory                   2, refused     2, refused
    two empty environments, first refuses    2, refused     2, refused
    two environments, different contents     2, refused     2, refused
    two external pointers                    2, refused     2, refused
    same closure, plain then compiled        2              1
    captured helper, plain then compiled     2              2
    same code, srcref and none               2              1
    captured helper, srcref and none         2              2
    load_all(ode) twice                      4, 5           3, 3
    ode unload() + library()                 6              3
    ode unloadNamespace() + library()        7              3
    load_all(sample) twice                   2, 3           2, 2
    sample unload() + library() after it     4              2
    sample unloadNamespace() + library()     5              2

The 2 providers after `load_all(sample)` are the documented bounded
extra: the provider captures `libname`, which `load_all()` and
`library()` pass differently. Every pair the rule keeps is a duplicate,
never a merge.

Cost (same script; base never keys anything):

    case                                       fix round 1         final
    register a closure over 49.6 MB            0.140 s, +50 MB     0.110 s, +50 MB
    register it again (same object)            0.530 s, +100 MB    0.120 s, +51 MB
    a small closure after it                   0.120 s, +50 MB     0.000 s, +0 MB
    10 small closures after it                 1.300 s, +150 MB    0.020 s, +1 MB
    library(eam, ode, sample) after it         0.720 s, +212 MB    0.280 s, +95 MB
                                               (base 0.250 s, +94 MB)

The one registration of the large closure still serializes it once
(0.11 s, a transient +50 MB). `dev/hygiene-regcost.R mine`: 1000
re-registrations of one test-environment closure 0.03 s (0.03 ms each),
3 -> 4 checks (the same object merges by `identical()`). The review's
reproducer (`regonce-mine-f2.log`): factory checks 2, fit refused;
providers 2; LOAD eam+ode+sample 0.17 s, +4.5 MB, checks 3, providers 1;
RELOAD ode 3, sample 1. `dev/hygiene-reload.R mine` (`reload-f2.log`):
3, 3, 3 checks and 1 provider throughout.

Reruns on the final tree (lane library, every `lib:` line names it):

- One file per process (`fix2-mine.log`, 75 of 75): core
  test-register-once 19, test-registry-helper 10, test-structure 128,
  test-compat-register 103, test-conditions 150, test-core-boundary 10,
  test-get-prior-route 82, test-prior-table-brms 96,
  test-weights-unbounded 160, test-weights-constant 8, test-codemeta 4;
  frmtmb.eam 33 files 1851 pass, frmtmb.ode 11 files 549, the 20
  frmtmb.sample files that touch frame checks or prior defaults 1000;
  all fail=0 err=0 warn=0.
- Forward sessions, own ports (`dev/hygiene-out/sess-fix2/`): frmtmb
  files=321/334 pass=24438 (5 more, the new blocks) fail=0 err=0 warn=0;
  frmtmb.sample files=86/91 pass=4238 fail=0 err=0 warn=0; no untagged
  LEAK line.
- R CMD check --as-cran, core: the first run died inside "checking
  tests" with no error text (the log stops at "Running 'testthat.R'"),
  during the machine-wide memory and page-file event; rerun once with at
  most 2 lane R processes: `Status: 1 NOTE` (V8), tests 27 min OK.

## Fix round 3 (short check, re-check after fix round 2)

- The walk of `registration_key()` stopped at any named environment,
  so a closure whose own or parent environment carried a user name
  (`attr(e, "name")`, `attach()`) was keyed by its code alone and a
  refusing check was lost. `registration_named_env()` now stops only at
  a namespace, the global environment, base or empty, and meeting any
  other named environment on the walk makes the key unique (no walk
  past it).
- Docs: `?frmtmb_register_frame_check` and NEWS say "held the same
  values at registration", name the named-environment rule, and qualify
  the reload claim: a top-level namespace function and an anonymous
  closure `.onLoad()` passes straight in merge on reload; a closure
  `.onLoad()` binds to a name first, or one beside a formula, is added
  again each reload (harmless duplicate).
- Nit: the in-place assertion now checks that the first `mk()`
  closure's environment is gone from the two slots, which fails on base
  (base keeps it), not only the slot order.

test-register-once.R with a new block "a named environment that is not
a namespace is never merged" (own env "cfg" with differing values, own
envs "cfgA"/"cfgB", parent named "cfg", and the refusal surviving):

    fix-round-2 build  pass=19 fail=4 (the three pairs merged, the fit ran)
    final build        pass=23 fail=0   (ro-f3-after.log)
    rellib-r21         pass=14 fail=9   (ro-f3-base.log; the reload merges
                       and the replacement check, as base appends)

`dev/reviews/hygiene-short-regkey2.R` on the final build
(`dev/hygiene-out/regkey2-f3.log`): the three NAMED pairs now give
"added 2; fit refused"; every other line is as in the short check's log
(CHAIN and STATE "added 1" are the accepted values-at-registration
cases; ONLOAD named closure 2, anonymous 1, beside a formula 2).
`hygiene-short-regkey.R` (`regkey-f3.log`): nested, local(), two
environments, two pointers all 2 and refused; ode reloads 3, 3, 3, 3;
sample providers 2, 2, 2, 2 (the documented libname extra); cost:
register the 49.6 MB closure 0.110 s +50 MB, again 0.110 s, later small
registrations 0.000 s, library(eam, ode, sample) 0.180 s +95 MB.
`dev/hygiene-reload.R mine` (`reload-f3.log`): checks 3, 3, 3, providers
1 throughout; the review's regonce reproducer: factory 2 and refused,
providers 2, reload ode 3, sample 1.

R CMD check --as-cran, core, final tree (R/structure.R's roxygen
changed): `Status: 1 NOTE` (V8), tests 20 min OK.
