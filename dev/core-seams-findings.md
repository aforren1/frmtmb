# Two core seams the frmtmb.ddm round surfaced

Worktree `frmtmb-wt-core-seams`, branch `wt-core-seams`, from `564e185`
(core 0.50.0). Nothing is committed.

Both defects are the same defect twice: a declaration that is accepted,
does nothing, and says nothing.

---

## Seam 1. A contributed compatibility rule that names nothing

### What was wrong

`frmtmb_register_compat()` stored `features` and `rules` and returned.
A rule side that named no feature simply matched no pair. The resolved
table then answered `untested` for those pairs from the `* x *`
default, which is exactly what it answers when nobody has said
anything. The registrant had no way to tell the two apart.

`frmtmb.ddm` lost two `wiener` rows for a whole round this way: it wrote
`mixture()` where the vocabulary says `mixture`, and it named `dec()`
before that spelling was a feature at all.

### What now happens

`frmtmb_register_compat()` calls `rules()` once at registration and
refuses anything it cannot resolve. It refuses BEFORE committing, so a
refused registration leaves the vocabulary as it found it rather than
half filled with the features of a rule set that never took.

Checked, in order: `features` shape; each kind against the closed set of
nine; a display name registered under a second kind; then, per rule,
each side (`"*"`, `kind:<kind>`, `group:<group>`, or a feature name), a
bare-name self-pair, and the status against the declared five.

### The refusal texts

Unknown feature, with a near miss (the `frmtmb.ddm` case, exactly):

> frmtmb_register_compat(): the rule `wiener x mixture()` names
> 'mixture()', and the registry has no such feature, so the rule would
> match no pair and be dropped without a word. Did you mean 'mixture'?
> Give the feature a kind in features = if this package supplies it, or
> name it in expects = if another package does. frm_compat_features()
> lists the vocabulary.

Unknown feature, nothing close (the `dec()` case before the term was
registered) - the same sentence without the guess:

> frmtmb_register_compat(): the rule `wiener x dec()` names 'dec()',
> and the registry has no such feature, so the rule would match no pair
> and be dropped without a word. Give the feature a kind in features =
> if this package supplies it, or name it in expects = if another
> package does. frm_compat_features() lists the vocabulary.

The near miss is found by matching the parentheses first (dropping or
adding `()` is the commonest miss, and it is the one that happened) and
falling back to `agrep()` at `max.distance = 0.2`, ranked by
`utils::adist()` and capped at two edits. Kept tight on purpose: a
confident wrong suggestion costs more than none, which is why `dec()`
gets no guess rather than being pointed at `se()`.

The four sibling refusals:

> ... the rule `gaussian x kind:fmaily` matches on 'kind:fmaily', and
> there is no such kind of feature. The kinds are: family, covstruct,
> aterm, special, autocor, mode, structure, method, grammar.

> ... the rule `gaussian x group:cdfs` matches on 'group:cdfs', and
> the registry has no such feature group. The groups are: cdf, ...

> ... the rule `gaussian x cens()` declares the status 'refuse', which
> is not one of: works, conditional, refused, broken, untested. The last
> of them is the point of the registry: an absent guard and a passing
> guard look the same from outside.

> ... the rule `wiener x wiener` names one feature on both sides. The
> resolved table holds unordered pairs of DISTINCT features, so this
> matches no pair either. State it in the note of a rule against
> something the feature really meets.

The last of those is the same failure by another route, and the core's
own rule table has been held to it by `test-compat.R` since the registry
was written; contributed rules were not.

### What it costs

Measured, not guessed. The check runs at `.onLoad()` of every
contributing package, so it is on their load time.

The dominant term is assembling the feature vocabulary, which
`frmtmb_compat_features_tbl()` builds by `rbind`-ing one data frame per
feature: about 37 ms on an idle machine, and pre-existing - `frm_compat()`
has always paid it per call. The first version of this check built it
twice per registration. It now builds it once and passes it down, which
took a 40-rule registration (the size of `frmtmb.ddm`'s) from 2.9 to 1.9
vocabulary builds. The residue is one build plus the loop over the
rules.

Vectorizing `frmtmb_compat_features_tbl()` itself would take the 37 ms
down to near nothing and would help `frm_compat()` more than it helps
registration, but it is a rewrite of a function outside this lane's
region and is left alone.

### `expects =`, and why it had to exist

There is one legitimate reason to name a feature the vocabulary does not
have: another package supplies it, and the rule becomes true when that
package loads. `frmtmb_register_compat(expects =)` declares those. They
are exempt from the check and are NOT added to the vocabulary, so the
rule lies dormant instead of putting a feature in the table that nothing
implements.

This is not hypothetical either. **`frmtmb.sample` was the second
instance of the same defect.** Its `hmm x frm_sample` and
`lca x frm_sample` rules name features owned by `frmtmb.latent`, which
`frmtmb.sample` only Suggests. In every session that loaded
`frmtmb.sample` without `frmtmb.latent` - which is every session that
does not use hidden Markov or latent class models - those two rules were
dropped in silence. Loading `frmtmb.sample` alone against the strict
core reproduces it as a refusal:

```
frmtmb_register_compat(): the rule `hmm x frm_sample` names 'hmm', and
the registry has no such feature, ...
```

so `extensions/frmtmb.sample/R/zzz.R` now declares
`expects = c("hmm", "lca")`. That file is outside the ownership list I
was given; the alternative was shipping a core change that breaks a
first-party extension on load, which is worse. It is one argument and a
comment, plus a stale sentence in the neighboring roxygen ("`hmm()` and
`lca()` are still frmtmb's" - they are `frmtmb.latent`'s).

### Decision: `frmtmb_register_aterm()` also registers the compat feature

**Yes.** The ddmvar lane's argument is right and it is structural, not a
convenience: a term the parser accepts and the table cannot describe is
a gap BY CONSTRUCTION. `frm()` would take `y | dec(x) ~ ...` while
`frm_compat("dec()", "cens()")` refused to answer, and no amount of care
by the registrant changes that the two facts came from one registration.
Three further reasons:

1. The kind is not a judgment call. A registered addition term is of
   kind `"aterm"`; there is nothing for the registrant to decide, so
   asking is asking them to repeat themselves.
2. The cost of the entry is nil. A new feature with no rules answers
   `untested` against everything, which is the honest answer and the one
   the registry exists to be able to give.
3. It removes the failure mode rather than reporting it. With the strict
   check of seam 1 alone, forgetting the feature becomes a refusal
   instead of a silent drop - better, but still a step the registrant
   can get wrong. Registering it at the seam that creates the term
   cannot be got wrong.

The duplicate concern is handled where it belongs: `features` entries
are matched against the vocabulary by display name, a repeat under the
SAME kind is a no-op, and a repeat under a different kind is refused
("One display name carries one kind"). So `frmtmb.ddm`'s explicit
`"dec()" = "aterm"` is untouched and cannot error, whether it registers
the term before or after the rules. `frmtmb_compat_features_tbl()`
de-duplicates by name as a backstop, because a second row for one
feature would give it a pair with itself and double every pair it has
with anything else.

One consequence to know: registration order now matters within a
`.onLoad()`. Rules naming `dec()` must be registered AFTER
`frmtmb_register_aterm("dec")` unless the registrant also spells the
feature in `features =`. `frmtmb.ddm` does spell it, so it works either
way; both help pages say so.

---

## Seam 2. `required_aterms` could only say ALL of

### What was wrong

`frmtmb_family(required_aterms =)` is a conjunction. Two `frmtmb.ddm`
families need EITHER `dec()` OR `vint1` for the same datum - which
boundary a trial ended at - and could not say so, so each hand-rolled
the refusal inside its response check.

### What now happens

`required_aterms` still takes a character vector, meaning all of. It now
also takes a LIST: an element of length one is required, an element of
length more than one is a set of spellings any one of which will do.

```r
required_aterms = list(c("dec", "vint1"), "vreal1")
# one of dec or vint1, and vreal1
```

The all-of refusal is byte-identical to before, deliberately: every
family that declares the old spelling, and the docs and tests that quote
it, see no change. One message template covers both, with the group
description computed:

All-of (unchanged):

> needs_vint: the density needs `vint1`, `vint2`, which nothing on
> this response supplies. Write the addition term: y |
> vint(<column>) + vint(..., <column>) ~ ...

Any-of unmet:

> either: the density needs one of `dec` or `vint1`, which nothing on
> this response supplies. Write the addition term: y | dec(<column>) ~
> ...

Mixed, plain requirements first so the sentence does not trail a bare
name off the end of a choice:

> either: the density needs `vreal1`, one of `dec` or `vint1`,
> which nothing on this response supplies. Write the addition term:
> y | vreal(<column>) + dec(<column>) ~ ...

The example formula writes the FIRST alternative of each unmet group,
because a formula has to pick one spelling to be a formula at all. Which
others there are is in the sentence above it. The help page says so.

### What adopting it would delete from frmtmb.ddm

The families were not touched, as instructed. On adoption:

* `wiener()`: the hand-rolled boundary check, which
  `R/wiener-family.R:445` itself calls "the one hand-rolled check left
  in this package". `required_aterms = list(c("dec", "vint1"))`.
* `gddm()`: the `is.null(up)` branch of `gd_check_response()`
  (`R/gddm.R` around line 1086, about twelve lines including the
  multi-line refusal), and its comment explaining why
  `required_aterms` could not express it. Its `req` computation over
  `terms` stays, because the drift components' `vreal` requirements are
  still a conjunction; only the boundary group is added to it, as
  `c("dec", "vint1")`.
* `lba()` gains nothing: `vint1` is its only route, so
  `required_aterms = "vint1"` is already exactly right.

Net: two refusals move from `valid_y` (after the frame is built, and
opaque to the framework) to a declaration `frm()` checks before assembly
and can report on. `frmtmb.ddm`'s NEWS entry at line 268, "One
hand-rolled check remains, and is not `required_aterms`'s fault",
becomes false and should be retired at adoption.

---

## Verification

Every count below was checked against the files on disk, not taken from
a log. The suite driver truncates its summary and re-runs every file; it
does not resume from a log.

* Full core suite, one file per process under `pkgload::load_all` with
  `NOT_CRAN=true`: **105 test files on disk, 105 result lines, 105
  distinct names, zero failures, zero errors, zero warnings.** The set
  difference against `ls tests/testthat/test-*.R` is empty in both
  directions, no name appears twice, and no line carries NO-SUMMARY, a
  crash, or a non-zero return code.

* The four extension suites, one process each under `test_local()`
  against this worktree's core installed into `cs-lib`. Files run
  against files on disk:

  | package | files on disk | files run | result |
  |---|---|---|---|
  | frmtmb.ode | 4 | 4 | 0 fail, 0 error |
  | frmtmb.latent | 4 | 4 | 0 fail, 0 error |
  | frmtmb.ddm | 14 | 14 | 0 fail, 0 error, 1 skip |
  | frmtmb.sample | 10 | 10 | 0 fail, 0 error, 2 skips |

  This is the test of the automatic-feature decision: `frmtmb.ddm`
  declares `"dec()" = "aterm"` explicitly AND registers the term, and
  meets a no-op rather than a duplicate.

* Each extension's `.onLoad()` also run ALONE in a fresh process, which
  is what the cumulative run hides: loading them in one session lets
  `frmtmb.latent` supply `hmm` and `lca` before `frmtmb.sample` asks
  for them, which is exactly the defect described above.

* `R CMD check --as-cran` with `_R_CHECK_CRAN_INCOMING_=false`:
  **`Status: 1 NOTE`, return code 0.** `checking tests ... [19m] OK`
  with `FAIL 0 | WARN 0 | SKIP 26 | PASS 5460`; vignettes rebuilt OK;
  PDF and HTML manuals built; every Rd section OK. The single NOTE is
  `Skipping checking math rendering: package 'V8' unavailable`, which is
  the machine and not the package.

  Put BOTH TinyTeX and pandoc on PATH before checking. Without TinyTeX
  (`~/AppData/Roaming/TinyTeX/bin/windows`) the run reports `1 ERROR,
  1 WARNING, 2 NOTEs`, every one of them downstream of `pdflatex is not
  available`, and the WARNING says "This typically indicates Rd
  problems", which here it does not. Without pandoc the run does not
  reach `check` at all: `R CMD build` dies while re-building all seven
  vignettes, before `check` starts. Pandoc is at
  `C:\Users\adf44\AppData\Local\Programs\Quarto\bin\tools` on this
  machine, and under the RStudio install at
  `C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools`.

### Reading these logs on a shared machine

This machine runs several lanes at once, and two things corrupted
earlier passes. Neither is a defect in the package; both will recur.

* A reviewer stopping processes by matching the session GUID killed an
  `R CMD check --as-cran` mid-tests and shells belonging to other lanes.
  A check that dies this way leaves a log ending in
  `Running 'testthat.R'` with no Status line: absence of a Status is not
  a pass and not a failure. Kill by the script's own name prefix.

* Cancelled runs of the suite driver survived their cancellation and
  kept appending to the one summary log, which read as 132 result lines
  for 105 files. **Always compare the log against
  `ls tests/testthat/test-*.R` by name, in both directions, before
  believing a summary.** One of those contended runs also failed
  `test-perf.R`, whose assertion is wall clock
  (`t_large < 100 * max(t_small, 0.01)`, 4.17 against 4.00) with four
  heavy R jobs on the machine. It passes in the clean run and inside
  `R CMD check`.

## Deliberate omissions

* `extensions/frmtmb.ddm` families unchanged, as instructed. The
  adoption is listed above rather than done.
* No `DESCRIPTION` version bump anywhere, so `frmtmb.sample`'s
  `expects =` requires the frmtmb in this repository rather than a
  released 0.50.0. Whoever cuts the release owns the
  `frmtmb (>= ...)` floor in `extensions/frmtmb.sample/DESCRIPTION`.
* The core's own rule table is not run through
  `compat_check_rules()` at load time. `test-compat.R` already asserts
  the same three properties over it, at test time, where the cost
  belongs.
* Registration does not check that a contributed rule WINS at least one
  pair, which `test-compat.R` does assert for the core table. Doing it
  would resolve every rule against every pair on every `.onLoad()`, and
  the answer depends on which other packages happen to be loaded, so it
  would be both expensive and unstable. The three cheap and stable
  causes of a dead rule - an unresolvable side, a self-pair, and a bad
  status - are refused instead.
* `expects` names are not checked against anything. There is nothing to
  check them against: the whole point is that the owner is absent. A
  misspelling inside `expects` is still a dangling rule, which is why
  the help page says to use it for that and nothing else.

---

# Punch round, 2026-09-05

Worked from `dev/review-core-seams.md` (718 lines). Its section numbering
governs; the nine punch-table rows are the contract. Nothing committed,
worktree HEAD still `564e185`, main checkout `C:/Users/adf44/source/r/frmtmb`
verified at `5dfdd84` and clean at the start and at the end, never written
to. All scratch work under `...\scratchpad\csp-*`; private library
`...\scratchpad\csp-lib`, explicit `--library=` and `lib=` throughout.

## 1. BLOCKING, closed. `frmtmb_register_aterm()` is atomic

`R/parse.R:117-119`. `compat_new_aterm_feature(name)` now runs BEFORE
`frmtmb_aterm_registry$reg[[name]] <- ...` rather than after it. It is the
last thing in the function that can refuse, so a refusal leaves both
registries as it found them, which is the property `frmtmb_register_compat()`
already had and this one did not.

The refusal itself moved into `compat_new_aterm_feature()`
(`R/compat.R:336-352`) so it can speak in the term registrant's own terms.
The generic "one display name carries one kind" is true and says nothing to
a term author; the collision is a grammar fact:

    frmtmb_register_aterm() cannot register 's', because 's()' already
    means something else in a formula: the compatibility vocabulary holds
    it as a 'special' feature. A formula writing s() could not say which
    was meant, so nothing has been registered. Choose another name.

The review's section (f) repro, re-run end to end against the punched core:

    before: 's' in aterm registry: FALSE
    register_aterm('s') -> Error: ... nothing has been registered ...
    AFTER: 's' in aterm registry: FALSE
    AFTER: kind of 's()' in vocabulary: special
    frm(bf(y | s(col) ~ x)) -> frm() refuses: Addition term `s()` is not
      supported (supported: weights(), trials(), cens(), trunc(), se(...

Before the fix that last line read `FIT OK frmtmb_fit`. The split state is
unreachable again.

Tests, `tests/testthat/test-compat-register.R`: "a refused
frmtmb_register_aterm() registers nothing at all" walks `s` (special), `ar`
(autocor) and `mm` (grammar), and asserts after EACH refusal that the name
is absent from `frmtmb_aterm_registry$reg`, that
`names(frmtmb_aterm_registry$reg)` is identical to what it was, that
`frm_compat_features()` is identical to what it was, and that
`registered_aterm_of("s")` is still NULL. A second test pins the wording.

## 2. A same-call duplicate under two kinds is refused

`R/compat.R:298-312`. `features <- features[!duplicated(names(features))]`
dropped the second entry in silence. It now refuses when the repeat carries
a DIFFERENT kind, and stays silent only when the kinds agree, which is the
half that has to stay quiet: a rule set may spell out a term
`frmtmb_register_aterm()` already declared.

    frmtmb_register_compat(features = c(zzq = "family", zzq = "aterm"))
    #> Error: ... registers 'zzq' as an 'aterm' feature, and the same call
    #>   already gives that name the kind 'family'. One display name
    #>   carries one kind.

The test asserts the refusal, that `zzq` did not reach the vocabulary, that
the same-kind repeat is silent, and that it leaves exactly one row.

## 3. An exact vocabulary member beats the approximate match

`R/compat.R:461-464`. The `hit` test in `compat_near_feature()` now runs
both parenthesis directions: `known == paste0(bad, "()")` was missing, so a
name MISSING its parentheses fell through to `agrep()`, which ranked `us`
(one edit) above `se()` (two).

    se     -> Did you mean 'se()'?      (was: Did you mean 'us'?)
    trials -> Did you mean 'trials()'?

The test pins both and asserts `'us'` does not appear in the `se` refusal.
The existing test that `dec()` gets NO guess still passes, so the tight
threshold is undisturbed.

## 4. A rules builder that throws names itself

`R/compat.R:198-206`, with `compat_registrant()` at `R/compat.R:244-250`.
`rules()` is wrapped so the failure names the argument and the registering
package. Verified through a real `.onLoad()`, not only a unit test: scratch
package `csppkgA` under `...\scratchpad\csp-work\csp-plant`, installed
`--no-test-load` into its own library so the refusal is seen at
`loadNamespace()`:

    .onLoad failed in loadNamespace() for 'csppkgA', details:
      error: frmtmb_register_compat() called the rules = builder that
      csppkgA gave it and it failed, which stops the whole registration:
      the builder runs once at registration so that a rule naming nothing
      is refused here rather than dropped later, so it must read nothing
      but its own arguments. It said: the builder read a file that was
      not there

Before, that was the bare line `the builder read a file that was not there`.
The package name is the half that matters out of `loadNamespace()`, where
there is no call stack to read. Tests cover the message and
`compat_registrant()` directly (a namespace, the global environment, and
frmtmb itself, which returns "" so the core never says "that frmtmb gave
it").

## 5. Floors

* `DESCRIPTION:3`, `Version: 0.50.0` becomes `0.51.0`. The check built
  `frmtmb_0.51.0.tar.gz`, so the bump is real and not only declared.
* `extensions/frmtmb.sample/DESCRIPTION:28`, `frmtmb (>= 0.46.0)` becomes
  `frmtmb (>= 0.51.0)`. The review measured this as a hard load failure
  (`unused argument (expects = c("hmm", "lca"))`) against a pristine
  `564e185` core, so the old floor permitted an installation that could not
  load.
* `extensions/frmtmb.sample/NEWS.md`, a bullet saying exactly that.
* No other extension floor touched, per the review's section 7 table.

## 6. `expects =` closed at both ends

The review's closing line: the registry re-resolves on every `frm_compat()`
call with no cache, so a forward-referenced rule goes live the moment its
owner loads whether or not it was declared, and `expects =` was never
checked against anything. Closed one argument wide, at the two moments
where an honest answer exists.

**At registration** (`R/compat.R:177-186`): a name the vocabulary already
holds is refused, because the declaration then exempts nothing.

    frmtmb_register_compat(..., expects = "cens()")
    #> Error: ... names 'cens()', and the registry already has that
    #>   feature, so the declaration exempts nothing. ...

**At query time** (`R/compat.R:1442-1445` and `1743-1747`):
`compat_unresolved_sides()` collects the bare-name rule sides no feature in
the session answers to, and `frm_compat()` attaches them as an `unresolved`
attribute, with class `frmtmb_compat` and a print method
(`R/compat.R:1752`) that names them under the rows.

The representation is the smallest honest one. The row CANNOT exist: there
is no such feature, so there is no pair, and inventing one would put a
feature in the table that nothing implements, which is the lie the registry
exists to refuse. The review's section 2 reaches the same conclusion and
calls the absent row the right choice. So the fact is reported beside the
table rather than dropped. The attribute and the class are attached ONLY
when the session holds one, so a session with no forward reference gets
exactly the plain data frame it got before, and the review's non-regression
result (`identical()` on the core-alone table) still holds.

Documented in `?frm_compat` under "Rule sides this session cannot resolve",
in `?frmtmb_register_compat`'s `expects` parameter, and on the
`frmtmb-extension-api` page.

Core test, `test-compat-register.R`: registers a rule naming `hmm` under
`expects = "hmm"`, then asserts no row names `hmm`, that
`attr(frm_compat("wiener"), "unresolved")` is `"hmm"`, that the object
carries the class, and that printing says so; then registers something that
supplies `hmm` and asserts the attribute, the class and the report all go
away and the pair answers `works`.

Extension test,
`extensions/frmtmb.sample/tests/testthat/test-sampling-ported.R`: with
`frmtmb.latent` loaded, neither `hmm` nor `lca` may remain in the unresolved
set, and both pairs must answer. That is what makes an `expects` typo
findable. A typo in `expects` alone is refused at registration; a typo in
the RULE, with the same typo in `expects`, survives registration and shows
up here. It ran rather than skipped: the sample suite went 881 to 884
passes with skips unchanged at 2, and 3 is exactly this test's expectation
count.

## 7. The pandoc note, corrected

`dev/core-seams-findings.md`, the prerequisite under "Reading these logs on
a shared machine". It said to put TinyTeX on PATH. Pandoc is required as
well, and without it `R CMD build` dies while re-building all seven
vignettes before `check` starts. Both paths are now given. This round's
check ran with `/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows` and
`/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools` on PATH, and
`which pandoc pdflatex` is logged at the top of `csp-check.log` as proof
both were found.

## 8. The remaining punch rows

* **Row 5, `compat_new_aterm_feature()` and `known`.** Applied
  (`R/compat.R:336-341`): the vocabulary is built once in that function and
  handed to both the collision check and `compat_new_features()`. **The
  row's premise did not hold, and this is the measurement rather than a
  guess.** There was never a second build to remove. On this machine a bare
  `frmtmb_compat_features_tbl()` is 16 ms and `frmtmb_register_aterm()` was
  and is 16 ms, which is exactly one build. The review measured 26 ms and
  22 ms on a busier machine and drew the same conclusion in its prose
  ("essentially one build") while its table row says a build is being
  wasted. What the plumbing does buy is that item 1's new collision check
  did not ADD a second build; without it `frmtmb_register_aterm()` would
  now cost two. Four terms still cost four builds, which needs a cache on
  the vocabulary, and that is a change to `frmtmb_compat_features_tbl()`
  outside this lane's region.
* **Row 6, the article.** `compat_kind_article()` at `R/compat.R:259`, used
  by both `compat_new_features()` refusals and by
  `compat_new_aterm_feature()`. `as a 'aterm'` is now `as an 'aterm'`, and
  `autocor` gets the same treatment. Tested both ways round.
* **Row 7, NEWS.** The `required_aterms` bullet in `NEWS.md` now names the
  two constructor behavior changes the review found in its section 4:
  `c("vint1", NA)` was accepted and is now refused, and `list()` was
  refused and is now accepted as no requirements.
* **Review section (e), the near-match pool.** Not a punch row
  ("defensible, undocumented"). Left as it is and now documented, in
  `?frmtmb_register_compat`'s `expects` parameter: the guess is searched
  over the vocabulary only, never over `expects`, because pointing at a
  feature that is not in the session would be its own confusion.
* **Review section 5, the frmtmb.ddm adoption.** NOT done, deliberately,
  and it is not a punch row. It is a change to another package's families
  that breaks four pinned expectations in two of its test files, which the
  review enumerates. Left for whoever adopts it, along with the two
  `test-gddm-family.R:32-33` pins the review adds to the lane's list.

## Verification

One process per file throughout. Every count below was reconciled against
the files on disk by name, in both directions, before being believed.

**The five named files, plus the new tests, one process each**

| file | fail | error | warn | skip | pass | review |
|---|---|---|---|---|---|---|
| `test-compat-register.R` | 0 | 0 | 0 | 0 | **87** | 38 |
| `test-compat.R` | 0 | 0 | 0 | 0 | 261 | 261 |
| `test-custom-family.R` | 0 | 0 | 0 | 0 | 41 | 41 |
| `test-message-uniqueness.R` | 0 | 0 | 0 | 0 | 6 | 6 |
| `test-bracket-access.R` | 0 | 0 | 0 | 0 | 7 | 7 |

The 49 added expectations are all in `test-compat-register.R`. Every other
named file matches the review exactly, message uniqueness included, so none
of the new refusal texts collides with an existing one.

**Full core suite, one file per process**

* **105 test files on disk, 105 result lines, 105 distinct names.** Set
  difference against `ls tests/testthat/test-*.R` empty in BOTH directions.
  No name twice. No `NO-SUMMARY`. No non-zero return code.
* **fail=0 error=0 warn=0**, skip=21, **pass=5524**. The review's baseline
  was 5475, and 5524 - 5475 = 49, exactly the expectations added, so
  nothing else moved.
* `test-perf.R` passed (3 pass) although it ran alongside the `R CMD check`
  test phase. The machine has 16 logical processors and was at 6 percent
  when the two jobs were started, which is two sequential R jobs, not the
  four heavy ones that broke it before.

**The four extension suites, one process each, against the punched core in
`csp-lib`**

| package | files on disk | files run | fail | error | skip | pass | review |
|---|---|---|---|---|---|---|---|
| frmtmb.ode | 4 | 4 | 0 | 0 | 0 | 211 | 211 |
| frmtmb.latent | 4 | 4 | 0 | 0 | 0 | 227 | 227 |
| frmtmb.ddm | 14 | 14 | 0 | 0 | 0 | 882 | 882 |
| frmtmb.sample | 10 | 10 | 0 | 0 | 2 | **884** | 881 |

Files run were reconciled against files on disk in both directions for each
package, all four empty. `frmtmb.sample` is +3 on the review, which is the
new expects test, and its skip count is unchanged, so the test ran rather
than skipping. `frmtmb.sample` installs and loads under the new
`frmtmb (>= 0.51.0)` floor against core 0.51.0.

Each extension also loaded ALONE in a fresh process against the punched
core, which the cumulative run hides: ode, latent, ddm and sample all load.

**R CMD check --as-cran**

`_R_CHECK_CRAN_INCOMING_=false`, TinyTeX and pandoc on PATH, tarball
`frmtmb_0.51.0.tar.gz`.

* **`Status: 1 NOTE`, return code 0.**
* The only NOTE is `Skipping checking math rendering: package 'V8'
  unavailable`, under `checking HTML version of manual`. Environmental, and
  the expected one.
* `checking tests ... [12m] OK`, with
  **`[ FAIL 0 | WARN 0 | SKIP 26 | PASS 5509 ]`**. The review recorded
  5460, and 5509 - 5460 = 49 again.
* `checking re-building of vignette outputs ... [172s] OK`.
* Relevant to the changed man pages and the new S3 method:
  `checking S3 generic/method consistency ... OK`,
  `checking for code/documentation mismatches ... OK`, every Rd check OK
  including the usage sections.

**Roxygen.** `roxygen2::roxygenise()` was run after each round of roxygen
edits and confirmed idempotent: the run after the last edit wrote
`frmtmb_register_aterm.Rd`, and an immediately repeated run wrote nothing.

**Timing, measured on a quiet machine** (median, warmed):

| | ms |
|---|---|
| bare `frmtmb_compat_features_tbl()` | 16 |
| `frmtmb_register_compat()`, 40 rules | 19 |
| `frmtmb_register_aterm()`, one term | 16 |
| `compat_unresolved_sides()` over 319 rules | 0.1 |
| `frm_compat()`, one pair | 170 |

The query-time report costs 0.1 ms against a 170 ms call, under a tenth of
a percent, and it is the only cost this round adds to a plain user's path.
`frm_compat()`'s own 170 ms is the uncached resolve the review's closing
question is about, and it is untouched here.

## Scope

`git -C C:/Users/adf44/source/r/frmtmb-wt-core-seams diff --name-only
564e185`, 19 files:

    DESCRIPTION  NAMESPACE  NEWS.md
    R/compat.R  R/families.R  R/frame.R  R/parse.R  R/structure.R
    extensions/frmtmb.sample/DESCRIPTION
    extensions/frmtmb.sample/NEWS.md
    extensions/frmtmb.sample/R/zzz.R
    extensions/frmtmb.sample/tests/testthat/test-sampling-ported.R
    man/frm_compat.Rd  man/frmtmb-extension-api.Rd  man/frmtmb_family.Rd
    man/frmtmb_register_aterm.Rd  man/frmtmb_register_compat.Rd
    tests/testthat/test-compat.R  tests/testthat/test-custom-family.R

Untracked, 3:

    dev/core-seams-findings.md  dev/review-core-seams.md
    tests/testthat/test-compat-register.R

The review left 13 modified and 3 untracked. This round adds six modified
files: `DESCRIPTION` and `NAMESPACE` (the version bump and the new S3
method), `man/frm_compat.Rd`, and the three `extensions/frmtmb.sample`
files that were not already in the lane's set.
