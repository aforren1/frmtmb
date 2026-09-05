# Review: wt-core-seams

Reviewer lane. Worktree `C:\Users\adf44\source\r\frmtmb-wt-core-seams`,
branch `wt-core-seams`, uncommitted, from `564e185`. Main checkout
`C:\Users\adf44\source\r\frmtmb` verified clean at `5dfdd84` at the
start and left untouched.

Private library:
`...\scratchpad\rcs-lib` (lane core + four extensions).
Pristine `564e185` core installed separately into `...\scratchpad\rcs-lib-pristine`,
extracted with `git archive` so the repository is not touched.

Verdict, punch list and the `expects =` question are at the end.

---

## 1. Planted misspellings and the sibling refusals

### 1a. Refusal texts, called directly

Fresh session, lane core only, no extension loaded. Every text below is
copied from the run, not from the findings document.

| probe | outcome |
|---|---|
| `wiener x mixture()` | refused, names the rule, the spelling, guesses `mixture` |
| `wiener x dec()` | refused, no guess |
| `hmn x frm_sample` with `expects = c("hmm","lca")` | refused, no guess |
| `gaussian x kind:fmaily` | refused, lists the nine kinds |
| `gaussian x group:cdfs` | refused, lists the 24 groups |
| `gaussian x cens()` status `refuse` | refused, lists the five statuses |
| `wiener x wiener` | refused as a self-pair |
| `features = c(wiener = "familly")` | refused, bad kind |
| `features = c(gaussian = "aterm")` | refused, second kind for one name |
| `rules = <data frame>` | refused, wants a function |
| `expects = 1:2` | refused, wants character |

The near-miss and no-guess texts are byte-identical to the two quoted
in `dev/core-seams-findings.md`, including the trailing
`frm_compat_features() lists the vocabulary.` The four sibling texts
match their quoted forms.

**The atomicity claim holds and is the part worth having.** After all
eleven refusals in one session,
`frmtmb:::frmtmb_compat_contrib$features` and `$rules` are both
length 0, and neither `wiener` nor `frm_sample` is in the vocabulary.
A refused registration leaves nothing behind.

**Why `dec()` gets no guess, confirmed mechanically.** `dec()` is two
edits from `se()` and from `mmc()` under `utils::adist()`, so the
2-edit cap alone would NOT have suppressed the guess. What suppresses
it is `agrep(max.distance = 0.2)`: 0.2 of a 5-character pattern allows
one edit, and both candidates need two. The `adist` cap is the second
filter, not the first. The findings document describes this correctly.

**Vocabulary as measured**: 107 core features, nine kinds, 24 groups,
five statuses. `hmm` and `lca` are NOT core features, which is the
premise of seam 1's second instance. `mixture` is a feature and
`mixture()` is not, which is the premise of the first.

### 1b. The same, through a real `.onLoad()`

Scratch copies of the two extensions under
`...\scratchpad\rcs-work\rcs-plant*`, installed with
`--no-test-load` into their own libraries so the refusal is observed at
`library()` rather than at install. The lane's tree is untouched.

| planted | result |
|---|---|
| ddm, `r("wiener", "mixture()", ...)` | `.onLoad` fails; text identical to 1a, guess `mixture` |
| sample, `r("hmn", "frm_sample", ...)` with `expects` left as-is | `.onLoad` fails, no guess |
| ddm with `"dec()" = "aterm"` REMOVED from `features` | `.onLoad` fails on `wiener x dec()`, no guess |

The third is the historical `frmtmb.ddm` case reproduced exactly, and it
is also the proof that the documented ordering caveat is load-fatal
rather than cosmetic: `ddm_compat_rules` names `dec()` and
`frmtmb_register_aterm("dec")` runs AFTER the registration, so the
explicit `features` entry is the only thing keeping the package
loadable. Reordering the two calls (fourth scratch copy) makes the
explicit entry unnecessary and the package loads. Both orders are
correct; the ordering rule is documented in `R/parse.R` and
`R/compat.R` roxygen and in `man/frmtmb_register_aterm.Rd`.

**One thing the near-miss search does not do.** `compat_near_feature()`
is called with `vocab`, not `c(vocab, expects)`, so a typo of a name
declared in `expects` can never be guessed at. `hmn` happens to have
nothing within two edits in the core vocabulary either, so this case
does not expose it, but a package declaring `expects = "gaussian2"` and
writing `gaussain2` gets the bare refusal. Minor, and arguably right
(a guess pointing at a feature that is not in the session is its own
confusion), but it is an unstated asymmetry. `R/compat.R`, the
`compat_near_feature(pat, vocab)` call inside `compat_check_rules()`.

---

## 2. The second instance, and what `expects =` produces

Four sessions, each a fresh process. Case A and B run the pristine
`564e185` core and the pristine `564e185` `frmtmb.sample`; C and D run
the lane's.

| session | `hmm` in vocabulary | `hmm x frm_sample` row | `frm_compat("hmm","frm_sample")` |
|---|---|---|---|
| A. pristine core, sample, NO latent | FALSE | **ROW ABSENT** | `Unknown feature: 'hmm'` |
| B. pristine core, sample + latent | TRUE | `works` | `works` |
| C. lane core, sample (`expects=`), NO latent | FALSE | ROW ABSENT | `Unknown feature: 'hmm'` |
| D. lane core, sample + latent | TRUE | `works` | `works` |

**A confirms the defect as described.** The registration was accepted,
`frm_sample` reached the vocabulary, and the two rules naming
`frmtmb.latent`'s features were dropped without a word. Nothing in the
session distinguishes that from nobody having said anything.

**C is the lane's choice: the row is absent but the reference is
declared.** It is the right choice. Adding `hmm` to the vocabulary so
the row could exist would put a feature in the table that nothing in
the session implements, and `frm_compat("hmm", "cens()")` would then
answer `untested` about a family that is not loaded, which is a worse
lie than the row being missing. Absent-and-declared says the true thing:
there is no `hmm` here, so there is no pair.

**The edit outside the ownership list was necessary, and this is the
evidence.** Two more sessions:

* Pristine `frmtmb.sample` (no `expects =`) against the LANE core:
  `.onLoad` **fails**. Shipping the core change alone breaks a
  first-party extension at load. The lane's justification stands.
* Lane `frmtmb.sample` (`expects =`) against the PRISTINE `564e185`
  core: `.onLoad` fails with
  `unused argument (expects = c("hmm", "lca"))`.
  See section 7; this is a shippable break, not a theoretical one.

---

## 3. Auto-registration of the addition-term feature

`frmtmb.ddm` loaded against the lane's core, as shipped (explicit
`"dec()" = "aterm"`, rules registered BEFORE the term) and in the
reordered scratch variant. Identical answers:

* `dec()` appears in the vocabulary **exactly once** (`sum(name == "dec()") == 1`),
  of kind `aterm`, in both orderings.
* ddm's explicit `features` entry does not error: the same-kind repeat
  is a no-op, and the vocabulary still holds one row afterwards.
* A second registration under a different kind refuses:
  `... registers 'dec()' as a 'special' feature, and the registry
  already holds that name under the kind 'aterm'. One display name
  carries one kind.`
* Calling `frmtmb_register_aterm("dec")` a second time is also a no-op
  and does not add a row.
* `frm_compat("dec()", "cens()")` and `frm_compat("wiener", "dec()")`
  both answer instead of refusing, which is the gap the decision closes.

In the as-shipped ordering the contributed feature blocks number 1 (the
explicit entry absorbed the auto one); in the reordered variant they
number 2 but still yield one row, which is what the `!duplicated(name)`
backstop in `frmtmb_compat_features_tbl()` is for.

---

## 4. Any-of groups

Ten probes, each a real `frm()` call on a 60-row Poisson custom family,
run against the lane's core and (for the all-of cases) the pristine one.

| probe | declaration | supplied | result |
|---|---|---|---|
| A1 | `c("vint1","vint2")` | nothing | refused, both named |
| A2 | `c("vint1","vint2")` | `vint(size)` | refused, `vint2` only |
| A3 | `c("vint1")` | `vint(size)` | fits |
| B0 | `list("vint1")` | nothing | refused, reads as the plain form |
| B1 | `list(c("vint1","vreal1"))` | `vint(size)` | fits |
| B2 | `list(c("vint1","vreal1"))` | `vreal(z)` | fits |
| B3 | `list(c("dec","vint1"))` | nothing | refused, "one of ... or ..." |
| C1 | `list(c("dec","vint1"), "vreal1")` | nothing | refused, plain first |
| C2 | same | `vint(size)` | refused, `vreal1` only, no "one of" |
| C3 | same | `vint(size) + vreal(z)` | fits |
| C4 | `list("vreal1", c("dec","vint1"))` | nothing | identical text to C1 |

C4 is worth having: the plain-first ordering is imposed by the message
builder, not by the order of declaration, so two families declaring the
same requirement in different orders produce the same sentence.

**The all-of refusal is byte-identical.** A1, A2 and A3 captured from
both cores and compared: `diff` empty, and both files hash to
`cd935b6b45ac62d9a2ccccd898bb5e56`. The claim is exact, not approximate.

### Two constructor behavior changes NEWS does not mention

| `required_aterms =` | pristine | lane |
|---|---|---|
| `c("vint1", NA)` | ACCEPTED | refused, "missing or empty" |
| `list()` | refused | ACCEPTED (no requirements) |
| `list("vint1")` | refused | accepted, plain requirement |
| `data.frame(a="vint1")` | refused | refused (the `!is.object` guard works) |

Both are improvements and neither is likely to bite, but the `NA`
tightening changes what an existing family may declare and is not in the
NEWS entry. Punch list, low.

---

## 5. Adoption dry run in `frmtmb.ddm`

Scratch copy at `...\scratchpad\rcs-work\rcs-adopt\frmtmb.ddm`. Three
edits, all in `R/wiener-family.R`:

1. `required_aterms = list(c("dec", "vint1")),` added to the
   `frmtmb::custom_family()` call.
2. The `if (is.null(up))` branch of `ddm_check_response()` deleted,
   13 lines.
3. The roxygen paragraph calling this "the one hand-rolled check left in
   this package" retired, 7 lines replaced by 4.

Net **-15 lines** in `wiener-family.R`. Installed and exercised:

| formula | before | after |
|---|---|---|
| `rt ~ 1` | `wiener: the decision indicator is missing. ...` | `wiener: the density needs one of` `dec` `or` `vint1` `, ... Write the addition term: rt | dec(<column>) ~ ...` |
| `rt | dec(upper) ~ 1` | fits | fits |
| `rt | vint(upper) ~ 1` | fits | fits |

`wiener()[["required_aterms"]]` is now `List of 1 $ : chr [1:2] "dec"
"vint1"`, and the refusal moves out of `valid_y` (after the frame is
built) into `assemble_frame()` (before it).

### Is the lane's deletion estimate right?

Yes, and for `gddm()` it is right to the line.

* **`wiener()`**: the lane names the check but gives no count. Measured:
  the whole `is.null(up)` branch goes, `ddm_indicator()` and the 0/1
  validity check stay. Correct.
* **`gddm()`**: the lane says "the `is.null(up)` branch of
  `gd_check_response()` (`R/gddm.R` around line 1086, about twelve
  lines including the multi-line refusal)". The branch is at
  **`R/gddm.R:1085`** and spans 16 lines, of which 4 are the comment
  explaining why `required_aterms` could not express it. Twelve lines of
  code. The estimate is accurate and the line number is off by one.
* **`lba()`**: the lane says it gains nothing. Confirmed,
  `R/lba.R:268` already reads `required_aterms = "vint1"`.
* The NEWS claim: **`extensions/frmtmb.ddm/NEWS.md:268`**, "One
  hand-rolled check remains, and is not `required_aterms`'s fault",
  exists at exactly that line and does become false. The lane is right
  to flag it.

### Refusal texts that ddm tests pin, and which ones change

Four expectations across two files break, all on the wiener adoption
alone. The lane does not enumerate them.

| file:line | pins | after adoption |
|---|---|---|
| `tests/testthat/test-family.R:112` | `"decision indicator is missing"` | **BREAKS** |
| `tests/testthat/test-family.R:115` | `"dec\\(decision\\)"` | **BREAKS** |
| `tests/testthat/test-defects.R:157` | `"decision indicator is missing"` | **BREAKS** |
| `tests/testthat/test-defects.R:159` | `"vint\\(upper\\)"` | **BREAKS** |
| `tests/testthat/test-gddm-family.R:158` | `"decision indicator is missing"` | unaffected; that is `gddm()`'s own check, which the wiener-only adoption does not touch |

The framework's replacement text contains none of those three strings.
It names `dec` and `vint1` (the aterm VALUES) rather than
`dec(decision)` and `vint(upper)` (the spellings), so the two
expectations that check the message points at a usable spelling lose
what they were checking: the new sentence writes `rt | dec(<column>) ~
...` and never mentions `vint` at all. That is the real cost of
adoption, and it is a fair trade only because the alternative spelling
is named in the `one of ... or ...` clause.

Two comment blocks also become false and should go with the
expectations: `test-family.R:109-111` ("frmtmb has no way for a family
to declare a required addition term") and `test-defects.R:154-156`.

**A separate pin that adoption of `gddm()` would break, which the lane
does not mention at all:** `tests/testthat/test-gddm-family.R:32-33`
asserts `expect_identical(gddm()[["required_aterms"]], character(0))`
and `expect_setequal(gddm(drift = ...)[["required_aterms"]], ...)`.
Adding the boundary group turns `req` into a list, so `expect_identical`
against `character(0)` and `expect_setequal` over a list both fail.
Whoever adopts `gddm()` owns those two lines as well.

---

## 6. Why `R/parse.R` and `R/structure.R` changed

Read hunk by hunk. Every hunk belongs to a seam.

* **`R/parse.R`, two hunks, both pure insertions, 14 lines total.**
  `@@ -45,0 +46,9 @@` is a roxygen paragraph on
  `frmtmb_register_aterm()` stating that the term joins the
  compatibility vocabulary as `"<name>()"` of kind `"aterm"`, and that
  the term must be registered before rules that name it.
  `@@ -91,0 +101,5 @@` is the `compat_new_aterm_feature(name)` call
  itself plus its comment. Both are seam 1's auto-registration decision;
  nothing else in the file moved.

* **`R/structure.R`, one hunk, 27 lines, all roxygen.** A new
  "What the registries hold you to" section on the
  `frmtmb-extension-api` page, one bullet per registry: what
  `frmtmb_register_compat()` refuses, what `frmtmb_register_aterm()`
  adds, and that `required_aterms` now takes choices. No code. It is
  documentation for both seams on the page an extension author reads,
  which is where it belongs.

* **`R/frame.R`, four hunks, all inside `assemble_frame()` lines
  1150-1176**: the `required_aterm_groups()` call, the plain-first
  ordering, the `needs` and `spell` vectors, and the two changed lines
  of the `stop()`. Seam 2 only.

### No overlap with the `mo` lane

Checked against the sibling worktree
`C:\Users\adf44\source\r\frmtmb-wt-mo-terms`, also from `564e185`.

* `mo_terms_in_brms_order()` does not exist in this lane's tree at all.
  The sibling defines it at **`R/frame.R:131`**, not in `parse.R`.
* Sibling `R/parse.R` hunks: `@@ -878 +878,3 @@`, one line inside
  `parse_linpred()`'s `mo` branch. This lane's parse.R hunks are at 46
  and 101, inside `frmtmb_register_aterm()`. Roughly 780 lines and two
  functions apart.
* Sibling `R/frame.R` hunks: `@@ -118,0 +119,33 @@` (the new function)
  and `@@ -1865,2 +1898,9 @@`, `@@ -1890,6 +1930,5 @@`. This lane's are
  1150-1176. Nearest approach is about 700 lines.

The only file both lanes edit at a shared point is `NEWS.md`, both
adding a development-version block at the top. That is the one merge
conflict to expect, and it is trivial.

---

## 7. The floor consequence, exactly

Measured, not reasoned. `extensions/frmtmb.sample` built from the
lane's source, loaded against the pristine `564e185` core in a library
containing only that one package:

    .onLoad failed in loadNamespace() for 'frmtmb.sample', details:
      call: frmtmb_register_compat(features = c(frm_sample = "method"),
            rules = sample_compat_rules,
      error: unused argument (expects = c("hmm", "lca"))

Its DESCRIPTION says `frmtmb (>= 0.46.0)`. Released core 0.50.0
satisfies that floor, so today's declaration permits an installation
that cannot load. The failure is not a refusal the registry wrote; it is
a raw `unused argument`, which tells the user nothing.

**What must change at consolidation.**

1. **`DESCRIPTION` (core), `Version: 0.50.0`.** Must rise. `expects =`
   is new API, the strict check is a behavior change, and
   `required_aterms` accepts a new type. Call it `0.51.0`; the NEWS
   block is already written as `# frmtmb (development version)`.
2. **`extensions/frmtmb.sample/DESCRIPTION`, `frmtmb (>= 0.46.0)` ->
   `frmtmb (>= 0.51.0)`.** MANDATORY, and it is the only mandatory one.
   Verified as a hard load failure above.

**What must NOT change, checked rather than assumed.** Each extension
was loaded ALONE in a fresh process against the lane's core, which is
the check a cumulative run hides:

| extension | floor today | loads alone on lane core | floor must rise? |
|---|---|---|---|
| frmtmb.ode | `>= 0.46.0` | yes | no |
| frmtmb.latent | `>= 0.48.0` | yes | no |
| frmtmb.ddm | `>= 0.49.0` | yes | **no, as it stands** |
| frmtmb.sample | `>= 0.46.0` | yes | **YES -> 0.51.0** |

`frmtmb.ddm` is the interesting one and the answer is no. It declares
`"dec()" = "aterm"` in `features =` itself, so its rules resolve on the
old core (which checked nothing) and on the new one (where the explicit
entry is a no-op against the auto-registered feature). The
auto-registration is a convenience it does not depend on. Nothing in
this repository yet relies on a term reaching the vocabulary without
being declared, so the auto-registration adds no floor to anything.

**But `frmtmb.ddm`'s floor rises the moment section 5's adoption lands**,
because `required_aterms = list(...)` is refused by the pristine
constructor with `names the addition terms ... as a character vector,
not a list of length 1`. If adoption and this lane ship together,
`extensions/frmtmb.ddm/DESCRIPTION` goes to `frmtmb (>= 0.51.0)` as
well. Keeping the two apart buys nothing, so the release owner should
plan on both.

**A third thing the release owner owns, which the findings document does
not raise.** The core is now strict about contributed rules, and the
strictness is enforced at `.onLoad()` of the CONTRIBUTOR. Any
third-party extension outside this repository that today registers a
rule naming a feature it does not supply will stop loading on upgrade,
with no deprecation window. Within the repository there is exactly one
such package and the lane fixed it; outside it, nobody can know. That is
an argument for the floor bump being a minor-version bump with the
behavior called out in NEWS (it is), not for softening the check.

## 8. Suites

Scope taken from `git diff --name-only 564e185` plus untracked, in the
worktree only. Main was clean at `5dfdd84` at the start and is clean and
unmoved at the end.

### Full core suite, one file per process

`pkgload::load_all()` under `NOT_CRAN=true`, one `Rscript` per file,
log truncated once at the start and reconciled against the disk by name
in both directions afterwards, as the findings document warns to do.

* **105 test files on disk, 105 result lines, 105 distinct names.**
* Set difference against `ls tests/testthat/test-*.R` empty in BOTH
  directions. No name appears twice. No `NO-SUMMARY`, no non-zero return
  code, no crash line.
* **fail=0 error=0 warn=0**, skip=21, pass=5475.

The five files named in the brief, each in its own process:

| file | result |
|---|---|
| `test-compat-register.R` | 38 pass, 0 fail |
| `test-compat.R` | 261 pass, 0 fail |
| `test-custom-family.R` | 41 pass, 0 fail |
| `test-message-uniqueness.R` | 6 pass, 0 fail |
| `test-bracket-access.R` | 7 pass, 0 fail |

`test-perf.R` also passed (3 pass) despite two other heavy R jobs on the
machine, which is the file the findings document flags as wall-clock
sensitive.

### The four extension suites, one process each

`testthat::test_local(load_package = "installed")` against the lane's
core installed in `rcs-lib`. Files run reconciled against files on disk.

| package | files on disk | files run | fail | error | skip | pass |
|---|---|---|---|---|---|---|
| frmtmb.ode | 4 | 4 | 0 | 0 | 0 | 211 |
| frmtmb.latent | 4 | 4 | 0 | 0 | 0 | 227 |
| frmtmb.ddm | 14 | 14 | 0 | 0 | 0 | 882 |
| frmtmb.sample | 10 | 10 | 0 | 0 | 2 | 881 |

Matches the findings document's table on every count that matters. My
`frmtmb.ddm` run reports 0 skips where the lane reports 1, which is the
machine having `RWiener` installed and not a disagreement about the
package.

The `frmtmb.ddm` number is also the control for section 5: **882 passes
here against 878 passes and 2 errors in the adopted build.** The
difference is exactly 4 expectations, which is exactly the four this
review predicted from reading the tests, in exactly the two files
predicted. Nothing else in `frmtmb.ddm` moves on adoption.

Each extension was additionally loaded ALONE in a fresh process against
the lane's core, which is the check a cumulative run hides. All four
load. (See section 2 for the two cross-version loads that do NOT.)

### `R CMD check --as-cran`

`_R_CHECK_CRAN_INCOMING_=false`, TinyTeX and pandoc on PATH.

* **`Status: 1 NOTE`, return code 0.**
* The only NOTE is under `checking HTML version of manual`:
  `Skipping checking math rendering: package 'V8' unavailable`. That is
  the machine, not the package. Environmental, as claimed.
* `checking tests ... [17m] OK`, with
  **`[ FAIL 0 | WARN 0 | SKIP 26 | PASS 5460 ]`** -- identical to the
  findings document's figures, reproduced independently.
* `checking re-building of vignette outputs ... [265s] OK`; PDF and HTML
  manuals built.
* Relevant to the four changed man pages:
  `checking for code/documentation mismatches ... OK` and
  `checking Rd \usage sections ... OK`, so the regenerated Rd matches
  the new `frmtmb_register_compat(features, rules, expects)` signature
  and the widened `required_aterms` documentation. Verified by hand as
  well: `man/frmtmb_register_compat.Rd` carries `expects` in both
  `\usage` and `\arguments`.

One correction to the reproduction notes, not to the package: the
findings document says to put TinyTeX on PATH before checking. Pandoc is
also required. Without it `R CMD build` fails while re-building all
seven vignettes, before `check` is ever reached. Pandoc is at
`C:\Users\adf44\AppData\Local\Programs\Quarto\bin\tools` on this
machine (and under the RStudio install).

### Scope and hygiene

* Scope from `git diff --name-only 564e185` plus untracked: 13 modified
  files, 2 untracked (the lane's findings and the new test file).
* Main checkout `C:\Users\adf44\source\r\frmtmb`: `5dfdd84`, `git status`
  empty, at the start of this review and again at the end. Never written
  to, never committed to.
* Worktree HEAD still `564e185`. Nothing committed.
* All scratch work under `...\scratchpad\rcs-*`. No file created by
  another lane was read destructively, truncated or deleted; no process
  was stopped.

---

## Code read: six things the lane's own verification would not catch

The last of them, (f), is the one that blocks.

### (a) A confidently wrong suggestion, which is the case the design set out to avoid

`compat_near_feature()` guesses `us` for the misspelling `se`:

    se           -> Did you mean 'us'?

`se()` IS in the vocabulary. The paren-exact branch is

    known == sub("\\(\\)$", "", bad) | paste0(known, "()") == bad

and both halves test the same direction: `bad` carrying parentheses the
vocabulary entry does not. The reverse, `bad` MISSING the parentheses,
falls through to `agrep()`. `agrep()` does find the right entry in most
such cases, because it matches the pattern as an approximate SUBSTRING
(verified: `cens`, `trials` and `weights` all guess correctly, and
`agrep("cens", "cens()", max.distance = 0.2)` matches despite an edit
distance of 2). But the `adist` ranking that follows prefers `us` at one
edit over `se()` at two. `se` is an addition term; `us` is an
unstructured covariance. Different kind, unrelated meaning.

This is the exact failure the findings document argues the tight
threshold prevents ("which is why `dec()` gets no guess rather than
being pointed at `se()`"). The threshold does prevent it for `dec()`. It
does not prevent it here, because the miss happens to be one edit from a
real entry.

One-line fix, `R/compat.R`, in `compat_near_feature()`: make the paren
branch bidirectional by adding `| known == paste0(bad, "()")` to the
`hit` test. `se` then resolves exactly and never reaches `agrep()`.
The other probes are undisturbed, and `dec()` and `hmn` stay unguessed.

### (b) A duplicate display name inside ONE `features =` call is dropped in silence

    frmtmb_register_compat(features = c(zzq = "family", zzq = "aterm"))
    # ACCEPTED. zzq is registered as "family"; the "aterm" declaration vanishes.

`compat_new_features()` refuses one name under two kinds when the second
comes from a DIFFERENT registration ("One display name carries one
kind", verified in section 3) but drops it without a word when both come
from the same vector, at

    features <- features[!duplicated(names(features))]

That is seam 1's own defect, one scope inward: a declaration accepted,
doing nothing, saying nothing. The line should refuse when duplicated
names carry different kinds and stay silent only when they agree. Punch
list, medium. It is small, but this is the lane that cannot have it.

### (c) A `rules` builder that throws now makes the package unloadable, with no context

Pristine called `rules()` lazily, at the first `frm_compat()`. The lane
calls it at registration, so a builder that errors takes down
`.onLoad()`. The error propagates raw:

    frmtmb_register_compat(rules = function() stop("boom"))
    #> boom

The registrant sees `boom` out of `loadNamespace()` with nothing naming
`frmtmb_register_compat(rules =)`. Every other refusal on this path is
carefully labeled. Wrapping the `rules()` call in a `tryCatch` that
prefixes the argument would match. Punch list, low.

### (d) `frmtmb_register_aterm()` rebuilds the vocabulary on every call

`compat_new_aterm_feature()` calls `compat_new_features()` without
passing `known`, so it takes the default and rebuilds the feature table.
Measured on this machine: a bare vocabulary build is 26 ms, a 40-rule
`frmtmb_register_compat()` is 48 ms (1.8 builds, which confirms the
findings document's "1.9"), and `frmtmb_register_aterm()` is 22 ms,
essentially one build. An extension registering four terms pays four.
Small, and it is the same fix the lane already applied one function
over. Punch list, low.

### (e) The near-match pool excludes `expects`

Noted in 1b. `compat_near_feature(pat, vocab)`, not
`c(vocab, expects)`. Defensible, undocumented.

---

### (f) THE ONE THAT BLOCKS: `frmtmb_register_aterm()` reports failure and registers anyway

Found by asking what the auto-registration does when the display name
`"<name>()"` already belongs to a feature of another kind. Thirteen
vocabulary entries are written with parentheses and are NOT of kind
`aterm`:

    s() t2() mo() mi_pred() gp_pred() cs_pred()   (special)
    ar() ma() arma() cosy() unstr()               (autocor)
    mm() mmc()                                    (grammar)

None of them is in `core_aterms`, so `frmtmb_register_aterm()` accepted
all thirteen names before this lane. Now:

    frmtmb_register_aterm("s")
    #> Error: frmtmb_register_aterm() registers 's()' as a 'aterm' feature,
    #>   and the registry already holds that name under the kind 'special'.
    #>   One display name carries one kind.

That refusal may well be right on the merits -- `y | s(col) ~ x` beside
a smooth `s(x)` is a genuine grammar collision. **The defect is that the
refusal is not atomic.** In `R/parse.R`, `frmtmb_register_aterm()` does:

    frmtmb_aterm_registry$reg[[name]] <- list(name = name, ...)   # line 101
    ...
    compat_new_aterm_feature(name)                                # line 105

The parser registration is committed first, and the compat registration
throws after it. Measured in one session on the lane's core:

    before: 's' in aterm registry: FALSE
    register_aterm('s') -> Error: ... already holds that name ...
    AFTER THE REFUSAL: 's' in aterm registry: TRUE
    AFTER THE REFUSAL: 's()' in compat vocabulary: TRUE  (kind: special)
    frm(bf(y | s(col) ~ x), family = gaussian(), data = d)
    #> FIT OK frmtmb_fit

The call said no and did yes. `frm()` now accepts `y | s(col) ~ x` as an
addition term, and `frm_compat()` describes `s()` as a smooth. That is
precisely the state the seam exists to make impossible: a term the
parser accepts and the table cannot describe -- reached, this time, by
the function that was supposed to close it.

It is also the property the lane was careful about one function over.
`R/compat.R` carries the comment "Nothing is committed until every rule
has resolved, so a refused registration leaves the vocabulary as it
found it rather than half filled". `frmtmb_register_aterm()` has no such
ordering.

In an extension this surfaces as an unloadable package whose
`.onLoad()` has already mutated a registry inside frmtmb's namespace.
Any later registration in that same `.onLoad()`, and any caller that
traps the error, runs on with the split state.

**Fix, one line moved, `R/parse.R` lines 99-105:** hoist
`compat_new_aterm_feature(name)` ABOVE the
`frmtmb_aterm_registry$reg[[name]] <- ...` assignment. Nothing after it
can throw, so that makes the seam atomic in the same sense
`frmtmb_register_compat()` already is.

Worth doing at the same time, since the message is what a registrant
will read: the collision is a formula-grammar fact, not a table fact,
and `... and the registry already holds that name under the kind
'special'` does not tell an author that `s()` is a smooth. A sentence
naming the conflict in `frmtmb_register_aterm()`'s own terms would.
Also `as a 'aterm'` should be `as an 'aterm'`, `R/compat.R` in
`compat_new_features()`.

---

## Non-regression: the core compatibility surface is unchanged

Both cores, no extension loaded, `frm_compat_features()` and
`frm_compat()` dumped and compared with `identical()`:

* 107 features, 5041 resolved pairs, on both.
* `identical(features)` TRUE. `identical(table)` TRUE.

The `!duplicated(out$name)` de-duplication and the
`rownames(out) <- NULL` added to `frmtmb_compat_features_tbl()` change
nothing about what the core alone reports. Whatever the strictness costs
a contributor, it costs a plain user nothing.

---

## VERDICT: GO WITH FIXES

> SUPERSEDED by the punch re-check below (2026-09-05): every row
> in this punch list is closed or withdrawn, and the updated
> verdict is GO. Kept as the record of what was found.

Both seams are real, both diagnoses are correct, and both fixes do what
the findings document says they do. Every quoted refusal text reproduced
byte for byte, the all-of message is provably unchanged
(`md5 cd935b6b...` on both trees), the atomicity claim holds under
eleven consecutive refusals, and the second instance reproduces exactly
as described on the pristine tree. The `frmtmb.sample` edit outside the
ownership list was necessary and is justified by evidence I reproduced:
without it, the core change makes a first-party extension unloadable.

One defect blocks, and it is in the part of the lane that had the least
scrutiny because it was framed as a decision rather than a change.

### Punch list

| # | severity | file:line | what |
|---|---|---|---|
| 1 | **blocking** | `R/parse.R:99-105` | `frmtmb_register_aterm()` commits to the parser registry BEFORE `compat_new_aterm_feature()` can refuse. For 13 names (`s`, `t2`, `mo`, `mi_pred`, `gp_pred`, `cs_pred`, `ar`, `ma`, `arma`, `cosy`, `unstr`, `mm`, `mmc`) the call now throws and registers anyway: `register_aterm("s")` errors, then `frm(bf(y \| s(col) ~ x))` FITS. Fix: move line 105 above line 99. |
| 2 | medium | `R/compat.R:224` | `features <- features[!duplicated(names(features))]` silently drops a same-call duplicate carrying a DIFFERENT kind. `c(zzq = "family", zzq = "aterm")` is accepted and half applied. The cross-call case is refused at line 229; make this one match. |
| 3 | medium | `R/compat.R:350-351` | `compat_near_feature()` answers `Did you mean 'us'?` for the misspelling `se`, when `se()` is in the vocabulary. Both halves of the `hit` test check the same paren direction. Fix: add `\| known == paste0(bad, "()")`. |
| 4 | low | `R/compat.R:166` | a `rules` builder that throws now takes down `.onLoad()` with a raw, unlabeled error. Wrap `rules()` so the message names `frmtmb_register_compat(rules =)`. |
| 5 | low | `R/compat.R:247` | `compat_new_aterm_feature()` does not take `known`, so every `frmtmb_register_aterm()` rebuilds the 107-row vocabulary (22 ms measured). Same fix the lane already applied to `frmtmb_register_compat()`. |
| 6 | low | `R/compat.R:229` | `registers 'x' as a 'aterm' feature` -> `as an`. |
| 7 | low | `NEWS.md` | the `required_aterms = c("x", NA)` tightening (accepted before, refused now) is not mentioned. |
| 8 | release | `DESCRIPTION`, `extensions/frmtmb.sample/DESCRIPTION` | core `0.50.0` -> `0.51.0`; sample floor `>= 0.46.0` -> `>= 0.51.0`. Verified hard load failure otherwise. See section 7. |
| 9 | trivial | `dev/core-seams-findings.md` | the prerequisite note says to put TinyTeX on PATH; pandoc is also required, and without it `R CMD build` fails on all seven vignettes before `check` starts. |

None of 2-7 is a reason to hold the lane; 1 is, because it puts the
registry into the exact state the lane exists to make unreachable.

### Every edit I made

One file, created by me, in the worktree:

* `C:\Users\adf44\source\r\frmtmb-wt-core-seams\dev\review-core-seams.md`
  (this document).

No other file in the worktree was touched. `git status` in
`frmtmb-wt-core-seams` shows the lane's own 13 modified and 2 untracked
files exactly as I found them, plus this one. The main checkout
`C:\Users\adf44\source\r\frmtmb` was verified clean at `5dfdd84` at the
start and again at the end, and never written to. Nothing was committed.
All scratch work is under `...\scratchpad\rcs-*`.

### Is `expects =` the right mechanism?

It is a symptom: the registry already re-resolves the whole table on
every `frm_compat()` call with no cache, so the `hmm x frm_sample` rule
goes live the moment `frmtmb.latent` loads without anyone re-registering
anything -- which means `expects` buys nothing at run time and exists
only to suppress a registration-time refusal, is never checked against
anything then or later, and so re-opens the silent-dangling-rule hole
one argument wide; the honest fix is to let the resolver report
unresolved rule sides at query time, when the session's package set is
actually known, and `expects` should be read as the note that says so.

---

# Punch re-check, 2026-09-05

Second pass by the reviewer lane, against the worker's record at
`dev/core-seams-findings.md:346`. Fresh install of the punched core and
all four extensions into `...\scratchpad\rcs-lib2` (core reports
`0.51.0`); the pre-punch install in `rcs-lib` kept for A/B measurement.
Main checkout verified clean at `5dfdd84` and unmoved. Worktree HEAD
still `564e185`, nothing committed. 19 tracked modified + 3 untracked,
as the worker states (the third untracked is this document).

## Punch rows, re-checked

| row | was | now | verified how |
|---|---|---|---|
| 1 | blocking | **CLOSED** | reproduced below |
| 2 | medium | **CLOSED** | same-call duplicate refused |
| 3 | medium | **CLOSED** | `se` -> `se()` |
| 4 | low | **CLOSED** | builder failure names package and argument |
| 5 | low | **WITHDRAWN, my premise was wrong** | measured below |
| 6 | low | **CLOSED** | article correct on both vowel kinds |
| 7 | low | **CLOSED** | NEWS covers the `NA` and `list()` changes |
| 8 | release | **CLOSED** | `0.51.0`, sample floor `>= 0.51.0` |
| 9 | trivial | **CLOSED** | pandoc named, `dev/core-seams-findings.md:288` |

Plus item 6 of the worker's own numbering, the `unresolved` reporting,
which was not on my punch list and is new work. Exercised end to end
below.

### Row 1, the blocking one: reproduced

`R/parse.R:117-119` now calls `compat_new_aterm_feature(name)` BEFORE
`frmtmb_aterm_registry$reg[[name]] <- ...`. Three collisions, one per
colliding kind, each in a fresh session:

| call | refusal | parser registry after | vocabulary after |
|---|---|---|---|
| `frmtmb_register_aterm("s")` | `'s()' already means something else in a formula ... holds it as a 'special' feature` | 0 entries, `s` absent | 107 rows, `s()` still `special` |
| `frmtmb_register_aterm("ar")` | `... as an 'autocor' feature` | 0 entries, `ar` absent | 107 rows, `ar()` still `autocor` |
| `frmtmb_register_aterm("mmc")` | `... as a 'grammar' feature` | 0 entries, `mmc` absent | 107 rows, `mmc()` still `grammar` |

`PARSER REGISTRY UNCHANGED: TRUE` and `VOCABULARY UNCHANGED: TRUE` on
all three. And the behavior that proved the defect is gone:

    frm(bf(y | s(col) ~ x), family = gaussian(), data = d)
    #> REJECTED: Addition term `s()` is not supported (supported:
    #>   weights(), trials(), cens(), trunc(), se(), vint(), vreal(), mi())

which is the pre-lane answer. A non-colliding name still registers
normally: `frmtmb_register_aterm("dec")` is accepted, lands in the
parser registry, and yields exactly one `dec()` vocabulary row.

The new message is also the right message. It is framed as a grammar
fact ("a formula writing `s()` could not say which was meant"), not as
a table fact, which is what I asked for.

### Row 5: WITHDRAWN. The worker is right and I was wrong

I claimed `frmtmb_register_aterm()` paid a vocabulary build it could
avoid, and called it "the same fix the lane already applied one function
over". That inferred a doubling from the missing `known` argument. There
was no doubling: `compat_new_features(known = frmtmb_compat_features_tbl())`
is a lazy default, evaluated once. My 22 ms figure was one build measured
while two other suites had the machine.

Measured properly this round -- registry reset between repetitions so
the vocabulary stays at 107 rows, 200 repetitions per cell so the ~15 ms
Windows clock tick is amortized rather than quantizing each call:

| | pre-punch 0.50.0 | punched 0.51.0 |
|---|---|---|
| bare `frmtmb_compat_features_tbl()` | 15.96 ms | 15.86 ms |
| `frmtmb_register_aterm()` | 17.75 ms | 18.29 ms |
| ratio | 1.11 builds | 1.15 builds |

One build either way, ~16 ms, exactly as the worker says. The plumbing
of `known` at `R/compat.R:336-341` earns its place by keeping the new
collision check from making it two, not by removing a second build that
was never there.

The worker's third point stands as well: four terms still cost four
builds, and that wants a vocabulary cache, which is outside this lane.
Recorded as a residual below rather than a punch row.

### Rows 2, 3, 4, 6: reproduced

**Row 2**, `R/compat.R:298-312`. `c(zzq = "family", zzq = "aterm")` is
now refused: *"registers 'zzq' as an 'aterm' feature, and the same call
already gives that name the kind 'family'"*. The harmless half stayed
silent: `c(zzq = "family", zzq = "family")` is accepted and leaves one
vocabulary row. The cross-call clash still refuses.

**Row 3**, `R/compat.R:461-464`, the paren test made three-way. Twelve
probes, all correct, and the one that was wrong is right:

    se        -> Did you mean 'se()'?     (was: 'us')
    cens      -> Did you mean 'cens()'?
    vint      -> Did you mean 'vint()'?
    trunc     -> Did you mean 'trunc()'?
    mixture() -> Did you mean 'mixture'?
    gausian   -> Did you mean 'gaussian'?
    dec()     -> (no guess)
    hmn       -> (no guess)

Both directions resolve exactly now; `dec()` and `hmn` still get no
guess, so the tightness the design wanted is intact.

**Row 4**, `R/compat.R:198-206` with `compat_registrant()` at `:244-250`.
From a namespace:

> frmtmb_register_compat() called the rules = builder **that
> frmtmb.sample gave it** and it failed, which stops the whole
> registration: ... It said: boom

The package clause fires only when the caller is a namespace, which is
correct -- from the global environment it degrades to "called the rules =
builder and it failed", still naming the argument.

**Row 6**, the article: `an 'aterm'`, `an 'autocor'`, `a 'special'`,
`a 'grammar'`. Correct on both vowel kinds.

## The worker's item 6: `expects =` reporting, end to end

Not one of my punch rows. This is the piece answering my closing
paragraph, and it is the most substantial thing in the round.

### Exercised in one session, in order

| step | class of `frm_compat()` | `unresolved` | `hmm x frm_sample` |
|---|---|---|---|
| core alone | `data.frame` | none | -- |
| + `frmtmb.sample` | `frmtmb_compat/data.frame` | `hmm`, `lca` | ROW ABSENT |
| + `frmtmb.latent` | `data.frame` | none | `works` |

The print method fires on the middle row:

> Unresolved rule sides in this session: 'hmm', 'lca'. A registered rule
> names each, declared through frmtmb_register_compat(expects =) as a
> feature another package supplies. No feature here answers to them, so
> those rules match no pair until the owning package is loaded.

Loading `frmtmb.latent` in the SAME session clears the attribute, drops
the class back to a plain data frame, and resolves the pair to `works`.
`frm_compat("hmm", "frm_sample")` then answers `works`.

**Non-regression confirmed**: the core alone gets exactly the plain
5041-row data frame it got before, no attribute, no class. The
`frmtmb_compat` class exists only in a session that actually holds a
forward reference.

`expects =` naming a feature the session already has is refused
(`R/compat.R:177-186`), tried three ways -- `gaussian` (a family),
`cens()` (an aterm) and `frm_sample` (a method contributed by the
loaded package itself):

> frmtmb_register_compat(expects =) names 'gaussian', and the registry
> already has that feature, so the declaration exempts nothing.

and a genuinely absent name is still admitted, still kept OUT of the
vocabulary, and shows up in `unresolved`.

### Is the `unresolved` attribute the right representation?

**Yes -- a status value could not exist, and the attribute is the only
honest home. But it is reported where users do not look.**

A status value would need a row to sit on, and there is none:
`sum(feature_a == "hmm" | feature_b == "hmm")` is **0**, because `hmm`
is not a feature of the session. Manufacturing the row means putting
`hmm` in the vocabulary, which is exactly the lie -- a feature in the
table that nothing implements -- that this whole lane exists to refuse.
Dangling-ness is a property of the RULE SET, not of any pair, and an
attribute on the table is where a property of the rule set belongs.

The one line of evidence on placement, which is where it falls short:

> `frm_compat("hmm", "frm_sample")` -- the only query a user with this
> problem would ever type -- still fails with
> `Unknown feature: 'hmm'. See frm_compat_features().`, never mentioning
> that a registered rule names `hmm` and is waiting for its package;
> and `frm_compat_features()`, the remedy that error points at, is a
> plain 108-row data frame carrying no `unresolved` attribute and no
> `hmm`, so following the pointer teaches nothing.

Meanwhile the note that would have said it sits at printed line 15,449
of 15,452 for the full table. It is legible on a subsetted query (the
attribute survives `[`, `head()` and `as.data.frame()`, which is good
and deliberate), but the two paths a stuck user actually walks are the
two that stay silent.

So: right representation, wrong reach. Two residual lines, both small.

## Suites, re-run against a fresh install of the punched core

One process per file, `pkgload::load_all()` under `NOT_CRAN=true`, into
the fresh `rcs-lib2` install (core `0.51.0`).

| suite | worker reports | I measure |
|---|---|---|
| `test-compat-register.R` | 87 (was 38) | **87**, 0 fail |
| `test-compat.R` | 261 | **261**, 0 fail |
| `test-custom-family.R` | 41 | **41**, 0 fail |
| `test-message-uniqueness.R` | 6 | **6**, 0 fail |
| `test-bracket-access.R` | 7 | **7**, 0 fail |
| `frmtmb.sample` suite | 884 | **884**, 10 files, 0 fail, 0 error, 2 skip |

Every count matches. The sample figure is my own pre-punch 881 plus the
3 the worker added, as claimed.

Full core suite as well, which the coordinator did not ask for but which
is the number worth checking:

* **105 files on disk, 105 result lines, 105 distinct names**, set
  difference empty in both directions, no duplicate, no `NO-SUMMARY`, no
  non-zero return code.
* **fail=0 error=0 warn=0**, skip=21, **pass=5524**.

5524 against my pre-punch 5475 is +49, which is the worker's additions,
exactly as reported.

I did NOT re-run `R CMD check --as-cran` this round; it was not in the
re-check list. The worker's `Status: 1 NOTE` / `PASS 5509` figures for it
are unverified by me, though the last round's independent run of the
pre-punch tree reproduced its equivalents exactly.

## Residual items

| # | severity | file:line | what |
|---|---|---|---|
| R1 | low-medium | `R/compat.R:1675-1679` | `check_features()` tests only `setdiff(x, ft$name)`, so `frm_compat("hmm", "frm_sample")` still dies with `Unknown feature: 'hmm'` and never says a registered rule is waiting on that name. This is the one call path a stuck user takes. Consulting `compat_unresolved_sides()` here would close it. |
| R2 | low | `R/compat.R`, `frm_compat_features()` | carries no `unresolved` attribute and no print method, so the remedy `See frm_compat_features().` that R1's error points at leads to a plain 108-row frame that says nothing about the pending rule. |
| R3 | low, outside this lane | `R/compat.R`, `frmtmb_compat_features_tbl()` | no cache: ~16 ms per build, paid once per `frm_compat()` call and once per `frmtmb_register_aterm()`. Four terms cost four builds. The worker says this needs a cache outside the lane and is right; recorded so it is not lost. |
| R4 | deferred, not a defect | `extensions/frmtmb.ddm` | the section 5 adoption was correctly NOT done: it is not a punch row, it breaks four pinned expectations (`test-family.R:112,115`, `test-defects.R:157,159`), and it would take ddm's floor to `>= 0.51.0`. Owner is whoever adopts, with section 5 of this document as the worked example. |

R1 and R2 are the same gap seen twice and would be one small change.
Neither blocks: the information exists, is correct, and is reachable --
it just is not offered on the path where it is needed.

## UPDATED VERDICT: GO

The blocking defect is closed, and closed properly: the fix is the
ordering change I asked for, the refusal now reads as the grammar fact
it is, and I reproduced atomicity on three different colliding kinds
with both registries verified unchanged and `frm()` rejecting
`y | s(col) ~ x` again.

All nine punch rows are closed or, in the case of row 5, withdrawn
because the worker showed my premise was wrong and my own A/B
measurement confirms it: `frmtmb_register_aterm()` pays one ~16 ms
vocabulary build, before and after, and always did.

The unrequested item 6 work is the best thing in the round. It answers
the question my first pass ended on -- `expects =` was an annotation
nothing ever reconciled -- by making the reconciliation happen at query
time, where the session's package set is finally known, and by refusing
an `expects =` that exempts nothing. The representation is right for a
structural reason: there is no row to hold a status, and inventing one
would be the lie the lane exists to refuse.

Nothing outstanding rises above low. Ship it, with R1 and R2 as a
follow-up whenever the compat surface is next opened.

### Every edit I made this round

None to the package. One file, appended to by me and by me only:
`C:\Users\adf44\source\r\frmtmb-wt-core-seams\dev\review-core-seams.md`.
Main checkout `C:\Users\adf44\source\r\frmtmb` verified at `5dfdd84`
with an empty `git status` at the start and end of this round, never
written to. Worktree HEAD still `564e185`. Nothing committed. All
scratch under `...\scratchpad\rcs-*`; no file or process belonging to
another lane was touched.
