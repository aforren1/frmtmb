# Addition-term declarations: exclusivity, derived compatibility rows,
# and a dead accessor

Lane `wt-aterms`, worktree `C:/Users/adf44/source/r/frmtmb-wt-aterms`,
branch `wt-aterms` off `main` at frmtmb 0.53.0 / frmtmb.eam 0.4.0.
Written incrementally as the work was measured.

Read first: `dev/structured-family-protocol.md`, section "The
addition-term allow-list (added 2026-09-05)", which records defect 1
under what the allow-list deliberately left out.

## Read this before the corrected numbers: the defect was in the instrument, not the product

Three numbers this lane published did not reproduce, and the punch
review was right about all three. A reader who sees a corrected count
will reasonably assume the shipped code was fixed. **It was not, and it
was never wrong.**

`compat_declared_displacements()` recognizes a derived row by comparing
its note for EXACT equality with `compat_derived_note()`. That is the
shipped predicate, it is what `test-compat.R` asserts against, and it
classifies the one ambiguous row correctly.

The throwaway probe I measured with did something else: it matched the
note PREFIX `"Refused by declaration:"`. `frmtmb.eam` hand-writes an
`rdm` x `dec()` note that opens with those same six words, so the probe
counted that hand-written row as derived, removed it from the baseline,
and reported it as a cell that moved. Measured both ways on the same
registry:

| predicate | `rdm` x `dec()` classified derived |
|---|---|
| exact note equality (shipped) | **0** |
| note prefix (my probe) | 1 |

One miscounted row inflated all three published numbers at once: the
cells that move, the `dec()` share of them, and the eam share. No shipped
behavior changed when the numbers were corrected, and no shipped
behavior needed to.

The lesson is not "check the arithmetic". It is that a measurement
instrument needs the same scrutiny as the thing it measures, and that a
predicate written twice, once for the product and once for the probe,
is a place where the two can disagree silently. The probe should have
called `compat_derived_note()` instead of restating a prefix of it.

## Defect 1: no declaration says two terms are mutually exclusive

### Measured, before any change

`wiener()` declares `accepts_aterms = c("dec", "vint", "weights")` and
`required_aterms = list(c("dec", "vint1"))`. Both spellings are on the
allow-list and the any-of group is satisfied by either, so supplying
BOTH passes every guard. `ddm_indicator()`
(`extensions/frmtmb.eam/R/wiener-family.R:441`) reads `dec` first and
falls back to `vint1`, so the second column is never read.

120 rows, `rt = 0.3 + rexp(n, 3)`, `up = rbinom(n, 1, 0.6)`, seed 2:

| formula | logLik |
|---|---|
| `rt \| dec(up) ~ 1` | -76.0486443897369 |
| `rt \| vint(up) ~ 1` | -76.0486443897369 |
| `rt \| dec(up) + vint(1 - up) ~ 1` | -76.0486443897369 |

The third is the defect at its sharpest: the two columns CONTRADICT each
other, the fit runs, and the log-likelihood is bit-identical to the
`dec()`-only model. Nothing warns.

## Defect 2: untested cells that the declarations already settle

### Audit, before any change

Every core family in the compatibility vocabulary was constructed from
`family_registry`, its `accepted_aterm_names()` read off, and the pair
looked up in `frm_compat()`. 36 families x 8 addition terms, minus the 8
cells `multinomial` could not answer for (its constructor needs `K`):

| table status | declaration | cells |
|---|---|---|
| refused | refuses | 123 |
| untested | refuses | 102 |
| works | accepts | 53 |
| conditional | accepts | 2 |
| (multinomial, unanswerable) | - | 8 |

**Zero disagreements.** No cell where the table says `works` or
`conditional` and the declaration refuses, and none where the table says
`refused` and the declaration accepts.

The 102 derivable cells are: `trials()` against the 32 families that are
not binomial-shaped, and `vint()` and `vreal()` against 35 families
each.

The `poisson` x `cens()` contradiction the previous round reported is
already gone at this head: `fam_poisson()` declares
`accepts_aterms = c("weights", "trunc")`, and the hand-written row
`r("cens()", "poisson", "refused", ...)` agrees with it.

Confirmed that a derived refusal is a FACT and not an inference. With
`gaussian()`:

```
frm(bf(y | vint(k) ~ x), data = d, family = gaussian())
#> Error: gaussian: the addition term `vint()` is not one this family
#>   reads, so writing it would change nothing about the fit. This
#>   family takes `cens()`, `mi()`, `se()`, `trunc()`, `weights()`.
```

`trials()` on `gaussian()` refuses identically.

### What is NOT derivable

Only refusals. A family that ACCEPTS a term has said nothing about
whether the pair works, and `untested` stays the honest answer:
`gddm` x `weights()` is declared-accepted and hand-written `untested`,
and that row is right as it stands.

## Defect 3: `structure_unit()` has no caller

`grep -rn structure_unit` over the whole monorepo, excluding `dev/`:

```
R/predict.R:2470          structure_unit_deviance(...)   <- different function
R/structure.R:849         structure_unit_deviance <- ... <- different function
R/structure.R:921         structure_unit <- function(st)
tests/testthat/test-structure.R:47,49                    <- the only callers
```

No package code calls it. The `unit` slot it reads is live: three
families declare one (`mixture()` at `R/families.R:3188`, `hmm()` at
`extensions/frmtmb.latent/R/hmm.R:496`, `rw_delta` at
`extensions/frmtmb.learn/R/family.R:203`) and
`extensions/frmtmb.sample/R/loo.R:161` reads `st[["unit"]]` directly,
with its own `%||% "a group"` fallback rather than through the accessor.

So the dead thing is the accessor, not the slot, and the documented
contract for the slot is false in core: `?frmtmb_structure` says of
`unit` that "the core quotes it where it must explain that a
per-OBSERVATION quantity does not exist", and core quotes it nowhere.

---

# What landed

## Defect 1: `frmtmb_family(exclusive_aterms =)`

A set of addition-term VALUES that say the same thing to the density, at
most one of which may be supplied.

| where | what |
|---|---|
| `R/families.R:341` | the new argument, defaulting to `list()` |
| `R/families.R:350` | validated against `required_aterms` at construction |
| `R/families.R:394` | stored on the family |
| `R/families.R:5060` | `exclusive_aterm_groups()`, the normalizer |
| `R/families.R:5080` | `check_exclusive_aterms()` |
| `R/families.R:5127` | `mixture_exclusive_aterms()` |
| `R/families.R:5149` | `check_exclusive_aterms_supplied()`, the refusal |
| `R/frame.R:1263` | the call site |
| `R/families.R:3052` | `mixture()` composes the sets |
| `extensions/frmtmb.eam/R/wiener-family.R:388` | `wiener()` adopts it |

Spelled in VALUES (`dec`, `vint1`), not in terms, so that it lines up
with `required_aterms`, which the same family writes as
`list(c("dec", "vint1"))`. A bare character vector is ONE set here and a
conjunction there: the two arguments say different things about a list
of names, so a shared convention would have made one of them
unwritable.

Placed immediately after the required-groups check at frame assembly,
not with the allow-list at the end. It is the same declaration read the
other way round (did the datum arrive at all, versus did it arrive
twice), and running it early keeps a family's own `valid_y` from
reporting on values the density was never going to read: `wiener()`'s
`ddm_check_response()` reads through `ddm_indicator()`, which prefers
`dec`, so it would have validated the winning column and said nothing
about the ignored one.

The message names the spelling to keep, and group order is that
precedence:

```
wiener: `dec` and `vint1` are spellings of the same datum for this
family, and this response supplies both. The density reads `dec` and
would ignore `vint1`, so the two could disagree with nothing to say so.
Keep dec(<column>) and drop vint(<column>).
```

### Interaction with any-of `required_aterms` groups

The two compose without either knowing about the other:

| declared | neither supplied | one supplied | both supplied |
|---|---|---|---|
| required any-of only | refused | fits | fits, second unread |
| required any-of + exclusive | refused | fits | refused |
| exclusive only | fits | fits | refused |

Together they read as "exactly one", which is `wiener()`'s contract. One
combination IS refused, at construction rather than at frame assembly: a
CONJUNCTION (`required_aterms = c("dec", "vint1")`, both needed) whose
values are also declared exclusive describes a family no model could
satisfy, so `check_exclusive_aterms()` stops and writes out the
alternative spelling.

Opt-in rather than implied by an any-of group, and `gddm()` is why: it
reads `dec()` and `vint1` together, as the boundary and the condition
index (`gd_indicator()`), so an implied rule would refuse the model that
family exists for. The comment at
`extensions/frmtmb.eam/R/gddm.R:1017` records that.

### Which eam families it applies to, measured

| family | accepts | reads two spellings of one datum? |
|---|---|---|
| `wiener()` | dec, vint, weights | YES: `ddm_indicator()` prefers `dec`, falls back to `vint1` |
| `gddm()` | dec, vint, vreal, weights | no: `dec` is the boundary, `vint1` the condition |
| `lba()` | vint, weights | no: does not accept `dec` |
| `rdm()` | vint, weights, cens, trunc | no: does not accept `dec` |
| `wiener_gng()` | dec, vreal, weights, cens | no: does not accept `vint` |

So `wiener()` alone, and the other four are unchanged. A family that
does not ACCEPT the second spelling never reaches the exclusivity
question: the allow-list refuses the term first.

`mixture()` composes the sets by INTERSECTION where it composes the
allow-list by union. The asymmetry is the point: a term reaches the
density if any component reads it, but two spellings are
interchangeable only if EVERY component treats them so.

The first implementation did not do that. It compared whole sets with
`setequal()`, which is stricter than the rule and fails OPEN: a
component declaring `c("dec", "vint1", "vint2")` beside one declaring
`c("dec", "vint1")` does treat `dec` and `vint1` as one datum, but the
sets are unequal, so the rule was dropped and the mixture then accepted
both spellings together, which is the defect this argument exists to
close. Not reachable from any shipped family, since `wiener()` is the
only one that declares a set, and covered by no test at all. Both fixed
2026-09-08 after the review found it (`R/families.R:5138`). Measured now:

| components | composed set |
|---|---|
| same set | `dec/vint1` |
| same set, order flipped | `dec/vint1` |
| `{dec,vint1}` vs `{dec,vint1,vint2}` | `dec/vint1` (was: none) |
| `{dec,vint1,vint2}` vs `{dec,vint1}` + `{vint2,vreal1}` | `dec/vint1` |
| one component declares none | none |
| disjoint sets | none |

The last two are the cases that must keep nothing, and they do.

## Defect 2: refusals derived from the declarations

| where | what |
|---|---|
| `R/compat.R:583` | `compat_aterm_rules(accepts, existing)`, exported |
| `R/compat.R:654` | `compat_already_refused()`, the deference |
| `R/compat.R:686` | `compat_derived_note()`, also how a derived row is recognized |
| `R/compat.R:714` | `compat_core_family_accepts()` |
| `R/compat.R:939` | `compat_hand_rules_tbl()`, split out and cached |
| `R/compat.R:1727` | `frmtmb_compat_rules_tbl()` = hand + derived + contributed |
| `R/compat.R:1893` | `compat_declared_displacements()`, the audit |
| `extensions/frmtmb.eam/R/ddm-shared.R:90` | `ddm_accepts`, one source of truth |
| `extensions/frmtmb.eam/R/zzz.R:171,246` | eam derives its own rows |

**Cells move from `untested` to `refused`, and nothing else moves, at
either package size.** Measured by resolving the whole pair table with
and without the derived block. The size matters and a core NEWS entry is
read with core alone, so both are given:

| configuration | pairs | untested -> refused | anything else |
|---|---|---|---|
| `frmtmb` alone | 5256 | **104** | none |
| `frmtmb` + `frmtmb.eam` | 5735 | **147** | none |

By term, with both loaded: `dec()` 36, `trials()` 37, `vint()` 36,
`vreal()` 38, summing to 147. Core alone has no `dec()` in the
vocabulary at all: `trials()` 32, `vint()` 36, `vreal()` 36, summing to
104. Every one is a family x addition-term pair; no other kind of pair
is touched, and 7 of the 147 belong to the eam families themselves.

The declarations refuse **295** family-by-addition-term cells in total.
147 of those are won by a derived row; the other **148** already had a
hand-written rule and keep its note.

CORRECTION, 2026-09-08, from the punch review. This section first
reported a flat "148 cells move", which is the DEFERRED half of that
295 split rather than the moved half, and with it `dec()` 37 and 8 eam
cells. All three came from one miscounted row in the measurement probe,
never from the shipped code; see "the defect was in the instrument, not
the product" at the top of this file.

### Three rules the derivation obeys

1. **Only the refusal is derived.** A family that ACCEPTS a term has
   said nothing about whether the pair works. `gddm` x `weights()` is
   declared-accepted and hand-written `untested`, and it stays
   `untested`.
2. **A pair already refused keeps its own note.** This is the second
   version of the code. The first derived a row for every non-accepted
   term, which is what the declaration supports, and it made **21**
   hand-written rules unreachable, 10 in core and 11 in `frmtmb.eam`.
   The biggest are `se()` x `kind:family` and `mi()` x `kind:family` at
   39 cells each, `trunc()` x `kind:family` at 24, `cens()` x
   `kind:family` at 13 and `cens()` x `group:discrete` at 11; the other
   16 win 1 to 4 cells each. Their notes say more than "the family did
   not list it" (the missing AD log-CDF, the residual variance only two
   families have), and a derived name x name row outranks the kind rule
   carrying them. So `compat_already_refused()` resolves each candidate
   pair against the hand-written rules first and skips the ones already
   refused. The derived block halved, 295 candidates to 147 rows, and
   every note survived: all 21 rules still win exactly the cells they
   won before. (This paragraph first said 8 rules, understating in the
   lane's own favour; corrected 2026-09-08 after the review measured
   21.)
3. **A displacement is a defect, not a resolution.**
   `compat_declared_displacements()` resolves the table with and without
   the derived block and reports any pair where a derived refusal takes
   a cell away from a rule that did not say `refused`. `test-compat.R`
   asserts it is empty. It recognizes contributed derived rows too, by
   recomputing `compat_derived_note()` from the pair, so an extension
   deriving its own rows is audited by the same check.

### One rule retired

`r("trials()", "kind:family", "untested", "trials() is meaningful only
for the binomial-type families.")` won no pair after the derivation: the
four families that take `trials()` are named by exact rules and the
other 32 refuse it by declaration. An unreachable rule is a claim nobody
can read, which `test-compat.R` already refuses, so it is gone and its
meaning was rewritten onto `trials() x group:trials_families`, which
does win. Not moved verbatim: the new note says which families take the
term and that every other refuses it by name. The only note change
anywhere else in the table follows from it, `beta_binomial` x `trials()`
and `zero_inflated_binomial` x `trials()` going from an empty note to
that sentence, because those two are the only families the group rule
covers without a name-level rule of their own.

### Cost

`frm_compat()` went from 0.19 s to 0.11 s, not up. The derivation adds
147 rules and one resolve over the 295 candidate pairs; caching
`compat_hand_rules_tbl()`, which the split forced, more than pays for
both.

## Defect 3: `structure_unit()` got its caller

`R/structure.R:924` (`structure_generic()`), `R/structure.R:945`
(`structure_unit()`).

Wired up rather than deleted, because the dead code was the accessor and
not the slot. `unit` is declared by three families and read by
`frmtmb.sample`, and `?frmtmb_structure` promises that "the core quotes
it where it must explain that a per-OBSERVATION quantity does not
exist". Deleting `structure_unit()` would have left that promise false
with nothing in core able to keep it.

The sentence has now been wrong twice, in OPPOSITE directions, and the
final wording asserts no mathematics at all:

```
main:  ... its likelihood does not factorize over the rows of the data,
       so the quantity this needs per row is not defined
first: ... its likelihood factorizes no finer than a hidden-Markov
       fix    sequence, so the quantity this needs per row is not defined
now:   ... the family replaces the rowwise likelihood with one the core
       evaluates whole, and it has not declared this capability. Its own
       unit for leaving data out is a hidden-Markov sequence
```

Main's version was false for a structure carrying no `loglik`, whose
likelihood IS rowwise. My first fix caught that and introduced the
opposite error: it read `unit` as a factorization claim, and `unit`
answers a different question. Its own `@param` at `R/structure.R:274-282`
and `dev/structured-family-protocol.md:200-206` both say so and name the
counterexample, `rw_delta`, whose finest factorization is the TRIAL and
whose leave-one-out unit is a whole subject.

The review made that sharper than a documentation contradiction: a
sibling lane has since declared `loglik_row` AND `loglik_group` on all
eight `frmtmb.learn` families, with the equality to `loglik` established
on the tape. Core would therefore have shipped a refusal telling a user
that a family "factorizes no finer than one subject's trial sequence"
while that same family declares, in the slots the protocol defines for
exactly that question, that it factorizes per row and per group. That is
a core message contradicting a shipped extension's own declaration.

The wording now reports two DECLARATIONS: the family hands the core one
whole-response likelihood, and it has not opted this capability in. The
unit is still quoted, because that is the slot's caller, but labelled
with the question it answers. The `loglik = NULL` branch is unchanged
and still says "the family declares it unsupported and gives no reason
of its own". The fallback for a family declaring no unit changed to "one
it has not named", which is what reads in the new frame.

NO TEST anywhere asserted this sentence, in core or in any extension,
which is how a falsehood survived a rewrite whose stated purpose was
removing one. `tests/testthat/test-structure.R` now has three, and what
they pin is the PROPERTY rather than the prose: for a family carrying
`loglik`, `loglik_row`, `loglik_group` and a per-subject `unit`, the
message must contain none of "factorize", "factorizes",
"factorization", "per row is not defined".

## The declaration-versus-row audit, in full

Two audits, one per package, both run BEFORE any change and both
re-run after.

### Core: 36 families x 8 addition terms

Every family in the compatibility vocabulary was constructed from
`family_registry`, its `accepted_aterm_names()` read off, and each pair
looked up in `frm_compat()`.

| table said | declaration says | cells | ruling |
|---|---|---|---|
| refused | refuses | 123 | agree, nothing to do |
| untested | refuses | 102 | **table wrong**, derive the refusal |
| works | accepts | 53 | agree |
| conditional | accepts | 2 | agree |
| (any) | unanswerable | 8 | `multinomial()` needs `K`; see below |

**No disagreement.** Not one cell where the table promised `works` or
`conditional` and the declaration refuses, and not one where the table
refused and the declaration accepts.

The eight unanswerable cells were `multinomial`, whose constructor
takes `K` and errors when called bare. Dropping the family would have
left its `vint()` and `vreal()` cells `untested` for no better reason
than an argument, so `compat_family_ctor_args` supplies `K = 2` with a
comment saying why any legal value serves: the object is read for its
`accepts_aterms` and nothing else, and no family's allow-list depends on
a constructor argument.

The one contradiction the previous round reported, `poisson` declaring
`cens()` while `cens()` is refused for discrete responses, is already
gone at this head. `fam_poisson()` declares
`accepts_aterms = c("weights", "trunc")` and the hand-written
`r("cens()", "poisson", "refused", ...)` agrees with it. Nothing
replaces it: there was no second contradiction to find.

### frmtmb.eam: 5 families x 9 addition terms

| family | hand-written rows | declaration | ruling |
|---|---|---|---|
| `wiener` | dec works, vint works, cens refused, trunc refused, weights works | dec, vint, weights | agree |
| `gddm` | dec works, vint works, vreal works, cens refused, trunc refused, weights UNTESTED | dec, vint, vreal, weights | agree |
| `lba` | dec refused, vint works, cens refused, trunc refused, weights UNTESTED | vint, weights | agree |
| `rdm` | vint works, dec refused, vreal refused, cens works, trunc conditional, weights works | vint, weights, cens, trunc | agree |
| `wiener_gng` | dec works, vint refused, vreal works, cens conditional, trunc refused, weights works | dec, vreal, weights, cens | agree |

**No disagreement here either.** The two rows in capitals are the case
worth naming: `gddm` x `weights()` and `lba` x `weights()` are
declared-ACCEPTED and hand-written `untested`, and they stay `untested`.
That is not a contradiction and the derivation must not touch it. The
declaration says the density can read the term; it says nothing about
whether the pair was ever run, and `untested` means exactly "nobody
ran it".

Two eam rows already cited the declaration in their notes before this
round (`lba` x `dec()` and `rdm` x `dec()`, both written after the
0.53.0 allow-list landed). Both are hand-written, both say more than the
derivation would, and both survive: `compat_already_refused()` sees the
pair is refused and skips it. `test-family.R` pins `rdm` x `dec()`'s
note for that reason.

## What was refused, and why

* **Deriving `works` from an allow-list.** A family that accepts a term
  has declared that the density can read it, not that the pair has ever
  been fitted. Turning acceptance into `works` would convert 53 measured
  cells and an unknown number of unmeasured ones into promises nothing
  checks, which is the failure the third status exists to prevent.

* **Implying exclusivity from an any-of `required_aterms` group.** It
  would have been a one-line rule with no new argument, and it is wrong
  for `gddm()`, which reads `dec()` and `vint1` together as two data.
  The protocol document already named that family as the obstacle; the
  measurement confirms it, since `gd_indicator()` returns both.

* **Deriving a refusal for every non-accepted term.** This was the first
  implementation and it is what the declaration literally supports. It
  cost 8 hand-written rules their reachability, among them
  `cens() x kind:family` and `se() x kind:family`, whose notes name the
  missing log-CDF and the residual variance. The declaration's note is
  true but poorer, and a name x name pair outranks the kind rule
  carrying the better sentence. Replaced with deference: derive only
  where the registry does not already refuse.

* **Deleting `structure_unit()`.** The accessor is dead; the slot it
  reads is not. Deleting it would have left `?frmtmb_structure`'s
  promise about `unit` unkeepable in core.

* **Touching `R/links.R`, the ordinal link switch in `R/families.R`, the
  `se()` and `cens()` gates in `R/frame.R`, and the band code in
  `R/predict.R`.** Sibling lanes own them. The `exclusive_aterms` check
  is a new call in `R/frame.R`'s response loop, 25 lines above the
  `mi()` block and well clear of the `se()`/`cens()` gates.

* **Composing `required_aterms` for `mixture()`.** Noticed while adding
  the exclusivity composition at `R/families.R:3052`: `mixture()` unions
  its components' `accepts_aterms` but passes no `required_aterms` at
  all, so `mixture(wiener(), lognormal())` does not inherit the
  boundary-indicator requirement. It is a real gap of the same family,
  but it is not one of this lane's three defects and closing it changes
  which models `mixture()` accepts. Left for whoever takes it, recorded
  here so it is not lost.

  The review confirmed the gap is wider than I wrote: such a mixture
  drops the EXCLUSIVITY as well as the requirement. It also confirmed
  this is not a regression, since that model behaves bit-identically
  before and after this lane, and that the direction of travel is right:
  `mixture(wiener(), wiener())` fits silently on main and is refused
  here.

## A note on the verification harness

The first two attempts at the one-process-per-file suite produced logs
that could not be audited, and the cause is worth recording because the
next lane will hit it.

Stopping a background task stops the task the harness tracks; it does
NOT kill that task's descendants. The `at-suite.sh` loops and their
`Rscript` children kept running after the stop, and a later run
truncating the same log then interleaved with them. The tell was a
torn line in the middle of the file:

```
ULT test-mean-fn.R NO-RESULT-LINE (process died)
```

`RES` and `ULT` split across two writers. 15 stray processes were still
alive, four of them `at-suite.sh` loops against two log paths.

Two fixes, both in the harness rather than in the procedure:

* The runner REFUSES to start when its log already exists. Two runs
  cannot share a file, so an interleave cannot happen even if a stray
  survives. Each run writes `at-<pkg>-run1.log`.
* Stray processes are found by their command line, which carries the
  `at-` prefix of this lane's scripts, and killed by PID. Nothing
  without that prefix is touched, because sibling lanes share the
  machine.

The counts reported below come from a single run of each suite against
one installed build, started after that cleanup. Nothing from the
earlier attempts is quoted.


## Verification (punch round, 2026-09-08)

Everything below is ONE run of each suite and ONE check per package,
against a single installed build, after every punch-list fix. The
earlier round's numbers are superseded.

### Suites, one process per test file, audited by name

| package | result lines | files on disk | duplicates | missing | fail | err | pass | skip |
|---|---|---|---|---|---|---|---|---|
| `frmtmb` | 131 | 131 | 0 | 0 | 0 | 0 | 5415 | 225 |
| `frmtmb.eam` | 17 | 17 | 0 | 0 | 0 | 0 | 1198 | 20 |

One testthat warning, in `test-prior-compat.R`, test "coef and group
narrow the classes that read them":
`Optimizer did not report convergence: singular convergence (7)`, from
that test's own `y ~ x + z + (x | g) + (z | h)` on 300 rows. The test
passes; it inspects which prior entries a class resolves to, not
convergence. Unreachable from this lane: that model carries no addition
term, no mixture and no structure.

### R CMD check --as-cran, WITH THE MANUAL

No `--no-manual`, no `--library=`, and both
`C:/Users/adf44/AppData/Roaming/TinyTeX/bin/windows` and
`C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools` on PATH.
The script prints where it found `pdflatex` and `pandoc` before starting,
because a run that cannot find them drops the manual stages silently.

| package | stages | manual PDF | tests | vignettes | Status |
|---|---|---|---|---|---|
| `frmtmb` | **58** | 546,166 bytes | OK 277 s | OK 347 s | **1 NOTE** |
| `frmtmb.eam` | **54** | 250,097 bytes | OK 433 s | OK | **1 NOTE** |

The only not-OK stage in either is:

```
* checking HTML version of manual ... NOTE
Skipping checking math rendering: package 'V8' unavailable
```

That comes from `R CMD check`'s own probe for `V8`, which is not
installed on this machine, and from no package. The stage count is
reported because it distinguishes a run that BUILT the manual from one
that skipped it: both runs here show `checking PDF version of manual`
and `checking HTML version of manual`, and both wrote the PDF. A run at
roughly 48 to 52 stages reporting zero notes has skipped the manual and
is not comparable.

NO example-timing note in either run, and the machine was busy: sibling
lanes held several R processes throughout, and the two chains ran
concurrently. An earlier quieter run of the same core tarball printed
`varCorr 4.47 0.50 5.07` under `checking examples ... OK`, which is
informational output below the note threshold rather than a NOTE. A
timing note is a property of machine load, not of a diff, and should be
checked against main before being attributed to one.

`checking for unstated dependencies in 'tests' ... OK` in BOTH runs with
no `R_BIOC_VERSION` set. The pin is deleted.

### Roxygen

Idempotent on both packages: two consecutive `roxygenise()` runs give
byte-identical output. Four generated files differ from HEAD:
`NAMESPACE`, `extensions/frmtmb.eam/NAMESPACE`, `man/frmtmb_family.Rd`,
`man/frmtmb_register_compat.Rd`.
