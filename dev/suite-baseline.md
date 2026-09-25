# The per-file suite baseline, and what it is for

`dev/suite-baseline.tsv` records one row per test file as of the 0.63.0
release, at frmtmb 0.63.0 and every extension's version of that release:
package, file, passing assertions, skips. It is a floor, not a target.

Every row is from the release run (`dev/release/suite.log`, 292 files).
13 files are new. One file's passing count fell against 0.62.0, by
design: frmtmb.learn's `test-counterfactual.R` went from 70 to 69,
because lane `wt-simnewdata` removed an assertion that `simulate()` has
no `newdata` argument, which it now has. No other file fell.

Read a tier log only when it POSTDATES every file it covers. Three lanes
of this round quoted a suite log written before their last edit, and in
each case the arithmetic gave the lane away: a file's count in the log
differed from its standalone count by exactly the assertions the last
edit added.

## Why it exists

A release audit that counts FILES cannot see a file that ran and
asserted less than it should. When a `test_that()` block throws, testthat
records the error and abandons the rest of that block, so every
assertion after the throw is silently not run. A runner that reports
only `pass` and `fail` then prints a clean line for a file that checked
eight fewer things than it claims to.

That happened during the 0.55.0 round, in a lane's own runner: a
rewritten message broke an `expect_error()` substring, the block
aborted, and the file reported PASS=51 FAIL=0 where it asserts 59. The
release harness here reports `error` separately and the tallies filter
on it, so an errored file is caught. A file whose count merely DROPPED
is not, and that is what this baseline is for.

## How to use it

After a suite run, join the run's per-file counts against this file and
report any file whose passing count fell. A drop is not automatically a
defect: deleting a test, merging two files, or gating a block all lower
a count legitimately. It is a question that has to be answered, which is
the point. A count that rose needs no explanation.

Regenerate it at each release, from the release run, once the suite is
green:

    awk '/^== /{pk=$2; next} /^RESULT/{split($3,a,"="); split($6,d,"=");
         printf "%s\t%s\t%s\t%s\n", pk, $2, a[2], d[2]}' <suites log> \
      > dev/suite-baseline.tsv

## What it does not cover

The gated tiers are not in it. Their counts are recorded in the release
commit message instead, and they are compared by hand against the
previous release, which is how the 0.55.0 round confirmed that rewriting
the student-t density moved no brms or BCM figure. Folding them in wants
a stable way to name a tier, since the same file appears gated and
ungated with different counts.

## What the 0.55.2 round added

The baseline earned its place twice, in opposite directions.

It caught a real drop: `frmtmb.sample/test-sampling-ported.R` fell from
208 to 205. The answer was legitimate, a three-assertion block removed
because the same ground is now covered properly by
`test-tmbstan-build-guard.R`, which constructs a broken build instead of
matching a synthesized string. A drop is a question, and this one had an
answer.

And it was itself found stale. The eam rows summed to 1376 over 19
files rather than the 1408 recorded, and `test-units.R` was missing
altogether, so a regression in that file could not have been seen. A
baseline is only a floor while it is regenerated at every release, which
is what the instruction above says and what did not happen last time.

One trap worth naming, because the release harness fell into it twice
in this round. A runner that reports only `pass` and `fail` cannot tell
a green run from a run where nothing executed: a harness bug that made
every file fail to load reported `pass=0 fail=0 err=0` for all eight
packages, and a second bug that left the packages unattached reported
plausible-looking damage in the one package the round had rewritten.
Read the file count and the load errors before reading the failures.

## What the round 3 release added

Regenerated at 220 rows and 13124 assertions, from 216 and 12794. Four
files are new: `frmtmb.eam/test-ndt-seam.R` and
`frmtmb.eam/test-gddm-conditions.R`, `frmtmb.learn/test-rlddm-ndt.R`,
and `frmtmb.latent/test-hmm-starts.R`. No file was removed.

Exactly one count fell, and it had an answer:
`frmtmb.sample/test-tmbstan-build-guard.R` from 24 to 23. The file
carries a conditional `succeed()` on the branch where the installed
StanHeaders cannot be compiled against, and pinning StanHeaders to
2.32.10 outside the user library takes that branch away. The assertion
is gone because the condition it stood in for cannot occur on this
machine any more.

That is worth keeping in mind when reading a drop here: a count can
fall because the ENVIRONMENT changed rather than because the tests did,
and the diff of the test file will show nothing.

## What the Phase 2 round added

Regenerated at 221 rows and 13159 assertions, from 220 and 13124. One
file is new, `frmtmb.spline/test-frailty.R` at 28, which replaces a
smoke test that asserted only finiteness. No file was removed.

Exactly one count fell and it had an answer before anyone asked:
`frmtmb.spline/test-surface.R` from 48 to 47. Item 2.5 removed three
assertions and added two, so the arithmetic closes at 48 - 3 + 2. The
lane's reviewer had already reconciled it independently, which is the
cheapest form this question can take: a drop explained in the review is
a drop nobody has to re-derive at release time.

## What the Phase 2.5 release added

Regenerated at 223 rows and 13247 assertions, from 221 and 13159. No
count fell. Two files are new, `frmtmb/test-generic-collision.R` at 47
and `frmtmb/test-scale-contract.R` at 37, and one rose,
`frmtmb.sample/test-draws-methods.R` from 97 to 101 for the two tests
that pin the shadowing note on draws. That is 47 + 37 + 4 = 88, the
whole of the increase.

This is the release where a silent drop was most likely and did not
happen. The lane's own round found `test-naming-collisions.R` losing 8
of its 31 assertions under the first design, because `hypothesis()`'s
note stopped firing once its generic belonged to brms. A file that
still runs green while asserting less is exactly what this baseline
exists to catch, and it holds its count here.

## What the frmtmb.sample generics release added

Regenerated at 224 rows and 13314 assertions, from 223 and 13247. No
count fell. One file is new, `frmtmb.sample/test-generic-collision.R` at
58, and one rose, `frmtmb/test-generic-collision.R` from 47 to 56, which
gained a behavioral block replacing one that had failed only on an
unused argument. 58 + 9 = 67, the whole of the increase.

`frmtmb.sample/test-loo.R` holding at 79 is the row worth reading here.
That file broke during the lane because `local_mocked_bindings()` cannot
mock an active binding: assigning to one calls it. It was repaired by
mocking the method instead, and an unchanged count is the evidence that
the repair restored the assertions rather than removing them.

## What the brms-matching release added

Regenerated at 225 rows and 13485 assertions, from 224 and 13314. One
file is new, `frmtmb/test-arg-refusal.R` at 118.

TWO counts fell, and both were derived line by line in review before
release rather than explained after it. `frmtmb/test-api-spellings.R`
went from 36 to 33 because three assertions checked `re.form` aliases
that item 2.5e removed outright. `frmtmb/test-brms-methods.R` went from
18 to 16 as a split of 3 to 2 and 2 to 1, where warnings that had been
asserted became errors.

A third drop was a REGRESSION and is not in this table because it was
fixed first. `frmtmb.learn/test-counterfactual.R` fell from 67 to 65
with one error. `simulate()` never had a `newdata` argument, and at
0.57.0 the name was silently swallowed by `...`, so a test passing
`newdata = d` reached an unrelated refusal and looked as though newdata
had been considered. The new dots refusal named it, the test was
correcting the bug, and it now stands at 68. The lane that caused it
had declared plainly that it did not run the five other extensions'
suites, and that is exactly where the release suite found it.

## What the 2.6c release changed

Regenerated at 231 rows and 15278 assertions, from 225 and 13485. Six
files are new: `frmtmb/test-brms-families.R`,
`frmtmb/test-brms-formula-priors.R`, `frmtmb/test-brms-names.R`,
`frmtmb.sample/test-brms-output.R`, `frmtmb.sample/test-brms-pins.R`
and `frmtmb.sample/test-prior-update.R`.

Three counts fell, and each equals the count the lane that changed the
file recorded in its own full run before the merge:

- `frmtmb/test-naming-collisions.R`, 31 to 29. wt-brmsnames removed the
  reserved-name shadow machinery and rewrote the file for its absence.
- `frmtmb/test-prior-compat.R`, 195 to 194, and
  `frmtmb/test-priors-autocor-classes.R`, 63 to 62. wt-priorform merged
  several print assertions into one `expect_identical()` each.

## What the 2.6b, 2.6e and tmbstan release changed

Regenerated at 262 rows and 15567 assertions, from 231 rows and 15278.
31 files are new and nothing fell below its previous count.

The new files are the census and condition tests item 2.6e generates
into all eight packages (`test-conditions-census.R` everywhere, and
`test-conditions.R` in core and each extension), and item 2.6b's ported
bin-1 tier of brms's own suite, which SKIPS here and runs in the gated
tier (`dev/release/run-gated.ps1`).

Read the gated line beside this one. Gated is 38 files and 3127
assertions at this release, from 23 files and 2496, and it compiled 109
Stan programs from an EMPTY cache against StanHeaders 2.39.1. That is
the evidence that dropping the StanHeaders pin is safe; the ungated
suite compiles no Stan program and cannot show it.

