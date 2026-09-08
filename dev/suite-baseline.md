# The per-file suite baseline, and what it is for

`dev/suite-baseline.tsv` records one row per test file as of frmtmb
0.55.0: package, file, passing assertions, skips. It is a floor, not a
target.

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
