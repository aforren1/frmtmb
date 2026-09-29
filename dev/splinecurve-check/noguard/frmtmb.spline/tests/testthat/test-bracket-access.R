# `[[` on the hazard containers, asserted on THIS package.
#
# `$` on a list partial-matches, so a `$` read of an absent slot on a
# container whose slot names collide returns a neighbor rather than
# NULL. frmtmb owns the rule, the container list and the scanner; see
# `frmtmb::frm_hazard_reads()` and dev/bracket-sweep.md.
#
# This file is deliberately short and carries no list of its own. It is
# here because the guard used to run only in frmtmb's suite, and
# frmtmb.coupling passed its own `R CMD check` and its own tests with
# 34 hazard reads before the release tally caught them.
#
# It scans the installed NAMESPACE, so unlike frmtmb's source-tree
# scanner it runs under `R CMD check` on a built tarball. Top-level
# package code that is not inside a function is the one thing it cannot
# see; frmtmb's suite still reads the sources of every extension.

test_that("no frmtmb.spline function reads a hazard container", {
  hits <- frmtmb::frm_hazard_reads("frmtmb.spline")
  expect_identical(
    hits, character(0),
    info = paste0("Use [[\"name\"]], or rename the local: these names ",
                  "are reserved. See ?frmtmb::frm_hazard_reads.\n",
                  paste(hits, collapse = "\n")))
})
