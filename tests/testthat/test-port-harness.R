# The ported brms suite's harness (helper-brms-suite.R) refuses an
# assertion that reads a name whose last assignment failed, because the
# name still holds an earlier object and brms's block would have stopped.
# A setup line that reads such a name and succeeds used to mark what it
# assigns as fresh, so the staleness lapsed one line later (the
# remeasure review, R7; dev/remeasure-bin2-stale-setup.R constructs it).

port_harness_rows <- function(body) {
  rec <- withr::local_tempfile(fileext = ".tsv")
  withr::local_envvar(FRMTMB_BRMSPORT_RECORD = rec,
                      FRMTMB_BRMSPORT_PKG = "frmtmb")
  body()
  x <- utils::read.delim(rec, header = FALSE, quote = "",
                         colClasses = "character")
  stats::setNames(as.list(x[[5L]]), x[[3L]])
}

test_that("a setup line that reads a stale name passes the staleness on", {
  held <- port_harness_rows(function() {
    brms_setup("x:1", a <- 1)
    brms_setup("x:2", a <- stop("the line brms's block would stop at"))
    # read directly: refused, as before
    brms_port("x:3", "pass", "", expect_equal(a, 1))
    # read through a setup line that itself succeeds
    brms_setup("x:4", b <- a + 1)
    brms_port("x:5", "pass", "", expect_equal(b, 2))
    # and the not-run form, a setup line that reaches brms's namespace
    brms_setup("x:6", prep <- list(data = list(Y = 1:3)))
    brms_setup_not_run("x:7", "prep")
    brms_setup("x:8", n <- length(prep$data$Y))
    brms_port("x:9", "pass", "", expect_equal(n, 3))
  })
  expect_identical(held[["x:3"]], "FALSE")
  expect_identical(held[["x:5"]], "FALSE")
  expect_identical(held[["x:9"]], "FALSE")
  # the setup line itself ran, and is recorded as having run
  expect_identical(held[["x:4"]], "TRUE")
})

test_that("a fresh assignment after the stale line clears it", {
  # the absent case: staleness must not outlive a successful reassignment
  # that reads nothing stale
  held <- port_harness_rows(function() {
    brms_setup("y:1", a <- stop("fails"))
    brms_setup("y:2", b <- a + 1)
    brms_setup("y:3", b <- 5)
    brms_port("y:4", "pass", "", expect_equal(b, 5))
    brms_setup("y:5", c0 <- 2)
    brms_port("y:6", "pass", "", expect_equal(c0 * 2, 4))
  })
  expect_identical(held[["y:4"]], "TRUE")
  expect_identical(held[["y:6"]], "TRUE")
})
