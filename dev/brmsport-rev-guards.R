# Reviewer (brmsport): hollow passes the harness's rules do not catch.
# Every case is recorded as 'pass' and is hollow, so a guard that
# worked would report 1 failure. "failures 0" means the hollow pass is
# accepted. The last two are controls that must report 1.
#   Rscript dev/brmsport-rev-guards.R > dev/brmsport-log/rev-guards.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/brmsport-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
Sys.unsetenv("FRMTMB_BRMSPORT_RECORD")
h <- new.env()
sys.source("tests/testthat/helper-brms-suite.R", envir = h)

run <- function(label, expr) {
  f <- tempfile(fileext = ".R")
  writeLines(c("test_that('g', {",
               "set.seed(91)",
               "d <- data.frame(y = rnorm(20), x = rnorm(20))",
               "fit <- frmtmb::frm(y ~ x, d)",
               expr, "})"), f)
  r <- as.data.frame(test_file(f, reporter = "silent",
                               env = new.env(parent = h)))
  got <- sum(r$failed) + sum(r$error)
  cat(sprintf("%-66s failures %d\n", label, got))
}

run("H1 regex expect_error met by an argument-name refusal",
    "brms_port('h:1', 'pass', '', expect_error(fitted(fit, bogus = 1), 'bogus'))")
run("H2 bare expect_error met by a 'cannot honor' argument refusal",
    "brms_port('h:2', 'pass', '', expect_error(fitted(fit, ndraws = 5)))")
run("H3 expect_equal on two NULLs read through fit$data",
    "brms_port('h:3', 'pass', '', expect_equal(fit$data$y, fit$data$x))")
run("H4 regex met by R's missing-argument error (brmsfit-methods:698)",
    "brms_port('h:4', 'pass', '', expect_error((function(group) group)(), 'group'))")
run("H5 stale object through a complex assignment x$a <- ",
    c("x <- list(a = 1)",
      "brms_setup('h:5', x$a <- stop('setup broke'))",
      "brms_port('h:6', 'pass', '', expect_equal(x$a, 1))"))
run("H6 missing function spelled 'of mode function was not found'",
    "brms_port('h:7', 'pass', '', expect_error(match.fun('no_such_fn_xyz')(), 'no_such'))")
run("H7 is.numeric() on an NA column (brmsfit-methods:394)",
    "brms_port('h:8', 'pass', '', expect_true(is.numeric(NA_real_)))")
run("control: bare expect_error met by 'has no argument' (want 1)",
    "brms_port('h:9', 'pass', '', expect_error(vcov(fit, cor = TRUE)))")
run("control: simple stale object (want 1)",
    c("x <- 1", "brms_setup('h:10', x <- stop('broke'))",
      "brms_port('h:11', 'pass', '', expect_equal(x, 1))"))
