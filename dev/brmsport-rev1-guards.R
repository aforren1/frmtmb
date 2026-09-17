# Reviewer (brmsport), recheck round 1: the new hollow rules and the
# own-words verdict, against the punched helper.
#   Rscript dev/brmsport-rev1-guards.R > dev/brmsport-log/rev1-guards.txt 2>&1
# "accepted" means a pass verdict reported 0 failures.
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
               "d <- data.frame(y = rnorm(20), x = 1:20, z = rnorm(20))",
               "fit <- frmtmb::frm(y ~ z, d)",
               expr, "})"), f)
  r <- as.data.frame(test_file(f, reporter = "silent",
                               env = new.env(parent = h)))
  got <- sum(r$failed) + sum(r$error)
  cat(sprintf("%-68s failures %d %s\n", label, got,
              if (got == 0) "ACCEPTED" else "rejected"))
}

cat("== own words\n")
run("O1 cov refusal on GAUSSIAN ma(x), a model brms fits",
    paste0("brms_port_own('o:1', 'ma[(][)]: only the residual-covariance ",
           "formulation is implemented, so the call needs cov = TRUE', '', ",
           "expect_error(frmtmb::frm(y ~ ma(x), d), 'Please set cov = TRUE'))"))
run("O2 own words reading a stale object",
    c("x <- fit",
      "brms_setup('o:2', x <- stop('setup broke'))",
      paste0("brms_port_own('o:3', \"Unknown dpar: 'inv'\", '', ",
             "expect_error(fitted(x, dpar = 'inv'), 'Invalid argument'))")))
run("O2 control: the same row as a plain pass verdict (want rejected)",
    c("x <- fit",
      "brms_setup('o:4', x <- stop('setup broke'))",
      paste0("brms_port('o:5', 'pass', '', ",
             "expect_error(fitted(x, dpar = 'inv'), 'Unknown dpar'))")))
run("O3 pattern '.' matches a designed refusal of another case",
    paste0("brms_port_own('o:6', '.', '', ",
           "expect_error(frmtmb::frm(y ~ z, d, family = poisson()), ",
           "'brms words'))"))
run("O4 own words met by an argument refusal (want rejected)",
    paste0("brms_port_own('o:7', 'bogus', '', ",
           "expect_error(fitted(fit, bogus = 1), 'brms words'))"))

cat("== under-reach: hollow passes still accepted\n")
run("U1 expect_null on a slot frmtmb lacks",
    "brms_port('u:1', 'pass', '', expect_null(fit$data$y))")
run("U2 expect_true(is.null()) on a slot frmtmb lacks",
    "brms_port('u:2', 'pass', '', expect_true(is.null(fit$data$y)))")
run("U3 expect_type() on NA",
    "brms_port('u:3', 'pass', '', expect_type(NA_real_, 'double'))")
run("U4 is.numeric(v) && length(v) == 1 on NA",
    paste0("brms_port('u:4', 'pass', '', ",
           "expect_true(is.numeric(NA_real_) && length(NA_real_) == 1))"))
run("U5 expect_length(<absent slot>, 0)",
    "brms_port('u:5', 'pass', '', expect_length(fit$data$y, 0))")
run("U6 argument refusal in other words ('Cannot interpret bf() argument')",
    paste0("brms_port('u:6', 'pass', '', ",
           "expect_error(frmtmb::bf(y ~ x, sigma1 = 'sigma2'), 'sigma'))"))
run("U7 stale through assign()",
    c("x <- 1", "brms_setup('u:7', assign('x', stop('broke')))",
      "brms_port('u:8', 'pass', '', expect_equal(x, 1))"))
run("U8 expect_equal(NULL, <variable holding NULL>)",
    c("e <- NULL",
      "brms_port('u:9', 'pass', '', expect_equal(fit$data$y, e))"))

cat("== over-reach: genuine passes now rejected\n")
run("V1 two genuinely NULL names (unnamed vectors)",
    "brms_port('v:1', 'pass', '', expect_equal(names(c(1, 2)), names(c(3, 4))))")
run("V2 an assertion ABOUT an argument refusal",
    paste0("brms_port('v:2', 'pass', '', ",
           "expect_error(fitted(fit, bogus = 1), 'has no argument'))"))
run("V3 is.numeric() on a genuinely NA estimate brms also reports as NA",
    "brms_port('v:3', 'pass', '', expect_true(is.numeric(NA_real_)))")
