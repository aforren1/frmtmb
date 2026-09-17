# Reviewer (brmsport): do some recorded passes fail on a broken frmtmb?
#   Rscript dev/brmsport-rev-mutants.R > dev/brmsport-log/rev-mutants.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/brmsport-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(testthat); library(frmtmb) })
h <- new.env()
sys.source("tests/testthat/helper-brms-suite.R", envir = h)
holds <- function(expr) {
  ok <- tryCatch({ force(expr); TRUE },
                 expectation_failure = function(e) FALSE,
                 error = function(e) FALSE)
  ok
}
say <- function(label, v) cat(sprintf("%-70s %s\n", label, v))

cat("== families:4/5 student(identity)$link: a constructor that ignores",
    "its link\n")
student_mut <- function(link = "identity", ...) frmtmb::student()
say("real   student(identity)$link == 'identity'",
    holds(expect_equal(frmtmb::student(identity)$link, "identity")))
say("MUTANT student(identity)$link == 'identity' (want FALSE if genuine)",
    holds(expect_equal(student_mut(identity)$link, "identity")))
say("real frmtmb honours a non-default: student(log)$link",
    frmtmb::student(log)$link)
say("real frmtmb honours a non-default: bernoulli(probit)$link",
    frmtmb::bernoulli(probit)$link)
say("real frmtmb honours a non-default: exponential(identity)$link",
    frmtmb::exponential(identity)$link)

cat("== standata:163 rcens: the shim writes the zeros\n")
set.seed(91)
dat <- data.frame(y = rnorm(9), c4 = c(sample(-1:1, 5, TRUE), rep(2, 4)))
fr <- frm(y | cens(c4, y + 2) ~ 1, data = dat, dry_run = "frame")
v <- h$brms_standata_view(fr)
say("real   rcens == c(rep(0, 5), y[6:9] + 2)",
    holds(expect_equal(v$rcens, as.array(c(rep(0, 5), dat$y[6:9] + 2)))))
fr2 <- fr
av <- fr2$aterm_values$y
av$cens_y2[av$cens != 2] <- 999
fr2$aterm_values$y <- av
v2 <- h$brms_standata_view(fr2)
say("MUTANT y2 = 999 on non-interval rows: still holds?",
    holds(expect_equal(v2$rcens, as.array(c(rep(0, 5), dat$y[6:9] + 2)))))
av$cens_y2[av$cens == 2] <- av$cens_y2[av$cens == 2] + 1
fr2$aterm_values$y <- av
v3 <- h$brms_standata_view(fr2)
say("MUTANT y2 shifted on interval rows (want FALSE)",
    holds(expect_equal(v3$rcens, as.array(c(rep(0, 5), dat$y[6:9] + 2)))))

cat("== standata:77 bernoulli TRUE/FALSE codes\n")
fb <- frm(y ~ 1, data = data.frame(y = rep(c(TRUE, FALSE), 5)),
          family = "bernoulli", dry_run = "frame")
say("real   Y == rep(1:0, 5)",
    holds(expect_equal(h$brms_standata_view(fb)$Y, as.array(rep(1:0, 5)))))
fb$y$y <- 1 - fb$y$y
say("MUTANT codes flipped (want FALSE)",
    holds(expect_equal(h$brms_standata_view(fb)$Y, as.array(rep(1:0, 5)))))

cat("== priors:139 as.brmsprior drops an unknown column\n")
bp <- as.brmsprior(data.frame(prior = "normal(0,1)", x = "test",
                              coef = c("a", "b")))
say("real   bprior$x is NULL", holds(expect_equal(bp$x, NULL)))
bpm <- as.data.frame(bp); bpm$x <- "test"
say("MUTANT column kept (want FALSE)", holds(expect_equal(bpm$x, NULL)))
