# Punch round 2 constructions: the two defects the recheck asked the lane
# to file (R1 and R5). frmtmb first, brms 2.23.0 second on the same call;
# brms only builds Stan code, so nothing compiles.
#
#   Rscript dev/brmsport-punch2-defects.R > dev/brmsport-log/punch2-defects.txt 2>&1
#
# Data: brm:110's own design, `dat` of tests.brm.R, seed 20260917.
.libPaths(c("C:/Users/adf44/source/r/brmsport-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
h <- new.env()
sys.source("tests/testthat/helper-brms-suite.R", envir = h)
cat("frmtmb", format(packageVersion("frmtmb")), "brms",
    format(packageVersion("brms")), "emmeans",
    format(packageVersion("emmeans")), "\n")

show <- function(label, expr) {
  out <- tryCatch({
    v <- suppressWarnings(suppressMessages(expr))
    paste(utils::capture.output(print(v)), collapse = " | ")
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("  %-44s %s\n", label, substr(gsub(" +", " ", out), 1, 230)))
}

set.seed(20260917)
dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))
dat$yc <- rpois(10, 3)

cat("\n== P1 ma(x) with brms's default cov = FALSE\n")
cat("  brm:110 is the poisson case; the gaussian and student cases are",
    "not in bin 1\n")
for (fam in c("gaussian", "student", "poisson")) {
  f <- if (fam == "poisson") yc ~ ma(x) else y ~ ma(x)
  famobj <- get(fam, envir = asNamespace("brms"))()
  show(sprintf("frmtmb frm(%s ~ ma(x), %s)", all.vars(f)[1], fam),
       class(frm(f, dat, family = get(fam)())))
  show(sprintf("brms make_stancode, %s", fam),
       sprintf("Stan code, %d characters",
               nchar(brms::make_stancode(f, dat, family = famobj))))
}
cat("  control: frmtmb with cov = TRUE on gaussian fits\n")
show("frmtmb frm(y ~ ma(x, cov = TRUE), gaussian)",
     logLik(frm(y ~ ma(x, cov = TRUE), dat)))

cat("\n== P2 emmeans hides frmtmb's refusal reason\n")
fit2 <- h$brms_fixture(2)
fit6 <- h$brms_fixture(6)
for (case in list(list("fit2 (nonlinear)", fit2, "Age"),
                  list("fit6 (multivariate)", fit6, "Age"))) {
  show(paste("emmeans::emmeans,", case[[1]]),
       emmeans::emmeans(case[[2]], case[[3]]))
  show(paste("frmtmb's own reason,", case[[1]]),
       emmeans::recover_data(case[[2]]))
}
cat("  control: a linear univariate fit reaches a grid\n")
fit_lin <- frm(y ~ x, dat)
show("emmeans::emmeans, y ~ x", nrow(summary(emmeans::emmeans(fit_lin, "x"))))
