# Reviewer, claim 4 follow-up: rate() under the identity link, with a
# plain column instead of I(x^2) (the harness names I() columns its own
# way), and a no-rate control through the same harness. Seed 306 data.
# Log: dev/aterms2-rev-log-03b-identity.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-aterms2"
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/aterms2-rev-stan-cache",
           FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
source(file.path(wt, "tests/testthat/helper-brms.R"))
run <- function(label, expr) {
  cat("\n==========", label, "\n")
  r <- tryCatch(expr, error = function(e) {cat("ERROR:", conditionMessage(e), "\n"); NULL},
                expectation_failure = function(e) {cat("EXPECTATION FAILED:", conditionMessage(e), "\n"); NULL})
  if (!is.null(r)) cat(sprintf("measured_const %.6g max_grad %.3g ours %.10g\n", r$measured_const, r$max_grad, r$ours))
  invisible(r)
}
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(306)
n <- 250
dr <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4), wt = runif(n, 0.5, 2))
dr$y <- rpois(n, exp(0.4 + 0.3 * dr$x) * dr$time)
dr$yn <- rnbinom(n, mu = (2 + 0.5 * dr$x^2) * dr$time, size = 3 * dr$time)
dr$x2 <- dr$x^2
f0 <- q(frm(yn ~ 1 + x2, data = dr, family = negbinomial("identity")))
run("control negbin identity, no rate", brms_lp_check(brms::bf(yn ~ 1 + x2), brms::negbinomial("identity"), dr, f0))
f1 <- q(frm(yn | rate(time) ~ 1 + x2, data = dr, family = negbinomial("identity")))
run("rate negbin identity", brms_lp_check(brms::bf(yn | rate(time) ~ 1 + x2), brms::negbinomial("identity"), dr, f1))
f2 <- q(frm(yn | rate(time) ~ 1 + x2, data = dr, family = geometric("identity")))
run("rate geometric identity", brms_lp_check(brms::bf(yn | rate(time) ~ 1 + x2), brms::geometric("identity"), dr, f2))
f3 <- q(frm(yn | rate(time) ~ 1 + x2, data = dr, family = poisson("identity")))
run("rate poisson identity", brms_lp_check(brms::bf(yn | rate(time) ~ 1 + x2), poisson("identity"), dr, f3))
cat("\nI(x^2) column names: frmtmb", colnames(f1$frame$linpreds[[1]]$X), "\n")
