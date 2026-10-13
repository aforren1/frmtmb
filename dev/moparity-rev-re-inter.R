# Reviewer: group-level mo() INTERACTION columns against brms 2.23.0's
# compiled program (joint density at frmtmb's estimates and modes):
# y ~ mo(x) * z + (mo(x) * z | g) and y ~ mo(x) * f2 + (mo(x) * f2 | g),
# where brms's Stan code is right (one covariate column).
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages({library(frmtmb); library(testthat)})
cat("lib:", find.package("frmtmb"), "\n")
td <- "C:/Users/adf44/source/r/frmtmb-wt-moparity/tests/testthat"
env <- new.env(parent = asNamespace("frmtmb"))
for (h in list.files(td, "^helper.*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = env)
}
options(frmtmb.brms_lp_report = TRUE)
set.seed(4)
n <- 400; n_g <- 25
dd <- data.frame(x = sample(0:3, n, TRUE), z = rnorm(n),
                 f2 = factor(sample(c("a", "b"), n, TRUE)),
                 g = factor(sample(paste0("g", 1:n_g), n, TRUE)))
eff <- c(0, 1, 1.5, 1.7)
bg <- rnorm(n_g, 0, 0.5); cg <- rnorm(n_g, 0, 0.3)
dd$y <- 1 + (1 + bg[dd$g]) * eff[dd$x + 1] +
  (0.5 + cg[dd$g]) * eff[dd$x + 1] * dd$z +
  0.4 * (dd$f2 == "b") * eff[dd$x + 1] + rnorm(n, 0, 0.5)
run <- function(fo) {
  cat("\n==", deparse1(fo), "\n")
  fit <- suppressWarnings(frm(bf(fo), family = gaussian(), data = dd))
  cat("ranef cols:", paste(colnames(ranef(fit)$g), collapse = " "), "\n")
  r <- tryCatch(
    test_that(deparse1(fo), env$brms_lp_check(
      brms::bf(fo), gaussian(), dd, fit, joint = TRUE,
      const = length(fit$frame$linpreds[[1]]$mo) * lgamma(3))),
    error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
  cat("test_that result:", isTRUE(r), "\n")
}
# first case passed in the earlier run (re-inter.txt, first block)
run(y ~ mo(x) * f2 + (mo(x) * f2 | g))
