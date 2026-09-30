# Reviewer, claim 2: a distributional sigma and mo() under subset()
# against brms, split so that the harness's missing translation rule for
# a multivariate simo_ (dev/aterms2-rev-log-03-brms.txt) is not hit.
# Seed 304 data. Log: dev/aterms2-rev-log-03e-sigma.txt
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
set.seed(304)
n <- 160
ds <- data.frame(x = runif(n, -2, 2), z = rnorm(n),
                 om = factor(sample(1:4, n, TRUE), ordered = TRUE),
                 s1 = rep(c(TRUE, TRUE, FALSE), length.out = n))
ds$y1 <- sin(ds$x) + as.integer(ds$om) * 0.3 + rnorm(n, 0, exp(0.2 * ds$z))
ds$y2 <- ds$z + rnorm(n)
f1 <- q(frm(bf(y1 | subset(s1) ~ x, sigma ~ z) + bf(y2 ~ z), data = ds,
            family = gaussian()))
run("mv subset, sigma ~ z", brms_lp_check(brms::bf(y1 | subset(s1) ~ x, sigma ~ z) +
      brms::bf(y2 ~ z) + brms::set_rescor(FALSE), gaussian(), ds, f1))
f2 <- q(frm(y1 | subset(s1) ~ mo(om) + x, data = ds))
run("univariate subset, mo()", brms_lp_check(brms::bf(y1 | subset(s1) ~ mo(om) + x),
      gaussian(), ds, f2))
# control: the same mo() model without subset(), on the subset's rows;
# a Dirichlet(1, 1, 1) simplex prior has log density log(2) = 0.6931
f3 <- q(frm(y1 ~ mo(om) + x, data = ds[ds$s1, ]))
run("control: no subset, mo(), rows of s1",
    brms_lp_check(brms::bf(y1 ~ mo(om) + x), gaussian(), ds[ds$s1, ], f3))
run("univariate subset, mo(), const = log(2)",
    brms_lp_check(brms::bf(y1 | subset(s1) ~ mo(om) + x), gaussian(), ds, f2,
                  const = log(2)))
