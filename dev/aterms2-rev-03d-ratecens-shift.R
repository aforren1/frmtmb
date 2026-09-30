# Reviewer, claim 4: rate() + cens() against brms once frmtmb's
# documented right-censoring convention (P(Y >= k), R/families.R, brms
# P(Y > k)) is aligned by giving frmtmb y + 1 on the right-censored
# rows. Seed 306 data. Log: dev/aterms2-rev-log-03d-ratecens-shift.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-aterms2"
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/aterms2-rev-stan-cache",
           FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
source(file.path(wt, "tests/testthat/helper-brms.R"))
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(306)
n <- 250
dr <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4), wt = runif(n, 0.5, 2))
dr$y <- rpois(n, exp(0.4 + 0.3 * dr$x) * dr$time)
dr$cc <- rep(c(0, 0, 1, -1), length.out = n)
df <- dr
df$y[df$cc == 1] <- df$y[df$cc == 1] + 1L
f <- q(frm(y | rate(time) + cens(cc) ~ x, data = df, family = poisson()))
r <- tryCatch(brms_lp_check(brms::bf(y | rate(time) + cens(cc) ~ x),
                            poisson(), dr, f),
              error = function(e) {cat("ERROR", conditionMessage(e), "\n"); NULL})
if (!is.null(r)) cat(sprintf("rate+cens shifted: measured_const %.6g max_grad %.3g\n",
                             r$measured_const, r$max_grad))
