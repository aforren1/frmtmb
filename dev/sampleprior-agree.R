# Lane sampleprior: frm_sample()'s posterior of a dispersion parameter
# against brms 2.23.0's under each one's default priors, on the four
# cases of dev/sampleprior-agree-common.R.
#
#   Rscript dev/sampleprior-agree.R <lane|r9> [case ...]
#
# Seed policy as test-loo.R's brms agreement: 2 chains, 2000 iterations,
# seed 3 on both sides. The comparison is of log(parameter), whose
# posterior is closer to normal; z is the difference of the two means
# over the pooled Monte Carlo standard error of the means.
a <- commandArgs(TRUE)
arm <- a[1]
which <- a[-1]
lib <- if (identical(arm, "r9")) character(0) else
  "C:/Users/adf44/source/r/wt-sampleprior-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sampleprior"
cache <- file.path(wt, "dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source(file.path(wt, "dev/sampleprior-agree-common.R"))
cat("arm", arm, "lib:", find.package("frmtmb.sample"), "\n")
cases <- agree_cases()
if (length(which)) cases <- cases[which]
chains <- 2L
iter <- 2000L
for (nm in names(cases)) {
  cs <- cases[[nm]]
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(agree_frm(cs, 3L, chains, iter),
                error = function(e) conditionMessage(e))
  t1 <- proc.time()[["elapsed"]]
  b <- agree_brms(cs, 3L, chains, iter, cache)
  if (is.character(f)) {
    cat(sprintf("%-12s %-6s frmtmb ERROR %s\n", nm, cs$par, f))
    next
  }
  lf <- log(f)
  lb <- log(b)
  mf <- agree_mcse(lf, chains)
  mb <- agree_mcse(lb, chains)
  z <- (mean(lf) - mean(lb)) / sqrt(mf^2 + mb^2)
  q <- function(v) paste(sprintf("%.3f", stats::quantile(v, c(0.05, 0.5,
                                                              0.95))),
                         collapse = " ")
  w90 <- function(v) diff(stats::quantile(v, c(0.05, 0.95)))
  wr <- unname(w90(lf) / w90(lb))
  rf <- posterior::rhat(matrix(lf, ncol = chains))
  rb <- posterior::rhat(matrix(lb, ncol = chains))
  cat(sprintf(paste0("%-12s log(%s): frmtmb mean %.4f sd %.4f mcse %.4f ",
                     "rhat %.3f q05/50/95 %s | brms mean %.4f sd %.4f ",
                     "mcse %.4f rhat %.3f q05/50/95 %s | z %.2f sd ratio ",
                     "%.3f 90%% width ratio %.3f | frmtmb %.0fs\n"),
              nm, cs$par, mean(lf), stats::sd(lf), mf, rf, q(lf), mean(lb),
              stats::sd(lb), mb, rb, q(lb), z, stats::sd(lf) / stats::sd(lb),
              wr, t1 - t0))
}
