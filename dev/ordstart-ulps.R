# Lane ordstart: the assertions of tests/testthat/test-mixture-starts.R
# over 24 ulp perturbations of x (x * (1 + 2e-15 k), k = 0..23), so that
# each holds on every rounding path, not one seed's.
#   cum_sratio seed 9 (step default, code 0): the fit has no step function, or
#     the fit-end warning says no start escaped one.
#   cum_cum seed 4 (proper default): with starts, no degenerate warning
#     and logLik at least the starts = 0 fit's less 1e-3.
#   Rscript dev/ordstart-ulps.R lane|base
args <- commandArgs(trailingOnly = TRUE)
arm <- if (length(args)) args[1] else "lane"
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("tests/testthat/test-mixture-starts.R")
eval(parse(text = src[grep("^oms_data <- ", src):
                        (grep("^test_that", src)[1] - 1L)]))
ok1 <- ok2 <- 0L
for (k in 0:23) {
  r <- oms_fit(oms_data(9, k = k), mixture(cumulative(), sratio()))
  dg <- frmtmb:::mixture_ord_degeneracy(r$fit, "y")
  step <- dg$component[(dg$sharpen > -0.1) %in% TRUE]
  pass <- if (length(step)) {
    any(grepl("None of 10 jittered starts", r$msgs, fixed = TRUE))
  } else TRUE
  ok1 <- ok1 + pass
  cat(sprintf("STEP k %2d ll %.4f step %s starts-warn %s pass %s\n", k,
              as.numeric(logLik(r$fit)), paste(step, collapse = ","),
              any(vapply(r$warnings, inherits, NA,
                         "frmtmb_mixture_step_start")), pass))
  d <- oms_data(4, k = k)
  fam <- mixture(cumulative(), cumulative())
  ctl0 <- if (arm == "lane") frmtmb_control(starts = 0) else frmtmb_control()
  off <- oms_fit(d, fam, control = ctl0)
  on <- oms_fit(d, fam)
  dg <- function(x) any(grepl("ended at a degenerate boundary", x$msgs,
                              fixed = TRUE))
  pass <- dg(off) || (!dg(on) && as.numeric(logLik(on$fit)) >=
                        as.numeric(logLik(off$fit)) - 1e-3)
  ok2 <- ok2 + pass
  cat(sprintf("PROPER k %2d off %.4f on %.4f off-degen %s on-degen %s pass %s\n",
              k, as.numeric(logLik(off$fit)), as.numeric(logLik(on$fit)),
              dg(off), dg(on), pass))
}
cat("SUMMARY step-default pass", ok1, "of 24; proper-default pass", ok2,
    "of 24\n")
