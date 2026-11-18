# Reviewer of lane ordstarts, build round. Probe 6 (question 2): a pair the
# build calls distinct that is the same model. With a symmetric link,
# sratio and cratio are one model: P(Y = j | Y >= j) = F(tau_j - eta) for
# sratio and 1 - F(eta - tau_j) = F(tau_j - eta) for cratio. So swapping
# them is label switching, and the swapped start should sit at the fit's
# objective. cloglog (asymmetric) is the control.
#   Rscript ordstarts-build-review-p6.R > ordstarts-short-p6.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
mk <- function(seed, n) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
q <- function(e) suppressWarnings(suppressMessages(e))
for (lk in c("logit", "probit", "cloglog")) {
  for (seed in 1:5) {
    d <- mk(seed, 300)
    fam <- mixture(sratio(link = lk), cratio(link = lk))
    t0 <- proc.time()[["elapsed"]]
    a0 <- q(frm(bf(y ~ x + z), family = fam, data = d,
                control = frmtmb_control(starts = 0)))
    t1 <- proc.time()[["elapsed"]]
    dd <- q(frm(bf(y ~ x + z), family = fam, data = d))
    t2 <- proc.time()[["elapsed"]]
    sw <- ns$mixture_ord_swap_starts(a0, "y"); nd <- length(ns$mixture_swap_distinct(a0, sw))
    dv <- if (length(sw)) a0$obj$fn(sw[[1]]$par) - a0$obj$fn(a0$opt$par) else
      NA
    ms <- dd$opt$mixture_starts
    cat(sprintf(paste0("%-7s seed %d | swaps built %d, distinct %d | objective at the",
                       " swapped start - at the fit %+.3g | D - A0 ll %+.2e",
                       " | swapped kept %s | cost D/A0 %.2f\n"),
                lk, seed, length(sw), nd, dv,
                as.numeric(logLik(dd)) - as.numeric(logLik(a0)),
                if (is.null(ms)) "-" else ms[["swapped"]],
                (t2 - t1) / (t1 - t0)))
  }
}
cat("DONE p6\n")
