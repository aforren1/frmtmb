# Review of lane ci701: the boundary-correlation fit of test-review-v29.R
# on ulp-moved data (seed 8000 + k, a different draw from the lane's
# 1000 + k), from the default start and from theta = (0, 0, 500).
# On the lane build each fit is made twice in the same process: as
# installed, and with se_edge_cor() stubbed to FALSE (0.70.0's check
# with D1 kept), so that the D2 effect is read without P/E-core noise.
#
#   Rscript dev/ci701-rev-batt.R <lane|base> <k0> <k1>
args <- commandArgs(TRUE)
arm <- args[1]; k0 <- as.integer(args[2]); k1 <- as.integer(args[3])
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
lane <- "C:/Users/adf44/source/r/wt-ci701-lib"
.libPaths(if (arm == "base") c(base, ul) else c(lane, base, ul))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "blas:", extSoftVersion()[["BLAS"]], La_library(), "\n")
set.seed(61)
ng <- 12
u <- stats::rnorm(ng, 0, 0.8)
per <- 5
dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                 x = stats::rnorm(ng * per))
dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.4)
dd$y2 <- -1 + 0.2 * dd$x + 0.6 * u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.3)
d0 <- dd
mdata <- function(k) {
  dd <- d0
  if (k > 0) {
    set.seed(8000 + k)
    for (v in c("x", "y1", "y2")) {
      dd[[v]] <- dd[[v]] * (1 + sample(0:4, nrow(dd), TRUE) * 2^-52)
    }
  }
  dd
}
f <- mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
          bf(y2 ~ x + (1 | q | id)) + gaussian())
# (ci701-rev-reSE.R) random-effect and fitted standard errors of the
# k = 4 data from t0 = 500 and t0 = 2000 (both code 0 on ref)
for (t0 in c(500, 2000)) {
  fit <- suppressMessages(suppressWarnings(
    frm(f, data = mdata(4), start = list(theta = c(0, 0, t0)))))

  fe <- suppressWarnings(fitted(fit))
  cat(sprintf("RES t0=%d code=%d t=%.1f\n", t0, fit$opt$convergence,
              fit$estimates$theta[3]))
  cat("  vcov se:", sprintf("%.5f", sqrt(diag(suppressWarnings(vcov(fit))))),
      "\n")
  cat("  fitted EstErr:", sprintf("%.5f", as.numeric(fe[1:3, "Est.Error", ])), "\n")
}
