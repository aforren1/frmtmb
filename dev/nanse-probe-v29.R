# The two suite fixtures whose assertions the lane changes, examined:
# test-review-v29.R's boundary correlation and test-bcm-mpt.R's MPT
# covariance block. Prints the lost set, the unit-diagonal spectrum,
# the eigenvectors of the removed directions, and base's sdreport SEs.
#   Rscript dev/nanse-probe-v29.R [lib]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/wt-nanse-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
show <- function(fit) {
  ns <- asNamespace("frmtmb")
  nm <- ns$outer_par_names(fit)
  sdr <- ns$sdr_of(fit)
  cat("lost:", paste(names(sdr$se_lost), sdr$se_lost, sep = "=",
                     collapse = " "), "\n")
  base <- suppressWarnings(RTMB::sdreport(fit$obj))
  cat("base SE:", paste(sprintf("%s=%.3g", nm,
                                suppressWarnings(sqrt(diag(base$cov.fixed)))),
                        collapse = " "), "\n")
  H <- ns$fit_outer_hessian(fit)$H
  D <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(D, D), symmetric = TRUE)
  cat("eig:", paste(signif(e$values, 3), collapse = " "), "\n")
  for (k in which(e$values < 1e-3 * max(e$values))) {
    v <- e$vectors[, k]
    cat(sprintf("  dir %d (%.3g):", k, e$values[k]),
        paste(sprintf("%s=%.3f", nm, v)[abs(v) > 0.01], collapse = " "),
        "\n")
  }
}
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
fit <- suppressWarnings(
  frm(mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
           bf(y2 ~ x + (1 | q | id)) + gaussian()), data = dd))
cat("== review-v29 boundary correlation, code", fit$opt$convergence, "\n")
show(fit)
