# Lane ci701: what the clamp pull-back (se_edge_cor()) changes on fits
# that already lost only their correlation at 0.70.0: y ~ x + (1 + x | g)
# with no slope variance (the shape of dev/ci701-vc-spread.R), seeds
# 1..20. Prints the correlation parameter, the lost set and every
# standard error vcov() shows, so the base and the lane build can be
# compared line by line.
#
#   Rscript dev/ci701-rs20-se.R <lib|base> [n_seeds]
args <- commandArgs(TRUE)
lib <- args[1]
ns <- if (length(args) > 1) as.integer(args[2]) else 20L
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (lib == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
for (s in seq_len(ns)) {
  set.seed(s)
  ng <- 20
  g <- factor(rep(seq_len(ng), each = 10))
  d <- data.frame(g = g, x = stats::rnorm(200))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(ng, 0, 0.7)[g] + stats::rnorm(200)
  fit <- suppressMessages(suppressWarnings(
    frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d)))
  lost <- frmtmb:::sdr_of(fit)$se_lost
  se <- suppressWarnings(sqrt(diag(vcov(fit, full = TRUE))))
  cat(sprintf("RES seed=%d t=%.4g lost=%s se=%s\n", s,
              fit$estimates$theta[3],
              if (length(lost)) paste(names(lost), lost, sep = ":",
                                      collapse = ",") else "none",
              paste(sprintf("%.6g", se), collapse = ",")))
}
