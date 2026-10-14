# Lane ci701: does confint_varcorr() lose every interval of a block when
# the standard-error check removes one of its parameters? VarCorr() reads
# the covariance through hyp_prop_var(), which gives NaN only to rows
# that move along a lost direction (lane nanse review, RB2);
# varcorr_trans_rows() reads sdr$cov.fixed, whose lost rows are NaN, and
# a product with it spreads NaN to every row of the block.
#
# Case: y ~ x + (1 + x | g) with no slope variance (the rs20 shape of
# dev/setier-rev-rs20.R), seeds 1..20, against VarCorr()'s standard
# errors of the same fit.
#
#   Rscript dev/ci701-vc-spread.R <lib|base> [n_seeds]
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
  bm <- 0L
  fit <- withCallingHandlers(
    frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d),
    frmtmb_boundary_fit = function(m) {
      bm <<- bm + 1L; invokeRestart("muffleMessage")
    },
    warning = function(w) invokeRestart("muffleWarning"))
  lost <- frmtmb:::sdr_of(fit)$se_lost
  ci <- suppressWarnings(confint_varcorr(fit))
  i1 <- which(ci$term == "(Intercept)" & ci$type == "sd")
  cat(sprintf(
    "RES seed=%d code=%d bmsg=%d lost=%s int_sd=%.4f [%.4f, %.4f]\n", s,
    fit$opt$convergence, bm,
    if (length(lost)) paste(names(lost), lost, sep = ":", collapse = ",")
    else "none",
    ci$estimate[i1], ci$lwr[i1], ci$upr[i1]))
}
