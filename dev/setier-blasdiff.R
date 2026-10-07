# Lane setier: the fits whose verdict differs between the BLAS builds in
# the lane's suite runs. Is it the fit that differs (the optimizer
# stopping elsewhere on a flat direction), or the verdict on the same
# point? Run with both Rscripts; compare the printed estimates.
#   Rscript dev/setier-blasdiff.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("BLAS probe", system.time({m <- matrix(1, 600, 600); m %*% m})[[3]],
    "\n")
show <- function(tag, fit) {
  lost <- ns$sdr_of(fit)$se_lost
  h <- fit$cache$hessian_fixed
  u <- fit$par_units %||% rep(1, length(fit$opt$par))
  cat(tag, "code", fit$opt$convergence, "logLik",
      format(as.numeric(logLik(fit)), digits = 12), "\n  par",
      format(fit$opt$par, digits = 6), "\n  lost",
      paste(names(lost), lost), "\n")
  if (!is.null(h)) {
    cat("  row max (opt units)",
        format(apply(abs(h$H * outer(u, u)), 1, max), digits = 3), "\n")
  }
}
q <- function(expr) suppressMessages(suppressWarnings(expr))
set.seed(79)
dd <- data.frame(t = factor(rep(1:4, 40)), g = factor(rep(1:20, each = 8)))
dd$y <- stats::rnorm(160, rep(stats::rnorm(20, 0, 0.7), each = 8), 1)
for (nm in c("cs", "homcs", "ar1", "hetar1")) {
  ff <- stats::as.formula(paste0("y ~ 1 + ", nm, "(t + 0 | g)"))
  show(paste("lkj", nm), q(frm(bf(ff) + gaussian(), data = dd)))
}
