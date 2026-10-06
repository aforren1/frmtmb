# Lane fixes, punch round: is the boundary runaway of a zero smooth SD
# driven by the scale of the null-space column? Seed 3, s(x1) + s(x2):
# the lane's frame with the two .fx columns multiplied by c, refitted
# through the internal fitting core (an exact reparameterization).
#   Rscript dev/fixes-sx-conv6.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
set.seed(3)
d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, scale = 2,
                                   verbose = FALSE))
fo <- bf(y ~ s(x1) + s(x2))
spec <- frm(fo, data = d, dry_run = "spec")
fr0 <- frm(fo, data = d, dry_run = "frame")
fx <- grep("[.]fx", colnames(fr0$linpreds[[1]]$X))
cat("fx column sds:", apply(fr0$linpreds[[1]]$X[, fx], 2, sd), "\n")
for (cc in c(0.04, 0.2, 1, 5, 25)) {
  fr <- fr0
  fr$linpreds[[1]]$X[, fx] <- fr$linpreds[[1]]$X[, fx] * cc
  fit <- suppressWarnings(frmtmb:::fit_assembled(
    spec, fr, frmtmb:::as_bform(fo, NULL), quote(frm()), REML = FALSE, start = NULL,
    control = frmtmb_control(), se = FALSE, lower = NULL, upper = NULL,
    prior = NULL, quadrature = FALSE))
  se <- suppressWarnings(fixef(fit)[, "Est.Error"])
  cat(sprintf("c = %5.2f: conv %d iters %3d min theta %8.2f finite SE %s logLik %.8f\n",
              cc, fit$opt$convergence, fit$opt$iterations,
              min(fit$estimates$theta), all(is.finite(se)),
              as.numeric(logLik(fit))))
}
