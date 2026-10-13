# Lane ci691: where nlminb stops on the test's xbeta phi 2e4 fit when the
# start moves by a few ulps, a stand-in for the runners' different libm
# rounding. Prints, per perturbation, the code, iterations, estimates,
# logLik and whether every standard error is finite.
#   Rscript dev/ci691-xbeta-perturb.R <lib> <phi>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r7"
phi <- if (length(args) > 1) as.numeric(args[2]) else 2e4
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r7",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "phi", phi, "\n")
set.seed(4200)
z <- stats::rbeta(1000, 0.047 * phi, 0.953 * phi)
d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
u <- suppressWarnings(frm(y ~ 1, family = xbeta(), data = d,
                          dry_run = "objective"))
p0 <- u$obj$par
cat("default start", format(p0, digits = 17), "\n")
for (k in 0:24) {
  st <- p0
  if (k > 0) {
    set.seed(k)
    st <- p0 * (1 + sample(c(-4:-1, 1:4), length(p0), TRUE) *
                  .Machine$double.eps)
  }
  ws <- character(0)
  fit <- withCallingHandlers(
    frm(y ~ 1, family = xbeta(), data = d, start = list(
      beta = unname(st[1]), betad = unname(st[2:3]))),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, error = function(e) cat("ERROR", conditionMessage(e), "\n"))
  se <- sqrt(diag(vcov(fit, full = TRUE)))
  cat(sprintf("k %2d code %d it %4d par %s ll %.6f se_finite %s warn %s\n",
              k, fit$opt$convergence, fit$opt$iterations,
              paste(sprintf("%9.4f", fit$opt$par), collapse = " "),
              as.numeric(logLik(fit)), all(is.finite(se)),
              paste(substr(ws, 1, 40), collapse = " | ")))
}
