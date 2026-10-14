# Lane illse: dev/ci691-xbeta-perturb.R on the lane library or
# rellib-r8. Where nlminb stops on the xbeta phi 2e4 fit when the start
# moves by a few ulps, and the SE verdict there: the code, the lost set,
# whether every SE is finite and the warnings.
#   Rscript dev/illse-xbeta-perturb.R <lane|base|dev> <phi>
args <- commandArgs(TRUE)
lib <- if (length(args) && args[1] == "lane") "C:/Users/adf44/source/r/wt-illse-lib" else "C:/Users/adf44/source/r/rellib-r8"
phi <- if (length(args) > 1) as.numeric(args[2]) else 2e4
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r8",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
if (length(args) && args[1] == "dev") {
  # the worktree as it is, without an install (pure R)
  suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-illse",
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
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
  # vcov() may deliver a deferred check; its warnings belong to the fit
  se <- withCallingHandlers(sqrt(diag(vcov(fit, full = TRUE))),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  lost <- frmtmb:::sdr_of(fit)$se_lost
  cat(sprintf(paste0("k %2d code %d it %4d par %s ll %.6f se_finite %s ",
                     "lost %s warn %s\n"),
              k, fit$opt$convergence, fit$opt$iterations,
              paste(sprintf("%9.4f", fit$opt$par), collapse = " "),
              as.numeric(logLik(fit)), all(is.finite(se)),
              if (length(lost)) paste(names(lost), lost, sep = ":",
                                      collapse = ",") else "-",
              paste(substr(ws, 1, 70), collapse = " | ")))
}
