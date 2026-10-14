# Lane illse: the xbeta phi 2e4 fit from perturbed start k (dev/illse-
# xbeta-perturb.R), its warnings in full and the unit-diagonal spectrum
# of the AD Hessian at the stop.
#   Rscript dev/illse-xbeta-k.R <k> [lane|base|dev]
args <- commandArgs(TRUE)
k <- as.integer(args[1])
arm <- if (length(args) > 1) args[2] else "dev"
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-illse-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "dev") {
  suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-illse",
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
phi <- 2e4
set.seed(4200)
z <- stats::rbeta(1000, 0.047 * phi, 0.953 * phi)
d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
u <- suppressWarnings(frm(y ~ 1, family = xbeta(), data = d,
                          dry_run = "objective"))
p0 <- u$obj$par
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
  })
cat("code", fit$opt$convergence, "\n")
cat(ws, sep = "\n---\n")
H <- fit$obj$he(fit$opt$par)
D <- sqrt(diag(H))
ev <- eigen(H / outer(D, D), TRUE, only.values = TRUE)$values
cat("\nspectrum", format(ev, digits = 4), "ratio", format(min(ev) / max(ev), digits = 4), "\n")
