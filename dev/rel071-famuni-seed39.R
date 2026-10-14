# The one row of lane famuni's 720-fit check that moves on the merged
# build: identity link, xi0 -0.99, n 1000, seed 39 (dev/famuni-r2-xiedge.R's
# seed rule). Its code is 1 with a warning on the lane and 0 on the
# merged build, at the same xi. Is that lane optstart's restart rule?
#   Rscript dev/rel071-famuni-seed39.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")), "\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
xi0 <- -0.99; n <- 1000; s <- 39
set.seed(30000 + 1000 * n + s + round(-xi0 * 100))
y <- rg(n, 0, 1, xi0)
ws <- character(0)
f <- withCallingHandlers(
  frm(y ~ 1, family = gen_extreme_value(link_xi = "identity"),
      data = data.frame(y = y)),
  warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") },
  message = function(m) invokeRestart("muffleMessage"))
cat("code", f$opt$convergence, f$opt$message, "\n")
cat("objective", format(f$opt$objective, digits = 15), "\n")
cat("par", format(f$opt$par, digits = 10), "\n")
if (!is.null(f$opt$restart)) {
  cat("restart:", paste(names(f$opt$restart), unlist(f$opt$restart), collapse = "; "), "\n")
}
cat("warnings:", length(ws), "\n"); for (w in ws) cat("  ", substr(w, 1, 120), "\n")
