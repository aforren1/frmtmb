# Reviewer probe (lane support): frmtmb.survival's cox() under a
# non-binding trunc(lb <= 0), the delayed-entry form with entry at 0.
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.survival)
})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(403)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- rexp(n, exp(0.5 * d$x)) + 0.01
d$cc <- as.numeric(runif(n) < 0.2)
d$lb0 <- 0
d$lbm <- -1
fitq <- function(fo) tryCatch(suppressWarnings(frm(fo, family = cox(), data = d)),
                              error = function(e) conditionMessage(e))
f0 <- fitq(y | cens(cc) ~ x)
for (fo in list(y | cens(cc) + trunc(lb = lb0) ~ x,
                y | cens(cc) + trunc(lb = lbm) ~ x)) {
  f1 <- fitq(fo)
  if (is.character(f1) || is.character(f0)) {
    cat(deparse(fo), "ERR", substr(if (is.character(f1)) f1 else f0, 1, 200), "\n")
    next
  }
  cat(deparse(fo), " dlogLik", as.numeric(logLik(f1)) - as.numeric(logLik(f0)),
      " dbeta", max(abs(fixef(f1)[, 1] - fixef(f0)[, 1])), "\n")
}
