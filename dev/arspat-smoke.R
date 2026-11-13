# Lane arspat: quick fits of the lifted combinations on the lane build.
# Output: dev/arspat-out/smoke.txt
lib <- Sys.getenv("ARSPAT_LIB", "C:/Users/adf44/source/r/wt-arspat-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(dirname(system.file(package = "frmtmb"))), "\n")
chk <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", conditionMessage(w), "\n")
           invokeRestart("muffleWarning")
         }),
       error = function(e) { cat("ERROR: ", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) print(r)
  invisible(r)
}
set.seed(2)
g <- rep(1:6, each = 8)
t <- rep(1:8, 6)
d <- data.frame(g = factor(g), t = t, x = rnorm(48))
d$y <- d$x + rnorm(48)
d$s <- runif(48, 0.2, 0.5)
chk("se + ar cov=TRUE", {
  f <- frm(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE), data = d)
  c(logLik = as.numeric(logLik(f)), f$opt$par)
})
chk("se(sigma=FALSE) + ar cov=TRUE",
    frm(y | se(s) ~ x + ar(t, g, cov = TRUE), data = d))
chk("se + cosy, student", {
  f <- frm(y | se(s, sigma = TRUE) ~ x + cosy(t, g), data = d,
           family = student())
  c(logLik = as.numeric(logLik(f)), f$opt$par)
})
d$ym <- d$y
d$ym[c(3, 20, 41)] <- NA
chk("mi + ar cov=TRUE", {
  f <- frm(ym | mi() ~ x + ar(t, g, cov = TRUE), data = d)
  c(logLik = as.numeric(logLik(f)), f$opt$par)
})
chk("mi + ar cov=TRUE, rows dropped by hand", {
  f <- frm(ym ~ x + ar(t, g, cov = TRUE), data = d)
  c(logLik = as.numeric(logLik(f)), f$opt$par)
})
chk("mi + fcor", {
  M <- 0.5^abs(outer(1:48, 1:48, "-"))
  f <- frm(ym | mi() ~ x + fcor(M), data = d, data2 = list(M = M))
  c(logLik = as.numeric(logLik(f)), f$opt$par)
})
