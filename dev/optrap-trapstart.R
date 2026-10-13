# Lane optrap: the 0.69.0 trap rebuilt on the lane by its start alone.
# The ri20 design with y x 1e-3 or x 1e3 and the group log sd started at
# 0 (an sd of 1 in the response's units, the old default): which fits
# still stop short, with what code, and what the fit says.
#   Rscript dev/optrap-trapstart.R <lib or ""> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:40
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", dirname(find.package("frmtmb")), "\n")
for (sc in c(1e-3, 1e3)) for (s in seeds) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- (1 + 0.5 * d$x + rnorm(100)) * sc
  w <- character(); m <- character()
  f <- withCallingHandlers(frm(y ~ x + (1 | g), data = d,
                               start = list(theta = 0)),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")},
    message = function(x) {m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")})
  l4 <- suppressMessages(lmer(y ~ x + (1 | g), data = d, REML = FALSE))
  dll <- as.numeric(logLik(f)) - as.numeric(logLik(l4))
  if (dll > -1e-4) next
  cat(sprintf("scale %g seed %d: code %d dll %.4f | %s | %s\n", sc, s,
              f$opt$convergence, dll, paste(substr(w, 1, 70), collapse = " ; "),
              paste(substr(m, 1, 40), collapse = " ; ")))
}
