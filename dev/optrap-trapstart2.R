# Lane optrap: a short fit at a nonzero optimizer code, built by the
# start (ri20 design, unit-scale response, group log sd started far
# below the response's unit).
#   Rscript dev/optrap-trapstart2.R <lib or "">
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", dirname(find.package("frmtmb")), "\n")
for (th in c(-7, -12, -20)) for (s in c(19, 22, 36, 37)) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(100)
  w <- character(); m <- character()
  f <- withCallingHandlers(frm(y ~ x + (1 | g), data = d,
                               start = list(theta = th)),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")},
    message = function(x) {m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")})
  l4 <- suppressMessages(lmer(y ~ x + (1 | g), data = d, REML = FALSE))
  dll <- as.numeric(logLik(f)) - as.numeric(logLik(l4))
  cat(sprintf("theta %g seed %d: code %d (%s) dll %.4f | %s | %s\n", th, s,
              f$opt$convergence, f$opt$message, dll,
              paste(substr(w, 1, 70), collapse = " ; "),
              paste(substr(m, 1, 40), collapse = " ; ")))
}
