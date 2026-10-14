# Reviewer of lane postfix, round 2: what the scan's accepted fits of
# the sinusoid stress cases are (dev/postfix-rev2-nlnew.R stress block).
.libPaths(c("C:/Users/adf44/source/r/wt-postfix-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
f1 <- bf(y ~ A * sin(6.283185307179586 * x / P), A + P ~ 1, nl = TRUE)
for (cs in list(c(3, 1001), c(7, 1002), c(1.3, 1001), c(9, 1002))) {
  P0 <- cs[1]; s <- cs[2]
  set.seed(s)
  xs <- seq(0, 15, length.out = 120)
  ds <- data.frame(x = xs, y = 2 * sin(6.283185307179586 * xs / P0) +
                     rnorm(120, 0, 0.6))
  msgs <- character(); warns <- character()
  r <- withCallingHandlers(frm(f1, data = ds),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  h <- suppressWarnings(frm(f1, data = ds, start = list(beta = c(2, P0))))
  cat(sprintf("\nP %.1f seed %d: conv %d logLik %.4f (hand %.4f), warnings %d\n",
              P0, s, r$opt$convergence, as.numeric(logLik(r)),
              as.numeric(logLik(h)), length(warns)))
  print(fixef(r))
  cat(substr(grep("scan", msgs, value = TRUE), 1, 400), "\n")
  cat("sd(y)", sd(ds$y), "  residual sd at fit",
      sd(ds$y - fitted(r)[, 1]), "\n")
}
