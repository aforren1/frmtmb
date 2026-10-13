# Lane optrap: test-se-tier.R's RB2 case, as the test runs it.
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (th in c(-12, -20, -30)) for (s in c(36, 19)) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(100)
  w <- character()
  f <- withCallingHandlers(frm(y ~ x + (1 | g), data = d,
                               start = list(theta = th)),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")})
  cat(s, f$opt$objective, exp(f$opt$par[4]), "\n")
  print(substr(w, 1, 120))
}
