.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib", "C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1); d <- data.frame(x = rnorm(50)); d$y <- 1 + d$x + (rexp(50)^(-0.1) - 1) / 0.1
for (f in list(bf(y ~ x), bf(y ~ x, xi ~ 1))) {
  fit <- frm(f, family = gen_extreme_value(), data = d)
  lp <- fit$frame$linpreds[[grep("xi", names(fit$frame$linpreds))]]
  cat("----", deparse(f$formula), length(f$pforms), "\n")
  str(lp[c("dpar", "center", "constant")], max.level = 1)
  cat("X cols:", colnames(lp$X), "\n")
  dp <- fit$spec$responses[[1]]$dpars[["xi"]]
  cat("spec dpars xi names:", names(dp), "\n")
  print(setdiff(names(lp), c("X")))
}
