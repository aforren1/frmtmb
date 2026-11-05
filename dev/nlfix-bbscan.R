# A beta_binomial body with no overdispersion in the data: which seeds send a
# restart to the numerical spike at a huge phi? Run on a build without
# nl_restart_credible(), to find the case its test needs.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
for (s in 1:12) {
  set.seed(s)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), ntr = 10L)
  d$y <- rbinom(n, 10, plogis(0.3 + 0.5 * d$x))
  m <- character(0)
  fit <- tryCatch(withCallingHandlers(
    frm(bf(y | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
        family = beta_binomial(), data = d),
    warning = function(c) invokeRestart("muffleWarning"),
    message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") }),
    error = function(e) NULL)
  cat(s, if (is.null(fit)) "ERR" else as.numeric(logLik(fit)),
      if (!is.null(fit)) fit$estimates$betad, length(m), "\n")
}
