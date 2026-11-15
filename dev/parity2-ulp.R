# Lane parity2: the new warning assertion of test-brms-surface-mixed.R
# (attr(k, "min") = 0 names the flat parameters) over y * (1 + 2e-15 k),
# k = 0..24, data seed 9. Prints, per k, whether the SE warning fired and
# the frame's min and D.
#   Rscript dev/parity2-ulp.R > dev/parity2-out/ulp.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(9)
d <- data.frame(k = rep(1:4, 10))
d$y <- 0.5 * d$k + stats::rnorm(40)
attr(d$k, "min") <- 0
y0 <- d$y
hits <- 0L
for (k in 0:24) {
  d$y <- y0 * (1 + 2e-15 * k)
  w <- character()
  fit <- withCallingHandlers(
    suppressMessages(frm(bf(y ~ mo(k)), family = gaussian(), data = d)),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    })
  se <- any(grepl("Standard errors are not available", w, fixed = TRUE))
  hits <- hits + se
  m <- fit$frame$linpreds[[1L]]$mo[[1L]]
  cat(sprintf("k=%2d se_warning=%s min=%s D=%s other=%d\n", k, se, m$min,
              m$D, sum(!grepl("Standard errors are not available", w,
                              fixed = TRUE))))
}
cat("SE warning on", hits, "of 25\n")
