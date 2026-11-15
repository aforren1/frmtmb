# Lane selimit2, backlog item 2829: vcov() and summary() under
# check_se = "ignore" on the rs20 design (y ~ x + (1 + x | g), slope sd
# 0), seeds 7 and 12. Prints every warning each call gives.
#   Rscript dev/selimit2-2829.R <lib or base> [seeds]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else c(7, 12)
ws <- function(label, expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) invokeRestart("muffleMessage"))
  cat(sprintf("  %-10s warnings %d %s\n", label, length(w),
              paste(substr(w, 1, 90), collapse = " | ")))
  v
}
for (s in seeds) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(120)
  cat("seed", s, "\n")
  f <- ws("frm", frm(y ~ x + (1 + x | g), data = d,
                     control = frmtmb_control(check_se = "ignore")))
  cat("  code", f$opt$convergence, "theta", signif(f$estimates$theta, 4),
      "\n")
  v <- ws("vcov", vcov(f))
  print(sqrt(diag(v)))
  s1 <- ws("summary", summary(f))
  v2 <- ws("vcov full", vcov(f, full = TRUE))
  print(diag(v2))
}
