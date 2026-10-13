# Lane seflags, item 1 under the reference BLAS: `y ~ x + s(x)`, whose
# smooth's fixed column is exactly a + b x, so the alias is (Intercept,
# x, s(x).fx1). The mean of x sets a, and so the intercept's share of
# the flat direction: below sqrt(se_flat_tol) the intercept keeps its
# standard error, and 0.69.0's prediction test read the direction
# without that share. Prints the lost set, the alias, and the warnings
# of fitted(newdata) for each mean.
#   Rscript dev/seflags-alias.R <lib> [re]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
re <- length(a) > 1 && a[2] == "re"
cat("lib", find.package("frmtmb"), "| R", R.home(), "| re", re, "\n")
ns <- asNamespace("frmtmb")
for (m in c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)) {
  set.seed(1)
  n <- 200
  x <- stats::rnorm(n)
  x <- x - mean(x) + m
  d <- data.frame(x = x, g = factor(rep(1:10, n / 10)))
  d$y <- 1 + 0.5 * d$x + sin(2 * d$x) + stats::rnorm(10, 0, 0.5)[d$g] +
    stats::rnorm(n, 0, 0.3)
  fo <- if (re) y ~ x + s(x) + (1 | g) else y ~ x + s(x)
  w <- character()
  fit <- withCallingHandlers(frm(fo, data = d),
    warning = function(z) { w <<- c(w, conditionMessage(z))
      invokeRestart("muffleWarning") },
    message = function(z) invokeRestart("muffleMessage"))
  lost <- suppressWarnings(ns$sdr_of(fit))$se_lost
  X <- fit$frame$linpreds[[1]]$X
  cn <- colnames(X)
  ifx <- grep("^s[(]x[)][.]?fx", cn)
  cf <- coef(stats::lm(as.numeric(X[, ifx[1]]) ~ d$x))
  nd <- d[c(1, 50, 100, 150, 200), ]
  pw <- character()
  pr <- withCallingHandlers(stats::fitted(fit, newdata = nd, re_formula = NA),
    warning = function(z) { pw <<- c(pw, conditionMessage(z))
      invokeRestart("muffleWarning") })
  cat(sprintf(paste0("mean %-6g code %d | fx1 = %.3g + %.4g x | lost %s |",
                     " fit warnings %d | fitted(newdata): %s | se %s\n"),
              m, fit$opt$convergence, cf[1], cf[2],
              paste(names(lost), lost, sep = ":", collapse = ","), length(w),
              if (length(pw)) substr(pw[1], 1, 40) else "silent",
              paste(signif(pr[, "Est.Error"], 3), collapse = " ")))
}
