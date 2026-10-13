# mo() with a factor and in group-level terms: fits, names, predictions.
# Usage: Rscript dev/moparity-mo-smoke.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args) && nzchar(args[[1]])) args[[1]] else
  "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(1)
n <- 240
dd <- data.frame(
  x = sample(0:3, n, TRUE),
  f = factor(sample(c("a", "b", "c"), n, TRUE)),
  g = factor(sample(paste0("g", 1:8), n, TRUE)),
  z = rnorm(n), w = rnorm(n)
)
eff <- c(0, 1, 1.5, 1.7)
bg <- rnorm(8, 0, 0.5)
dd$y <- 1 + (1 + (dd$f == "b") - (dd$f == "c") + bg[dd$g]) * eff[dd$x + 1] +
  rnorm(n, 0, 0.5)
show <- function(fo) {
  cat("\n==", deparse1(fo), "\n")
  r <- tryCatch({
    fit <- frm(bf(fo), family = gaussian(), data = dd)
    print(round(fixef(fit), 4))
    cat("variables:", variables(fit), "\n")
    if (length(fit$frame$re_blocks)) print(VarCorr(fit))
    p1 <- fitted(fit)[1:3, 1]
    p2 <- fitted(fit, newdata = dd[1:3, ])[, 1]
    cat("fitted in sample vs newdata, max diff:", max(abs(p1 - p2)), "\n")
    cat("logLik", format(as.numeric(logLik(fit)), digits = 10), "\n")
    fit
  }, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
  invisible(r)
}
show(y ~ mo(x) * f)
show(y ~ mo(x):f)
show(y ~ f + mo(x):f)
show(y ~ z + mo(x):z)
show(y ~ mo(x) * z * w)
show(y ~ mo(x) + (mo(x) | g))
show(y ~ mo(x) + (1 + mo(x) || g))
show(y ~ mo(x) * z + (mo(x) | g))
show(y ~ (mo(x) | g))
show(y ~ mo(x) * f + (mo(x) | g))
