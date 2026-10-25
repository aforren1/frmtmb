.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1); n <- 9; W <- matrix(1/8, n, n); diag(W) <- 0
d <- data.frame(y = rnorm(n), x = rnorm(n))
for (f in c("y ~ x + sar(W)", "y ~ x + sar(W, type = 'error')", "y ~ x + bf(y ~ x)")) {
  if (grepl("bf", f)) next
  r <- tryCatch({frm(as.formula(f), data = d, data2 = list(W = W)); "FIT"},
                error = function(e) conditionMessage(e))
  cat(f, "->", r, "\n")
}
r <- tryCatch({frm(bf(y ~ x, autocor = ~ sar(W)), data = d, data2 = list(W = W)); "FIT"},
              error = function(e) conditionMessage(e))
cat("bf(autocor = ~ sar(W)) ->", r, "\n")
cat("frmtmb.spatial loaded:", isNamespaceLoaded("frmtmb.spatial"), "\n")
suppressMessages(library(frmtmb.spatial))
r <- tryCatch({f <- frm(y ~ x + sar(W), data = d, data2 = list(W = W)); paste("FIT rho", format(confint_varcorr(f)$estimate))},
              error = function(e) conditionMessage(e))
cat("after library(frmtmb.spatial):", r, "\n")
