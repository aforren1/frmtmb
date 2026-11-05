# sar() with a body, an extension's residual term: the body's maximum and
# fitted values against the linear spelling, and what brms builds. Seed 9101.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
set.seed(9101)
n <- 100
W <- matrix(0, n, n)
for (i in 2:n) W[i, i - 1] <- W[i - 1, i] <- 1
W <- W / rowSums(W)
d <- data.frame(x = rnorm(n))
e <- solve(diag(n) - 0.4 * W, rnorm(n, 0, 0.5))
d$y <- 2 + 0.8 * d$x + e
fl <- frm(bf(y ~ x, autocor = ~ sar(W)), data = d, data2 = list(W = W))
fn <- frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE, autocor = ~ sar(W)),
          data = d, data2 = list(W = W),
          start = list(beta = c("b0_(Intercept)" = 2, "b1_(Intercept)" = 0.5)))
cat("logLik linear, nl:", format(as.numeric(logLik(fl)), digits = 10),
    format(as.numeric(logLik(fn)), digits = 10), "\n")
cat("max |fitted difference|:",
    max(abs(fitted(fl)[, "Estimate"] - fitted(fn)[, "Estimate"])), "\n")
r <- tryCatch(nchar(suppressMessages(brms::stancode(
  brms::bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE, autocor = ~ sar(W)),
  data = d, data2 = list(W = W)))), error = function(e) conditionMessage(e))
cat("brms stancode:", r, "\n")
