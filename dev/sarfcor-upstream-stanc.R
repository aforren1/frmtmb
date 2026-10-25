# Upstream repro (dev/upstream-bugs.md): brms 2.23.0 writes Stan code
# for sar() with weights(), se(), cens() or trunc() that Stan's own
# parser rejects, so the model fails at compile time rather than with a
# message. The plain model type-checks. Brms and rstan only, no C++
# compile. Rscript dev/sarfcor-upstream-stanc.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages({library(brms); library(rstan)})
cat("brms", format(packageVersion("brms")), "rstan",
    format(packageVersion("rstan")), "\n")
k <- 4
g <- expand.grid(r = 1:k, c = 1:k)
n <- nrow(g)
A <- outer(1:n, 1:n, function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- A / rowSums(A)
set.seed(5)
d <- data.frame(y = rnorm(n), x = rnorm(n), w = runif(n, 0.5, 1.5),
                se = 0.3, cen = 0L, lb = -10)
for (f in list(y ~ x + sar(W), y | weights(w) ~ x + sar(W),
               y | se(se, sigma = TRUE) ~ x + sar(W),
               y | cens(cen) ~ x + sar(W), y | trunc(lb = lb) ~ x + sar(W))) {
  code <- stancode(f, data = d, data2 = list(W = W))
  r <- tryCatch({stanc(model_code = code, verbose = FALSE); "type-checks"},
                error = function(e) {
                  paste("stanc ERROR:",
                        substr(gsub("\n", " ", conditionMessage(e)), 1, 200))
                })
  cat(deparse1(f), ":", r, "\n")
}
