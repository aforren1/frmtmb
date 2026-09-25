# Is predict()'s cs() defect newdata-only, or does predict() drop cs()
# offsets entirely? Compare predict() in sample with fitted() (exact).
#   RV_LIB=base Rscript dev/simnewdata-review/rv-cs-insample.R
source("dev/simnewdata-review/rv-prelude.R")
set.seed(5)
n <- 300
d <- data.frame(x = rnorm(n))
d$yo <- factor(cut(2 * d$x + rlogis(n), c(-Inf, -0.5, 0.5, Inf),
                   labels = FALSE), ordered = TRUE)
fo <- frm(bf(yo ~ cs(x)), family = sratio(), data = d)
nd <- 20000
p <- unclass(predict(fo, ndraws = nd))[1:3, ]
f <- fitted(fo)
ex <- sapply(1:3, function(k) f[1:3, "Estimate", k])
cat("x of rows 1..3:", format(d$x[1:3], digits = 3), "\n")
cat("exact (fitted()):\n"); print(round(ex, 4))
cat("predict(), in sample:\n"); print(round(p, 4))
cat("observed marginal freq:", format(prop.table(table(d$yo)), digits = 3), "\n")
cat(sprintf("max |z| in-sample predict vs exact: %.1f\n",
            max(abs(p - ex) / sqrt(ex * (1 - ex) / nd))))
s <- as.matrix(simulate(fo, nsim = nd, seed = 1))[1:3, ]
fr <- t(apply(s, 1, function(r) table(factor(r, levels = 1:3)) / nd))
cat(sprintf("max |z| simulate() in sample vs exact: %.1f\n",
            max(abs(fr - ex) / sqrt(ex * (1 - ex) / nd))))
