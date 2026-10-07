# Reviewer of lane setier, final check: dev/setier-rev2-window.R's expb
# seed 6, which the lane loses and base keeps. Is the direction real
# curvature? The exact Hessian, its smallest direction, the symmetric
# second differences at the probe's step and half of it, and an exact
# profile of the objective along the direction.
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(6)
d <- data.frame(x = 1 + runif(60) * 1e-5)
d$y <- 2 * exp(0.5 * d$x) + rnorm(60, 0, 1)
f <- suppressWarnings(frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
                          data = d, start = list(beta = c(2, 0.5))))
p <- f$opt$par
H <- f$obj$he(p)
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
cat("par", signif(p, 6), "| unit-diag ev", signif(e$values, 3), "\n")
k <- which.min(e$values)
v <- e$vectors[, k] / D
f0 <- f$obj$fn(p)
step <- sqrt(4e-3 / e$values[k])
for (m in c(0.125, 0.25, 0.5, 1, 2)) {
  cc <- f$obj$fn(p + m * step * v) + f$obj$fn(p - m * step * v) - 2 * f0
  cat(sprintf("step x %.3f: symmetric second difference %.4g (quadratic %.4g)\n",
              m, cc, 4e-3 * m^2))
}
cat("move along it (natural units) at the full step:", signif(step * v, 3),
    "\n")
cat("base-style SE (solve(H)):", signif(sqrt(diag(solve(H))), 3), "\n")
