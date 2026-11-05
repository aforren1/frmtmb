# R7: why marginaleffects' SE differs from fitted()'s on a new row, on the review's
# fit (dev/review-latres-m6b.R, seed 167): its delta method moves get_coef(), which
# holds beta and not theta. Compares with the closed form exp(b0 + b1 x + v / 2).
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(warn = -1)
set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:30))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(30, arima.sim(list(ar = .5), 5, sd = .5)))
d$y <- rpois(n, exp(.8 + .4 * d$x + e))
f <- suppressMessages(frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = poisson(), data = d))
nd <- data.frame(t = 1:3, g = factor(rep("new", 3)), x = c(-1, 0, 1))
fv <- fitted(f, newdata = nd)
me <- marginaleffects::predictions(f, newdata = nd)
th <- f$estimates$theta
v <- exp(2 * th[1]) * (1 + th[2]^2)
V <- as.matrix(vcov(f, full = TRUE))
X <- cbind(1, nd$x)
m <- as.numeric(exp(X %*% f$estimates$beta + v / 2))
pb <- 1:2
pt <- match(c("theta_1", "theta_2"), rownames(V))
beta_only <- sapply(1:3, function(i) m[i] * sqrt(drop(t(X[i, ]) %*% V[pb, pb] %*% X[i, ])))
full <- sapply(1:3, function(i) {
  g <- numeric(nrow(V))
  g[pb] <- m[i] * X[i, ]
  g[pt] <- m[i] * c(v, exp(2 * th[1]) * th[2])
  sqrt(drop(t(g) %*% V %*% g))
})
print(data.frame(fitted = fv[, 1], predictions = me$estimate, se_fitted = fv[, 2],
                 se_closed_form = full, se_marginaleffects = me$std.error,
                 se_beta_only = beta_only))
print(names(frmtmb:::get_coef.frmtmb_fit(f)))
