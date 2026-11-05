# A mixture with mu1 ~ gp(x): the translator has no rule for brms's theta, so the
# joint is read against a hand-written density at frm()'s estimate instead, and the
# fit against brms's program by the same hand. Seed 20261102.
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261102)
n <- 80L
x <- round(runif(n, 0, 6), 3)
cl <- rbinom(n, 1, 0.5)
y <- ifelse(cl == 1, sin(2 * x), 3) + 0.2 * rnorm(n)
d <- data.frame(x = x, y = y)
fit <- suppressWarnings(frm(bf(y ~ 1, mu1 ~ gp(x), mu2 ~ 1), family = mixture(gaussian, gaussian), data = d))
print(fit$opt$convergence)
est <- fit$estimates
print(names(est))
print(est$beta); print(est$betad); print(est$theta)
bk <- fit$frame$re_blocks[[1L]]
th <- est$theta[bk$theta_idx]
u <- sort(unique(x)); J <- match(x, u)
K <- exp(-outer(u, u, "-")^2 / (2 * exp(2 * th[2])))
diag(K) <- diag(K) + 1e-12
f <- exp(th[1]) * as.vector(t(chol(K)) %*% est$b[bk$b_idx])
dp <- frmtmb:::eval_dpars(fit)[[1L]]
print(names(dp))
ours <- -as.numeric(fit$obj$env$f(fit$obj$env$last.par.best))
m1 <- dp$mu1; m2 <- dp$mu2; s1 <- dp$sigma1; s2 <- dp$sigma2; t1 <- dp$theta1
t1 <- plogis(t1)
lp <- sum(log(t1 * dnorm(y, m1, s1) + (1 - t1) * dnorm(y, m2, s2)))
hand <- lp + sum(dnorm(est$b[bk$b_idx], log = TRUE))
cat(sprintf("NCGPM ours %.8f hand %.8f diff %.3e; field vs mu1 max diff %.3e\n", ours, hand, ours - hand,
            max(abs(m1 - (mean(est$beta[1]) + f[J])))))
