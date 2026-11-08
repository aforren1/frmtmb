# Item 1: fitted() on a new random-effect level. Reference: Monte Carlo over
# the fit's own RE distribution at the estimates (u ~ N(0, Sigma_hat)), seed
# below, 2e6 draws.
source("dev/readers-setup.R")
set.seed(4101)
J <- 40; n <- 400
g <- factor(sample(seq_len(J), n, TRUE))
u <- rnorm(J, 0, 1)
x <- runif(n, -1, 1)
dd <- data.frame(x = x, g = g)
dd$y <- rpois(n, exp(0.5 + 0.3 * x + u[g]))
fit <- frm(bf(y ~ x + (1 | g)), family = poisson(), data = dd)
sdg <- exp(fit$estimates$theta[1])
vc <- VarCorr(fit)
cat("sd(g):", format(vc$g$sd[1, 1], digits = 6), "\n")
nd <- data.frame(x = c(-0.5, 0, 0.8), g = "new")
f <- fitted(fit, newdata = nd, allow_new_levels = TRUE)
print(f)
fe <- fixef(fit)[, "Estimate"]
eta <- fe[1] + fe[2] * nd$x
s <- vc$g$sd[1, 1]
cat("closed form exp(eta + s^2/2):", format(exp(eta + s^2 / 2), digits = 6),
    "\n")
cat("median exp(eta):", format(exp(eta), digits = 6), "\n")
uu <- rnorm(2e6, 0, s)
mc <- sapply(eta, function(e) c(mean = mean(exp(e + uu)),
                                sd = sd(exp(e + uu))))
cat("MC mean:", format(mc[1, ], digits = 6), "\n")
cat("MC sd over u alone:", format(mc[2, ], digits = 6), "\n")
cat("ratio closed/fitted:", format(exp(eta + s^2 / 2) / f[, "Estimate"],
                                   digits = 6), "\n")
# a random slope, correlated: the new level's eta variance is z' Sigma z
set.seed(4102)
u2 <- MASS::mvrnorm(J, c(0, 0), matrix(c(0.64, 0.2, 0.2, 0.36), 2))
dd$y2 <- rpois(n, exp(0.5 + 0.3 * x + u2[g, 1] + u2[g, 2] * x))
fit2 <- frm(bf(y2 ~ x + (1 + x | g)), family = poisson(), data = dd)
vc2 <- VarCorr(fit2)$g
S <- vc2$cov
if (is.null(S)) S <- diag(vc2$sd[, 1]^2)
if (length(dim(S)) == 3L) S <- S[, 1, ]
print(S)
f2 <- fitted(fit2, newdata = nd, allow_new_levels = TRUE)
print(f2)
fe2 <- fixef(fit2)[, "Estimate"]
eta2 <- fe2[1] + fe2[2] * nd$x
Z <- cbind(1, nd$x)
v2 <- rowSums((Z %*% S) * Z)
cat("closed form slope:", format(exp(eta2 + v2 / 2), digits = 6), "\n")
# bernoulli logit: no closed form; MC over u
set.seed(4103)
dd$yb <- rbinom(n, 1, plogis(0.2 + 0.5 * x + 1.5 * u[g]))
fit3 <- frm(bf(yb ~ x + (1 | g)), family = bernoulli(), data = dd)
s3 <- VarCorr(fit3)$g$sd[1, 1]
f3 <- fitted(fit3, newdata = nd, allow_new_levels = TRUE)
print(f3)
fe3 <- fixef(fit3)[, "Estimate"]
eta3 <- fe3[1] + fe3[2] * nd$x
uu3 <- rnorm(2e6, 0, s3)
cat("sd:", s3, " MC mean:", format(sapply(eta3, function(e)
  mean(plogis(e + uu3))), digits = 6), "\n")
cat("MC sd over u alone:", format(sapply(eta3, function(e)
  sd(plogis(e + uu3))), digits = 6), "\n")
# predict(): brms posterior_predict draws y given u; its mean is the same
p1 <- predict(fit, newdata = nd, allow_new_levels = TRUE)
print(p1)
