# Reviewer of lane postfix: the data the brms and frmtmb.sample arms of
# the conditional_effects(method = "predict") comparison share.
set.seed(20261013)
n <- 100
xp <- runif(n, -1, 1)
dp <- data.frame(x = xp, y = rpois(n, exp(1 + 0.8 * xp)))
xb <- runif(n, -1, 1)
db <- data.frame(x = xb, y = rbinom(n, 1, plogis(0.3 + 1.2 * xb)))
