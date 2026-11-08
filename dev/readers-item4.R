# Item 4 as filed: a latent on mu, a me() value on hu alone, hurdle_poisson.
# Reference: Monte Carlo over the two draws at the estimates, 2e6 draws, seed
# 11. Data seed 9204.
source("dev/readers-setup.R")
ns <- asNamespace("frmtmb")
set.seed(9204)
d <- expand.grid(t = 1:25, s = factor(1:8))
d$x <- rnorm(nrow(d)); d$zo <- rnorm(nrow(d)); d$sz <- 0.5
e <- ave(rnorm(nrow(d)), d$s, FUN = function(z) stats::filter(z, 0.5,
                                                               "recursive"))
zl <- d$zo + rnorm(nrow(d), 0, 0.5)
lam <- exp(0.4 + 0.3 * d$x + 0.4 * e)
cnt <- rpois(nrow(d), lam)
while (any(cnt == 0)) { i <- cnt == 0; cnt[i] <- rpois(sum(i), lam[i]) }
d$y <- ifelse(runif(nrow(d)) < plogis(-0.5 + 1.2 * zl), 0L, cnt)
f <- suppressWarnings(frm(bf(y ~ x + ar(time = t, gr = s, cov = TRUE),
                             hu ~ me(zo, sz)), data = d,
                          family = hurdle_poisson()))
nd <- data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1), zo = c(-1, 1.5),
                 sz = 0.5)
print(fitted(f, newdata = nd))
lp <- function(nd, dp) as.numeric(frm_linpred(f, newdata = nd, dpar = dp,
                                              type = "link"))
em <- lp(nd, "mu"); eh <- lp(nd, "hu")
bh <- lp(transform(nd, zo = zo + 1), "hu") - eh
sl <- ns$latres_row_sd(f, "y", NULL)
set.seed(11)
z1 <- rnorm(2e6); z2 <- rnorm(2e6)
mc <- sapply(1:2, function(i) {
  l <- exp(em[i] + sl * z1)
  v <- (1 - plogis(eh[i] + bh[i] * 0.5 * z2)) * l / (1 - exp(-l))
  c(mean = mean(v), sd = sd(v))
})
cat("latent sd", sl, " me slope on hu", bh[1], "\n")
cat("MC mean:", format(mc[1, ], digits = 6), " MC sd (draws):",
    format(mc[2, ], digits = 6), "\n")
mm <- ns$fitted_me_moments(f, nd, NULL, NULL, NULL, FALSE)
cat("route mean:", format(mm$mean, digits = 6), " route draw sd:",
    format(sqrt(mm$var), digits = 6), "\n")
