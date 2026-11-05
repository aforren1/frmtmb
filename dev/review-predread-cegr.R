# M4: brms on the calls of conditional_effects() that the lane refuses on me(gr = ).
# Data seed 7901, as in review-predread-misc.R.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7901)
n <- 300
d <- data.frame(z = rnorm(n), g = factor(rep(sprintf("g%02d", 1:30), each = 10)))
gl <- rnorm(30)[d$g]
d$sg <- 0.4
d$xg <- gl + rnorm(30, 0, 0.4)[d$g]
d$yg <- 1 + 1.2 * gl - 0.5 * d$z + rnorm(n, 0, 0.5)
G <- rev_pair("bf(yg ~ me(xg, sg, gr = g) + z)", "gaussian()", d, ndraws = 200)
ce <- brms::conditional_effects(G$fx, "xg", resolution = 4)[[1]]
cat("brms CE xg estimate:", signif(ce$estimate__, 5), "\n")
cat("brms grid g:", as.character(ce$g), " cond sg:", ce$sg[1], "\n")
cat("brms CE default names:", names(brms::conditional_effects(G$fx, resolution = 4)), "\n")
cat("brms CE posterior_predict:", signif(brms::conditional_effects(G$fx, "xg", method = "posterior_predict", resolution = 4)[[1]]$estimate__, 5), "\n")
cat("frmtmb coefficients:", signif(G$fit$estimates$beta, 5), "\n")
