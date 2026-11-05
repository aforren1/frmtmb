# Item 5: on the share mutant, what posterior_epred() returns by margin.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7701)
n <- 400
d <- data.frame(lat = rnorm(n))
d$sx <- 0.5
d$xo <- d$lat + rnorm(n, 0, d$sx)
d$y1 <- 1 + 1.5 * d$lat + rnorm(n, 0, 0.4)
d$y2 <- -1 + 1.2 * d$lat + rnorm(n, 0, 0.4)
fit <- quiet(frm(frmtmb::mvbf(bf(y1 ~ me(xo, sx)), bf(y2 ~ me(xo, sx)), rescor = FALSE),
                 family = gaussian(), data = d))
ds <- rev_draws_at(fit, 200L)
nd <- data.frame(xo = c(-1, 0, 1), sx = c(0.5, 1, 2))
set.seed(1)
fe <- frmtmb.sample::posterior_epred(ds, newdata = nd)
cat("dim", dim(fe), "\nsd by row, resp 1:", round(apply(fe[, , 1], 2, sd), 4), "\nsd by row, resp 2:", round(apply(fe[, , 2], 2, sd), 4), "\n")
cat("first rows resp 1:", round(fe[1:3, 1, 1], 4), " resp 2:", round(fe[1:3, 1, 2], 4), "\n")
cat("what the mv path evaluates:\n")
print(body(frmtmb.sample::posterior_epred))
one <- frmtmb.sample:::draws_fit_at(ds, 1L)
a <- frm_linpred(one, newdata = nd, resp = "y1")
b <- frm_linpred(one, newdata = nd, resp = "y2")
cat("one draw y1:", round(a, 4), " y2:", round(b, 4), "\n")
cat("beta:", round(fit$estimates$beta, 4), "\n")
cat("latent implied by y1:", round((a - fit$estimates$beta[1]) / fit$estimates$beta[2], 4),
    " by y2:", round((b - fit$estimates$beta[3]) / fit$estimates$beta[4], 4), "\n")
