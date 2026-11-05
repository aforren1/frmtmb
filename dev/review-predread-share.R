# Item 5: is sharing the latent "one line"? The lane says drop resp and dpar from the
# store's key. Run on the mutant me_share (REV_LIB=.../rev-predread-lib/me_share), which
# is exactly that line. Data seed 7701, as in review-predread-mejoint.R.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7701)
n <- 400
d <- data.frame(lat = rnorm(n))
d$sx <- 0.5
d$xo <- d$lat + rnorm(n, 0, d$sx)
d$y1 <- 1 + 1.5 * d$lat + rnorm(n, 0, 0.4)
d$y2 <- -1 + 1.2 * d$lat + rnorm(n, 0, 0.4)
d$ys <- 0.5 + 1.0 * d$lat + rnorm(n, 0, exp(-0.8 + 0.5 * d$lat))
fit <- quiet(frm(frmtmb::mvbf(bf(y1 ~ me(xo, sx)), bf(y2 ~ me(xo, sx)), rescor = FALSE),
                 family = gaussian(), data = d))
ds <- rev_draws_at(fit, 4000L)
nd <- data.frame(xo = c(-1, 0, 1), sx = c(0.5, 1, 2))
fmt <- function(v) paste(sprintf("%.4f", v), collapse = " ")
rowcor <- function(a) vapply(seq_len(dim(a)[2L]), function(j) cor(a[, j, 1], a[, j, 2]), 0)
set.seed(1)
cat("posterior_epred  cor by row:", fmt(rowcor(frmtmb.sample::posterior_epred(ds, newdata = nd))), "\n")
cat("posterior_predict cor by row:", fmt(rowcor(frmtmb.sample::posterior_predict(ds, newdata = nd))), "\n")
pm <- predict(fit, newdata = nd, summary = FALSE, ndraws = 4000, propagate_error = FALSE)
cat("predict(fit)     cor by row:", fmt(rowcor(pm)), "\n")
b1 <- fit$estimates$beta[grep("mexosx", names(fit$estimates$beta))[1L]]
b2 <- fit$estimates$beta[grep("mexosx", names(fit$estimates$beta))[2L]]
s_1 <- exp(fit$estimates$betad[1L])
s_2 <- exp(fit$estimates$betad[2L])
cat("target if shared, predictive:", fmt(b1 * b2 * nd$sx^2 / sqrt((b1^2 * nd$sx^2 + s_1^2) * (b2^2 * nd$sx^2 + s_2^2))), "\n")
cat("target if shared, epred: 1 1 1\n")
fs <- quiet(frm(bf(ys ~ me(xo, sx), sigma ~ me(xo, sx)), family = gaussian(), data = d))
dss <- rev_draws_at(fs, 4000L)
set.seed(2)
pp <- frmtmb.sample::posterior_predict(dss, newdata = nd)
cat("two dpars, predict skew by row:", fmt(apply(pp, 2, function(q) mean((q - mean(q))^3) / sd(q)^3)), "\n")
