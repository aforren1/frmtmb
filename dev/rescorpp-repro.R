## Reviewer: frmtmb.sample posterior_predict() and log_lik() of a
## set_rescor(TRUE) fit. Seed 3, 300 rows, residual correlation 0.8.
##   Rscript dev/review-mointerop-rescorpp.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-common.R")
set.seed(3); n <- 300
x <- rnorm(n); e <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .8, .8, 1), 2))
d <- data.frame(x = x, y1 = 1 + x + e[, 1], y2 = -1 + 0.5 * x + e[, 2])
f <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE), data = d)
show("rescor", rescor_matrix(f))
ds <- suppressWarnings(suppressMessages(frm_sample(
  f, chains = 1, iter = 600, seed = 1, refresh = 0)))
show("posterior_predict(ds) dim", dim(posterior_predict(ds)))
p1 <- posterior_predict(ds, resp = "y1")
p2 <- posterior_predict(ds, resp = "y2")
show("cor of draws across responses, rows 1..5",
     round(sapply(1:5, function(i) cor(p1[, i], p2[, i])), 3))
pc <- predict(f, summary = FALSE, ndraws = 300)
show("core predict cor rows 1..5",
     round(sapply(1:5, function(i) cor(pc[, i, 1], pc[, i, 2])), 3))
## log_lik(ds) against the bivariate normal density at each draw
ll <- show("log_lik(ds) dim", dim(log_lik(ds)))
L <- log_lik(ds)
dr <- as.data.frame(posterior::as_draws_df(ds))
show("draw columns", grep("^(b_|sigma|rescor)", names(dr), value = TRUE))
ref <- t(sapply(seq_len(nrow(L)), function(s) {
  m1 <- dr$b_y1_Intercept[s] + dr$b_y1_x[s] * x
  m2 <- dr$b_y2_Intercept[s] + dr$b_y2_x[s] * x
  s1 <- dr$sigma_y1[s]; s2 <- dr$sigma_y2[s]
  r <- dr$`rescor__y1__y2`[s]
  z1 <- (d$y1 - m1) / s1; z2 <- (d$y2 - m2) / s2
  -log(2 * pi * s1 * s2 * sqrt(1 - r^2)) -
    (z1^2 - 2 * r * z1 * z2 + z2^2) / (2 * (1 - r^2))
}))
show("log_lik(ds) vs bivariate normal: max rel diff",
     max(abs(L - ref) / abs(ref)))
