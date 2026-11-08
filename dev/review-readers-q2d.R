# Review of lane readers: fitted()'s Q columns at a two-direction new row
# (me_band_dir(): the quantiles along the direction the value moves most,
# widened to the whole variance) against the quantiles of the mean over the
# fit's own draws. Hurdle (1 | p | g) on mu and hu as case 4 of
# dev/review-readers-shapes.R (data seed 5204), and the same row with the
# correlation's sign forced by a larger hu sd (data seed 5205).
# Monte Carlo seed 61, 1e6 draws per row.
source("dev/review-readers-setup.R")
run <- function(seed, S) {
  set.seed(seed)
  J <- 30; n <- 600
  d <- data.frame(x = runif(n, -1, 1), g = factor(sample(J, n, TRUE)))
  U <- MASS::mvrnorm(J, c(0, 0), S)
  lam <- exp(0.6 + 0.4 * d$x + U[d$g, 1])
  cnt <- rpois(n, lam)
  while (any(cnt == 0)) cnt[cnt == 0] <- rpois(sum(cnt == 0), lam[cnt == 0])
  d$y <- ifelse(runif(n) < plogis(-0.4 + U[d$g, 2]), 0L, cnt)
  f <- frm(bf(y ~ x + (1 | p | g), hu ~ (1 | p | g)),
           family = hurdle_poisson(), data = d)
  vc <- VarCorr(f)$g
  sd2 <- vc$sd[, 1]
  rho <- vc$cor[1, 1, 2]
  Sh <- diag(sd2) %*% matrix(c(1, rho, rho, 1), 2) %*% diag(sd2)
  nd <- data.frame(x = c(-0.5, 0.5), g = "new")
  ff <- fitted(f, newdata = nd, allow_new_levels = TRUE)
  em <- as.numeric(frm_linpred(f, newdata = nd, re_formula = NA))
  eh <- as.numeric(frm_linpred(f, newdata = nd, re_formula = NA,
                               dpar = "hu"))
  set.seed(61)
  Z <- MASS::mvrnorm(1e6, c(0, 0), Sh)
  cat("sd", round(sd2, 3), "rho", round(rho, 3), "\n")
  for (i in 1:2) {
    l <- exp(em[i] + Z[, 1])
    v <- (1 - plogis(eh[i] + Z[, 2])) * l / (-expm1(-l))
    q <- quantile(v, c(0.025, 0.975))
    cat(sprintf(paste("row %d: Est %.4f MC %.4f | Q2.5 %.4f MC %.4f | Q97.5",
                      "%.4f MC %.4f | coverage of [Q2.5, Q97.5] under MC %.3f\n"),
                i, ff[i, "Estimate"], mean(v), ff[i, "Q2.5"], q[1],
                ff[i, "Q97.5"], q[2],
                mean(v >= ff[i, "Q2.5"] & v <= ff[i, "Q97.5"])))
  }
}
run(5204, matrix(c(0.64, -0.3, -0.3, 0.49), 2))
run(5205, matrix(c(0.36, 0.4, 0.4, 1.44), 2))
