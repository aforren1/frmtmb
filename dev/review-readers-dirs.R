# Review of lane readers: rows with three or more correlated directions, on
# rev-readers-lib (REV_ARM=mutant), which is the lane's source with switches
# read from REV_MUT at run time; unset, it is the lane.
# A: a 3-predictor row whose three me() calls each reach two predictors
#    (chol_rows()'s cross term is live): the lane, the chol_cross mutant, and
#    Monte Carlo (seed 51, 2e6 draws). Data seed 6201.
# B: a mixture of two gaussians with five predictors, each reading its own
#    me() call: the frmtmb_draw_dims warning, how often it fires per call,
#    and what the truncation costs against Monte Carlo (seed 52, 2e6 draws).
#    Data seed 6202.
source("dev/review-readers-setup.R")
dm <- function(f, nd, dp) as.numeric(frm_linpred(f, newdata = nd, dpar = dp))
slope <- function(f, nd, dp, v) {
  nd1 <- nd
  nd1[[v]] <- nd1[[v]] + 1
  dm(f, nd1, dp) - dm(f, nd, dp)
}

cat("== A. hurdle_negbinomial, mu ~ me1 + me2, shape ~ me2 + me3, hu ~ me1 + me3\n")
set.seed(6201)
n <- 2000
d <- data.frame(l1 = rnorm(n), l2 = rnorm(n), l3 = rnorm(n), s1 = 0.3,
                s2 = 0.3, s3 = 0.3)
for (k in 1:3) d[[paste0("x", k)]] <- d[[paste0("l", k)]] + rnorm(n, 0, 0.3)
mu <- exp(0.7 + 0.6 * d$l1 + 0.5 * d$l2)
sh <- exp(1 + 0.4 * d$l2 - 0.4 * d$l3)
cnt <- rnbinom(n, mu = mu, size = sh)
while (any(cnt == 0)) {
  i <- cnt == 0
  cnt[i] <- rnbinom(sum(i), mu = mu[i], size = sh[i])
}
d$y <- ifelse(runif(n) < plogis(-0.3 + 0.6 * d$l1 + 0.6 * d$l3), 0L, cnt)
fa <- frm(bf(y ~ me(x1, s1) + me(x2, s2), shape ~ me(x2, s2) + me(x3, s3),
             hu ~ me(x1, s1) + me(x3, s3)),
          family = hurdle_negbinomial(), data = d)
nd <- data.frame(x1 = c(0, 1), x2 = c(0, -0.5), x3 = c(0.5, 0), s1 = 1.2,
                 s2 = 1.2, s3 = 1.2)
B <- sapply(c("mu", "shape", "hu"), function(dp) {
  sapply(c("x1", "x2", "x3"), function(v) slope(fa, nd[1, ], dp, v))
})
print(round(B, 4))
Sys.setenv(REV_MUT = "")
lane <- ns$fitted_me_moments(fa, nd, NULL, NULL, NULL, FALSE)
Sys.setenv(REV_MUT = "chol_cross")
mut <- ns$fitted_me_moments(fa, nd, NULL, NULL, NULL, FALSE)
Sys.setenv(REV_MUT = "")
set.seed(51)
N <- 2e6
mc <- t(vapply(seq_len(nrow(nd)), function(i) {
  Z <- matrix(rnorm(3 * N), N) * 1.2
  em <- dm(fa, nd[i, ], "mu") + Z %*% B[, "mu"]
  es <- dm(fa, nd[i, ], "shape") + Z %*% B[, "shape"]
  eh <- dm(fa, nd[i, ], "hu") + Z %*% B[, "hu"]
  m <- exp(em); k <- exp(es)
  v <- (1 - plogis(eh)) * m / (-expm1(-k * log1p(m / k)))
  c(mean(v), sd(v) / sqrt(N), sd(v))
}, numeric(3)))
cat("row  MC mean (se)          lane mean  rel        mutant mean rel        MC sd     lane sd   mutant sd\n")
for (i in seq_len(nrow(nd))) {
  cat(sprintf("%d   %.5f (%.5f)  %.5f %9.2e  %.5f %9.2e   %.4f   %.4f   %.4f\n",
              i, mc[i, 1], mc[i, 2], lane$mean[i], lane$mean[i] / mc[i, 1] - 1,
              mut$mean[i], mut$mean[i] / mc[i, 1] - 1, mc[i, 3],
              sqrt(lane$var[i]), sqrt(mut$var[i])))
}

cat("\n== B. mixture(gaussian, gaussian), five predictors with me()\n")
set.seed(6202)
n <- 3000
d <- data.frame(matrix(rnorm(5 * n), n))
names(d) <- paste0("l", 1:5)
for (k in 1:5) {
  d[[paste0("x", k)]] <- d[[paste0("l", k)]] + rnorm(n, 0, 0.3)
  d[[paste0("s", k)]] <- 0.3
}
th2 <- plogis(0.2 + 0.8 * d$l5)
comp2 <- runif(n) < th2
d$y <- ifelse(comp2,
              3 + 0.8 * d$l2 + rnorm(n, 0, exp(-0.3 + 0.3 * d$l4)),
              -2 + 0.8 * d$l1 + rnorm(n, 0, exp(-0.5 + 0.3 * d$l3)))
fb <- tryCatch(frm(bf(y ~ 1, mu1 ~ me(x1, s1), mu2 ~ me(x2, s2),
                      sigma1 ~ me(x3, s3), sigma2 ~ me(x4, s4),
                      theta2 ~ me(x5, s5)),
                   family = mixture(gaussian(), gaussian()), data = d),
               error = function(e) {
                 cat("fit error:", conditionMessage(e), "\n")
                 NULL
               })
if (!is.null(fb)) {
  cat("dpars in order:", names(fb$spec$responses[[1]]$dpars), "\n")
  nd <- data.frame(x1 = 0, x2 = 0, x3 = 0, x4 = 0, x5 = 0.5,
                   s1 = 1, s2 = 1, s3 = 1, s4 = 1, s5 = 2)
  r <- with_warn(fitted(fb, newdata = nd))
  cat("fitted():", r$value[, "Estimate"], " Est.Error",
      r$value[, "Est.Error"], "\n")
  cat("warnings per fitted() call:", length(r$warnings), " classes:",
      paste(wclass(r$warnings), collapse = " ; "), "\n")
  r2 <- with_warn(conditional_effects(fb, effects = "x1",
                                      int_conditions = list(x1 = 0),
                                      conditions = nd[, -1]))
  cat("warnings per conditional_effects() call:", length(r2$warnings), "\n")
  dps <- c("mu1", "mu2", "sigma1", "sigma2", "theta2")
  vs <- c(mu1 = "x1", mu2 = "x2", sigma1 = "x3", sigma2 = "x4",
          theta2 = "x5")
  sx <- c(1, 1, 1, 1, 2)
  e0 <- sapply(dps, function(p) dm(fb, nd, p))
  bb <- sapply(dps, function(p) slope(fb, nd, p, vs[[p]]))
  cat("slopes:", round(bb, 3), "\n")
  set.seed(52)
  N <- 2e6
  Z <- matrix(rnorm(5 * N), N)
  et <- sweep(Z, 2, bb * sx, `*`)
  et <- sweep(et, 2, e0, `+`)
  t2 <- plogis(et[, 5])
  v <- (1 - t2) * et[, 1] + t2 * et[, 2]
  cat("MC mean over all five draws:", mean(v), "+-", sd(v) / sqrt(N),
      " sd", sd(v), "\n")
  cat("value with theta2's draw at 0:",
      mean((1 - plogis(e0[5])) * et[, 1] + plogis(e0[5]) * et[, 2]), "\n")
}
