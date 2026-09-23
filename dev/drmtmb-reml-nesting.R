# drmTMB REML integrates beta_sigma out only when sigma has a random
# effect (dev/drmtmb-log/reml-sigma.log). If so, sigma ~ z is not the
# sd -> 0 limit of sigma ~ z + (1 | g) under its REML criterion, and the
# larger model can report a LOWER restricted likelihood than the model
# it contains. frmtmb integrates beta_mu only in both, so it nests.
source("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-agree-lib.R")
one <- function(seed) {
  set.seed(seed)
  ng <- 30; n <- ng * 8
  g <- factor(rep(seq_len(ng), each = 8))
  x <- rnorm(n); z <- rnorm(n)
  # No group variation in sigma: the sigma random effect is truly zero.
  y <- 1 + 0.5 * x + rnorm(ng, 0, 0.6)[g] + rnorm(n, 0, exp(-0.2 + 0.3 * z))
  d <- data.frame(y, x, z, g)
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  dA <- q(drmTMB::drmTMB(dbf(y ~ x + (1 | g), sigma ~ z), data = d,
                         family = gaussian(), REML = TRUE))
  dB <- q(drmTMB::drmTMB(dbf(y ~ x + (1 | g), sigma ~ z + (1 | g)), data = d,
                         family = gaussian(), REML = TRUE))
  fA <- q(frm(bf(y ~ x + (1 | g), sigma ~ z), data = d, family = gaussian(),
              REML = TRUE))
  fB <- q(frm(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), data = d,
              family = gaussian(), REML = TRUE))
  dAml <- q(drmTMB::drmTMB(dbf(y ~ x + (1 | g), sigma ~ z), data = d,
                           family = gaussian()))
  dBml <- q(drmTMB::drmTMB(dbf(y ~ x + (1 | g), sigma ~ z + (1 | g)),
                           data = d, family = gaussian()))
  c(seed = seed,
    drm_REML_B_minus_A = as.numeric(logLik(dB)) - as.numeric(logLik(dA)),
    frm_REML_B_minus_A = as.numeric(logLik(fB)) - as.numeric(logLik(fA)),
    drm_ML_B_minus_A = as.numeric(logLik(dBml)) - as.numeric(logLik(dAml)),
    drm_REML_sd_sigma = exp(dB$opt$par[["log_sd_sigma"]]),
    frm_REML_sd_sigma = exp(fB$opt$par[[4]]),
    drm_REML_A_minus_frm_REML_A = as.numeric(logLik(dA)) -
      as.numeric(logLik(fA)))
}
out <- t(sapply(1:10, one))
print(signif(out, 6))
cat("\nseeds where drmTMB REML logLik(B) < logLik(A) - 1e-6:",
    sum(out[, "drm_REML_B_minus_A"] < -1e-6), "of", nrow(out), "\n")
cat("seeds where frmtmb REML logLik(B) < logLik(A) - 1e-6:",
    sum(out[, "frm_REML_B_minus_A"] < -1e-6), "of", nrow(out), "\n")
cat("seeds where drmTMB ML logLik(B) < logLik(A) - 1e-6:",
    sum(out[, "drm_ML_B_minus_A"] < -1e-6), "of", nrow(out), "\n")
