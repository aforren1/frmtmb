# Review of lane illse2: the REML log-likelihood shift on the lane's own
# REML designs (test-se-reparam-reml.R's data), against lme4 and 0.71.0.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-rev-remlll2.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
reml_data <- function(seed, deg) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.5)[d$g] +
    stats::rnorm(200, 0, 0.3)
  for (k in seq_len(deg)) d[[paste0("x", k)]] <- d$x^k
  d
}
ldT <- function(f) {
  bw <- f$beta_whiten
  if (is.null(bw)) return(0)
  sum(vapply(bw$blocks, function(b) {
    determinant(b$T, logarithm = TRUE)$modulus[[1]]
  }, 0))
}
for (deg in 3:6) {
  d <- reml_data(1, deg)
  fo <- stats::as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                                "+ (1 | g)"))
  f <- suppressWarnings(frm(fo, data = d, REML = TRUE))
  l <- suppressWarnings(lme4::lmer(fo, data = d, REML = TRUE))
  cat(sprintf("deg %d %s | REML logLik frm %.8f lme4 %.8f | diff %+.6f | log|det T| %+.6f | AIC frm %.4f\n",
              deg, arm, as.numeric(logLik(f)), as.numeric(logLik(l)),
              as.numeric(logLik(f)) - as.numeric(logLik(l)), ldT(f),
              AIC(f)))
}
