# Lane adhess, fix round 1: families whose log density is piecewise
# linear, or has no curvature, in some parameter. Per case: the
# variances vcov() shows (how many are 0, negative, NaN), whether any
# confint() row has zero width, the se_lost reasons, and the warnings.
# Seeds 1 to 3.
#   ADHESS_ARM=lane|base Rscript dev/adhess-kinksweep.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
mk <- list(
  asym_laplace_free = function(d) frm(y ~ x, family = asym_laplace(),
                                      data = d),
  asym_laplace_q50 = function(d) frm(bf(y ~ x, quantile = 0.5),
                                     family = asym_laplace(), data = d),
  asym_laplace_q90 = function(d) frm(bf(y ~ x, quantile = 0.9),
                                     family = asym_laplace(), data = d),
  zi_asym_laplace = function(d) {
    d$y[seq(1, nrow(d), by = 10)] <- 0
    frm(y ~ x, family = zero_inflated_asym_laplace(), data = d)
  },
  huber = function(d) frm(y ~ x, family = huber(), data = d),
  huber_k0.01 = function(d) frm(y ~ x, family = huber(k = 0.01), data = d),
  poisson_all0 = function(d) {
    d$c <- 0L
    frm(c ~ x, family = poisson(), data = d)
  },
  bernoulli_all1 = function(d) {
    d$b <- 1L
    frm(b ~ x, family = bernoulli(), data = d)
  },
  cens_all_right = function(d) {
    d$cen <- "right"
    frm(y | cens(cen) ~ x, data = d)
  }
)
for (nm in names(mk)) for (seed in 1:3) {
  set.seed(seed)
  n <- 400
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n)
  r <- cap(mk[[nm]](d))
  fit <- r$v
  if (is_err(fit)) {
    cat(sprintf("%-18s s%d ERROR %s\n", nm, seed, substr(fit, 1, 100)))
    next
  }
  rv <- cap(vcov(fit, full = TRUE))
  v <- diag(rv$v)
  rc <- cap(confint(fit))
  ci <- rc$v
  zw <- if (is_err(ci)) NA else
    sum(ci[, "upr"] - ci[, "lwr"] <= 0, na.rm = TRUE)
  lost <- ns$sdr_of(fit)$se_lost
  cat(sprintf(paste0("%-18s s%d code %d | var 0 %d, <0 %d, NaN %d of %d | ",
                     "ci zero-width %s | lost %s | warn fit %s, vcov %s, ",
                     "confint %s\n"),
              nm, seed, fit$opt$convergence, sum(v == 0, na.rm = TRUE),
              sum(v < 0, na.rm = TRUE), sum(is.na(v)), length(v), zw,
              if (length(lost)) paste(unique(lost), collapse = ",") else "-",
              wtag(r$w), wtag(rv$w), wtag(rc$w)))
}
