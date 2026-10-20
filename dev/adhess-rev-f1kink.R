# Re-check of fix round 1, finding 1: sdr_nan_nonpos(). Per fit: which
# variances are NaN, whether every entry of vcov() outside the NaN rows
# is bitwise the raw inverse's, zero-width confint() rows, warnings that
# escape vcov(), confint() and summary(), and the lost reasons.
#   REV_ARM=lane Rscript dev/adhess-rev-f1kink.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-common.R")
mk <- list(
  al_free = function(d) frm(y ~ x, family = asym_laplace(), data = d),
  al_q50 = function(d) frm(bf(y ~ x, quantile = 0.5),
                           family = asym_laplace(), data = d),
  al_q90 = function(d) frm(bf(y ~ x, quantile = 0.9),
                           family = asym_laplace(), data = d),
  zi_al = function(d) {
    d$y[seq(1, nrow(d), by = 10)] <- 0
    frm(y ~ x, family = zero_inflated_asym_laplace(), data = d)
  },
  huber = function(d) frm(y ~ x, family = huber(), data = d),
  huber_k0.01 = function(d) frm(y ~ x, family = huber(k = 0.01), data = d))
for (nm in names(mk)) for (seed in 1:5) {
  set.seed(seed)
  n <- 400
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n)
  r <- cap(mk[[nm]](d))
  fit <- r$v
  if (inherits(fit, "rev_err")) {
    cat(nm, seed, "ERROR", fit, "\n")
    next
  }
  sdr <- sdr_of(fit)
  raw <- sdr$cov_fixed_raw %||% sdr$cov.fixed
  rv <- cap(vcov(fit, full = TRUE))
  V <- rv$v
  bad <- is.na(diag(V))
  same <- identical(unname(V[!bad, !bad]), unname(raw[!bad, !bad]))
  prop <- sdr$cov_fixed_prop
  rc <- cap(confint(fit))
  ci <- rc$v
  zw <- if (inherits(ci, "rev_err")) NA else
    sum(ci[, "upr"] - ci[, "lwr"] <= 0, na.rm = TRUE)
  rs <- cap(utils::capture.output(summary(fit)))
  sline <- grep("nocurv|curvature|not available", rs$v, value = TRUE)
  lost <- sdr$se_lost
  cat(sprintf(paste0("%-11s s%d code %d | NaN %d of %d | raw diag %s | ",
                     "kept entries identical %s | ci zero-width %s | lost %s",
                     " | warn fit %d vcov %d confint %d summary %d\n"),
              nm, seed, fit$opt$convergence, sum(bad), length(bad),
              paste(signif(diag(raw), 3), collapse = "/"), same, zw,
              if (length(lost)) paste(names(lost), lost, sep = "=",
                                      collapse = ",") else "-",
              length(r$w), length(rv$w), length(rc$w), length(rs$w)))
  if (length(c(rv$w, rc$w, rs$w))) {
    cat("   escaped:", substr(c(rv$w, rc$w, rs$w), 1, 120), sep = "\n   ")
  }
  if (seed == 1 && length(sline)) cat("   summary:", substr(sline, 1, 150),
                                      sep = "\n   ")
}
