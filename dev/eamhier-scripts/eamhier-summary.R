# Lane eamhier: read the per-replicate records and report the table.
#
# Every number in dev/eamhier-findings.md comes out of this script, and
# every coverage carries a Wilson interval, because 60 replicates put a
# standard error of about 2.8 points on a coverage near 95 and a bare
# percentage invites a reader to distinguish 91.7 from 95.0.
#
# Run:  Rscript --vanilla dev/eamhier-scripts/eamhier-summary.R <dir>...

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
source("dev/eamhier-scripts/eamhier-summary-read.R")

# The coverage AND the two numbers that say what a miss is made of: the
# ratio of the mean reported standard error to the spread of the
# estimates over replicates (below 1 is a narrow interval), and the bias
# in units of that spread.
cover_line <- function(tag, lo, hi, truth, est = NULL) {
  ok <- is.finite(lo) & is.finite(hi)
  k <- sum(lo[ok] <= truth & hi[ok] >= truth)
  n <- sum(ok)
  w <- wilson(k, n)
  extra <- ""
  if (!is.null(est)) {
    e <- est[ok]
    se <- (hi[ok] - lo[ok]) / (2 * stats::qnorm(0.975))
    extra <- sprintf("  se/sd %5.3f  bias/sd %6.3f",
                     mean(se) / stats::sd(e),
                     (mean(e) - truth) / stats::sd(e))
  }
  cat(sprintf("  %-22s %3d/%3d = %5.1f%%  (%4.1f, %5.1f)  truth %-9g%s\n",
              tag, k, n, 100 * k / n, 100 * w[1L], 100 * w[2L], truth,
              extra))
  invisible(c(k = k, n = n))
}

# The Monte Carlo error on a mean over replicates. A bias is only a bias
# if it is larger than this.
est_line <- function(tag, x, truth = NA_real_) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) {
    cat(sprintf("  %-22s (none)\n", tag))
    return(invisible(NULL))
  }
  m <- mean(x)
  se <- stats::sd(x) / sqrt(n)
  if (is.finite(truth)) {
    cat(sprintf(
      "  %-22s mean %10.5f  mcse %8.5f  sd %8.5f  bias %9.5f  z %6.2f  n %d\n",
      tag, m, se, stats::sd(x), m - truth, (m - truth) / se, n))
  } else {
    cat(sprintf(
      "  %-22s mean %10.5f  mcse %8.5f  sd %8.5f  min %9.4f  max %9.4f  n %d\n",
      tag, m, se, stats::sd(x), min(x), max(x), n))
  }
  invisible(m)
}

d <- read_records(args)
# `ndt` is drawn as 0.25 * exp(u), u normal with standard deviation
# 0.12, so 0.25 is the MEDIAN of the per-subject non-decision times and
# their mean is 0.25 * exp(0.12^2 / 2). A population estimate is scored
# against the mean, and getting that wrong would read as a 2 ms bias
# that is arithmetic rather than estimation.
tr <- list(mu0 = 0.4, mu_cond = 0.9, bs = 1.4, ndt = 0.25,
           ndt_mean = 0.25 * exp(0.12^2 / 2),
           sd_mu = 0.35, sd_lbs = 0.20, sv = 0.4)

key <- paste(d$arm, d$bound, d$ns, d$nt, sep = "/")
for (k in unique(key)) {
  x <- d[key == k, , drop = FALSE]
  cat("\n== ", k, "  n = ", nrow(x), " replicates ==\n", sep = "")
  bad <- x$status != "ok"
  cat(sprintf("  fits: ok %d, error %d | conv!=0 %d | pdHess FALSE %d",
              sum(!bad), sum(bad), sum(x$conv != 0, na.rm = TRUE),
              sum(!(x$pdHess == "TRUE"), na.rm = TRUE)))
  cat(sprintf(" | any NaN se %d | flat %d | extreme_theta %d\n",
              sum(x$n_bad_se > 0, na.rm = TRUE),
              sum(x$n_flat > 0, na.rm = TRUE),
              sum(x$n_extreme_theta > 0, na.rm = TRUE)))
  cat(sprintf("  diagnose unbounded_dpar: %s\n",
              paste(names(table(x$unbounded)), table(x$unbounded),
                    sep = ":", collapse = " ")))
  cat(sprintf("  fit seconds: median %.1f, max %.1f\n",
              stats::median(x$fit_s, na.rm = TRUE),
              max(x$fit_s, na.rm = TRUE)))
  ok <- x$status == "ok"
  x <- x[ok, , drop = FALSE]
  if (!nrow(x)) next

  cat(" -- estimates --\n")
  est_line("mu0", x$mu0, if (x$arm[1L] == "B") 1.1 else tr$mu0)
  if (any(is.finite(x$mu_cond))) {
    est_line("mu_cond", x$mu_cond, tr$mu_cond)
  }
  est_line("bs = exp(lbs)", exp(x$lbs), tr$bs)
  if (any(is.finite(x$lsv))) est_line("sv = exp(lsv)", exp(x$lsv), tr$sv)
  vcn <- strsplit(x$vc_names[1L], ";", fixed = TRUE)[[1L]]
  cat("   variance components, in VarCorr order: ",
      paste(vcn, collapse = " | "), "\n", sep = "")
  for (j in seq_along(vcn)) {
    est_line(paste0("sd", j, " (", vcn[j], ")"), x[[paste0("sd", j)]])
  }
  est_line("ndt population, s", x$ndt_hat,
           if (x$arm[1L] == "B") 0.25 else tr$ndt_mean)
  est_line("ndt bias, ms", x$ndt_bias_ms, 0)
  est_line("ndt per-subject rmse ms", x$ndt_rmse_ms)
  est_line("ndt spread hat, ms", x$ndt_sd_hat_ms)
  est_line("ndt spread true, ms", x$ndt_sd_true_ms)
  est_line("ndt corr with truth", x$ndt_corr)
  est_line("logLik", x$logLik)
  est_line("min margin to floor ms", x$ndt_margin_min_ms)
  cat(sprintf("  subjects below their own floor: %d of %d replicates ",
              sum(x$n_below_own_floor == x$ns), nrow(x)))
  cat("have all subjects inside\n")

  cat(" -- Wald coverage, nominal 95 --\n")
  cover_line("mu0", x$mu0_lo, x$mu0_hi,
             if (x$arm[1L] == "B") 1.1 else tr$mu0, x$mu0)
  if (any(is.finite(x$mu_cond))) {
    cover_line("mu_cond", x$mu_cond_lo, x$mu_cond_hi, tr$mu_cond,
               x$mu_cond)
  }
  cover_line("log bs", x$lbs_lo, x$lbs_hi, log(tr$bs), x$lbs)
  if (any(is.finite(x$lsv))) {
    cover_line("log sv", x$lsv_lo, x$lsv_hi, log(tr$sv), x$lsv)
  }
  if (x$arm[1L] != "B") {
    cover_line("sd(mu|s)", x$sd1_lo, x$sd1_hi, tr$sd_mu, x$sd1)
    cover_line("sd(log bs|s)", x$sd2_lo, x$sd2_hi, tr$sd_lbs, x$sd2)
  } else {
    cover_line("sd(log bs|s)", x$sd1_lo, x$sd1_hi, 0.25, x$sd1)
  }
  # NOT a calibrated population interval, and it is reported to say so:
  # `ndt_hat` is the fitted population FRACTION times the mean of the
  # floors the data produced, and its standard error carries the
  # fraction's uncertainty only. The per-subject error below is the
  # quantity item 1.0a's note says to read.
  cover_line("ndt population",
             x$ndt_hat - stats::qnorm(0.975) * x$ndt_se,
             x$ndt_hat + stats::qnorm(0.975) * x$ndt_se,
             if (x$arm[1L] == "B") 0.25 else tr$ndt_mean)
  # against the mean of the draw this replicate actually used, which is
  # a different estimand from the population value and moves with the
  # 30 subjects drawn
  lo <- x$ndt_hat - stats::qnorm(0.975) * x$ndt_se
  hi <- x$ndt_hat + stats::qnorm(0.975) * x$ndt_se
  k <- sum(lo <= x$ndt_pop_true & hi >= x$ndt_pop_true)
  w <- wilson(k, nrow(x))
  cat(sprintf(
    "  %-22s %3d/%3d = %5.1f%%  (%4.1f, %5.1f)  the draw's own mean\n",
    "ndt vs drawn mean", k, nrow(x), 100 * k / nrow(x),
    100 * w[1L], 100 * w[2L]))

  # What a miss is MADE of. z is the error over the standard error the
  # fit reported: sd(z) near 1 with a matching robust spread and a
  # kurtosis near 3 says the interval is the right shape and the wrong
  # width, which is a different defect from a heavy tail. This is the
  # decomposition ?lca uses for the same shape of finding.
  cat(" -- the standardized error, z = (est - truth) / reported se --\n")
  zline <- function(tag, est, lo, hi, truth) {
    ok <- is.finite(lo) & is.finite(hi)
    se <- (hi[ok] - lo[ok]) / (2 * stats::qnorm(0.975))
    z <- (est[ok] - truth) / se
    n <- length(z)
    kur <- mean((z - mean(z))^4) / stats::sd(z)^4
    cat(sprintf(
      paste0("  %-22s sd(z) %5.3f  IQR/1.349 %5.3f  kurtosis %5.2f",
             "  |z|>1.96 %d of %d\n"),
      tag, stats::sd(z), stats::IQR(z) / 1.349, kur,
      sum(abs(z) > stats::qnorm(0.975)), n))
  }
  zline("mu0", x$mu0, x$mu0_lo, x$mu0_hi,
        if (x$arm[1L] == "B") 1.1 else tr$mu0)
  if (any(is.finite(x$mu_cond))) {
    zline("mu_cond", x$mu_cond, x$mu_cond_lo, x$mu_cond_hi, tr$mu_cond)
  }
  zline("log bs", x$lbs, x$lbs_lo, x$lbs_hi, log(tr$bs))
  if (any(is.finite(x$lsv))) {
    zline("log sv", x$lsv, x$lsv_lo, x$lsv_hi, log(tr$sv))
  }
  if (x$arm[1L] != "B") {
    zline("sd(mu|s)", x$sd1, x$sd1_lo, x$sd1_hi, tr$sd_mu)
    zline("sd(log bs|s)", x$sd2, x$sd2_lo, x$sd2_hi, tr$sd_lbs)
  }

  # What the SHIPPED tier asserts on this row, over the same
  # replicates: tests/testthat/test-scale.R expects
  # scale_z(ndt_hat, ndt_se, 0.25) < 4, at one seed. A tier assertion
  # that holds on 59 of 60 seeds is a seed-dependent assertion.
  if (x$arm[1L] != "B") {
    z <- abs(x$ndt_hat - 0.25) / x$ndt_se
    cat(sprintf(
      "  tier assertion z < 4 on ndt: %d of %d pass, z median %.2f, max %.2f\n",
      sum(z < 4), length(z), stats::median(z), max(z)))
    z2 <- abs(x$ndt_hat - x$ndt_pop_true) / x$ndt_se
    cat(sprintf(
      "  the same z against the draw's own mean: %d of %d under 4, max %.2f\n",
      sum(z2 < 4), length(z2), max(z2)))
    # The replacement this lane proposes for that assertion: the
    # per-subject error against the between-subject spread the same run
    # measured. A ratio under 1 says the fit tells the subjects apart.
    rr <- x$ndt_rmse_ms / x$ndt_sd_true_ms
    cat(sprintf(
      "  per-subject rmse / true spread: %d of %d under 1, max %.3f\n",
      sum(rr < 1), length(rr), max(rr)))
  }
}
