# Lane eamhier: the single-level control for the `sv` question.
#
# Same 12,000 rows, same two conditions, same truths, NO hierarchy. If
# `sv` and the condition effect recover here and not in arm C, the
# hierarchy is implicated; if they miss here too, the design is.
#
# Two things are reported that a coverage alone does not carry. The
# ratio of the mean reported standard error to the spread of the
# estimates over replicates says whether a miss is a BIAS or a narrow
# interval. The correlation between the estimated `sv` and the
# estimated condition effect says whether the two trade off, which is
# the mechanism dev/ndt-findings.md proposed from one seed.
#
# Run:  Rscript --vanilla dev/eamhier-scripts/eamhier-summary-S1.R <dir>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
source("dev/eamhier-scripts/eamhier-summary-read.R")

d <- read_records(args)
d <- d[d$status == "ok", , drop = FALSE]
cat("single-level sv control: ", nrow(d), " replicates, n = ",
    d$n[1L], " rows, sv true ", d$sv_true[1L], "\n", sep = "")
cat("  seeds ", min(d$seed), " to ", max(d$seed), "\n", sep = "")
cat(sprintf("  convergence 0: %d | pdHess TRUE: %d | any NaN se: %d\n",
            sum(d$conv == 0), sum(d$pdHess == "TRUE"),
            sum(d$n_bad_se > 0)))
cat(sprintf("  fit seconds: median %.1f, max %.1f\n",
            stats::median(d$fit_s), max(d$fit_s)))

row <- function(tag, est, lo, hi, truth) {
  n <- length(est)
  se <- (hi - lo) / (2 * stats::qnorm(0.975))
  k <- sum(lo <= truth & hi >= truth)
  w <- wilson(k, n)
  cat(sprintf(
    paste0("  %-12s mean %8.4f  bias %8.4f (mcse %6.4f)  sd %7.4f",
           "  mean se %7.4f  se/sd %5.3f  cover %2d/%2d = %5.1f%%",
           " (%4.1f, %5.1f)\n"),
    tag, mean(est), mean(est) - truth, stats::sd(est) / sqrt(n),
    stats::sd(est), mean(se), mean(se) / stats::sd(est), k, n,
    100 * k / n, 100 * w[1L], 100 * w[2L]))
}

cat(" -- estimates and Wald coverage, nominal 95 --\n")
row("mu0", d$mu0, d$mu0_lo, d$mu0_hi, 0.4)
row("mu_cond", d$mu_cond, d$mu_cond_lo, d$mu_cond_hi, 0.9)
row("log bs", d$lbs, d$lbs_lo, d$lbs_hi, log(1.4))
row("log sv", d$lsv, d$lsv_lo, d$lsv_hi, log(d$sv_true[1L]))
row("ndt, s", d$ndt_hat, d$ndt_hat - stats::qnorm(0.975) * d$ndt_se,
    d$ndt_hat + stats::qnorm(0.975) * d$ndt_se, 0.25)
cat(sprintf(
  "  sv on the response scale: mean %.4f, median %.4f, %.4f to %.4f\n",
  mean(exp(d$lsv)), stats::median(exp(d$lsv)),
  min(exp(d$lsv)), max(exp(d$lsv))))

# A COLLAPSED fit: `sv` has run to the log link's floor and the Wald
# interval it reports is not a statement about anything. The rule is a
# ratio to what the run itself measured, never a constant: a reported
# standard error more than five times the median one over the same
# replicates. Such a fit "covers" whatever it is asked about, so the
# coverage above is reported again with them removed.
lse <- (d$lsv_hi - d$lsv_lo) / (2 * stats::qnorm(0.975))
coll <- lse > 5 * stats::median(lse)
cat(sprintf("  log sv standard error: median %.4f, max %.4f\n",
            stats::median(lse), max(lse)))
cat(sprintf("  collapsed fits (se > 5x median): %d of %d, sv = %s\n",
            sum(coll), nrow(d),
            paste(formatC(exp(d$lsv[coll]), digits = 4,
                          format = "g"), collapse = " ")))
if (any(coll)) {
  e <- d[!coll, , drop = FALSE]
  cat("  with the collapsed fits removed:\n")
  row("  log sv", e$lsv, e$lsv_lo, e$lsv_hi, log(e$sv_true[1L]))
  row("  mu_cond", e$mu_cond, e$mu_cond_lo, e$mu_cond_hi, 0.9)
}

cat(" -- the trade-off, over replicates --\n")
# Spearman as well as Pearson, because one collapsed `sv` is an outlier
# on the log scale large enough to set the Pearson value on its own.
cat(sprintf("  cor(sv hat, mu_cond hat)  %7.4f pearson  %7.4f spearman\n",
            stats::cor(d$lsv, d$mu_cond),
            stats::cor(d$lsv, d$mu_cond, method = "spearman")))
cat(sprintf("  cor(sv hat, mu0 hat)      %7.4f pearson  %7.4f spearman\n",
            stats::cor(d$lsv, d$mu0),
            stats::cor(d$lsv, d$mu0, method = "spearman")))
cat(sprintf("  cor(sv hat, ndt hat)      %7.4f pearson  %7.4f spearman\n",
            stats::cor(d$lsv, d$ndt_hat),
            stats::cor(d$lsv, d$ndt_hat, method = "spearman")))
miss <- d$mu_cond_lo > 0.9 | d$mu_cond_hi < 0.9
if (any(miss)) {
  cat(sprintf("  replicates that miss on mu_cond: %s\n",
              paste(d$seed[miss], collapse = " ")))
  cat(sprintf("  their sv: %s\n",
              paste(formatC(exp(d$lsv[miss]), digits = 3,
                            format = "f"), collapse = " ")))
  cat(sprintf("  mean sv on the misses %.4f against %.4f on the rest\n",
              mean(exp(d$lsv[miss])), mean(exp(d$lsv[!miss]))))
}
