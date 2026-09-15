# REVIEW of lane eamhier, ranked items 2 and 6.
#
# 6: where do "power 0.178" and "435 replicates" come from, and are they
#    the right instrument for a coverage count against a KNOWN nominal
#    rate?
# 2: is 0.694 at 14 replicates consistent with 0.913 at 60, or did
#    something other than the count change? The reference distribution
#    is se/sd over random 14-subsets of the same 60 records, which is
#    what "a ratio of two spreads estimated from 14 replicates is itself
#    noisy" has to mean if it means anything.
#
# Run: Rscript --vanilla dev/eamhierrev-power.R

rev_read <- function(dir) {
  fs <- sort(list.files(dir, pattern = "[.]tsv$", full.names = TRUE))
  rows <- lapply(fs, function(f) {
    kv <- strsplit(readLines(f, warn = FALSE)[[1L]], "\t",
                   fixed = TRUE)[[1L]]
    out <- as.list(trimws(sub("^[^=]*=", "", kv)))
    names(out) <- sub("=.*$", "", kv)
    out
  })
  nms <- unique(unlist(lapply(rows, names)))
  d <- do.call(rbind, lapply(rows, function(r) {
    r <- r[nms]; names(r) <- nms
    as.data.frame(lapply(r, function(x) if (is.null(x)) NA else x),
                  stringsAsFactors = FALSE)
  }))
  for (nm in nms) d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
  d
}

cat("== 6. where 0.178 and 435 come from ==\n")
za <- stats::qnorm(0.975); zb <- stats::qnorm(0.80)
p0 <- 0.95; p1 <- 0.90; n <- 60
one_pow <- stats::pnorm((abs(p1 - p0) - za * sqrt(p0 * (1 - p0) / n)) /
                          sqrt(p1 * (1 - p1) / n))
two_pow <- stats::pnorm((abs(p1 - p0) -
                           za * sqrt((p0 * (1 - p0) + p1 * (1 - p1)) /
                                       n)) /
                          sqrt((p0 * (1 - p0) + p1 * (1 - p1)) / n))
cat(sprintf("  ONE-sample z test against a KNOWN 0.95, n=60 : power %.4f\n",
            one_pow))
cat(sprintf(paste0("  TWO-sample z test, two proportions, n=60 ",
                   "each: power %.4f  <- the lane's 0.178\n"),
            two_pow))
n1 <- (za * sqrt(p0 * (1 - p0)) + zb * sqrt(p1 * (1 - p1)))^2 /
  (p1 - p0)^2
n2 <- (za + zb)^2 * (p0 * (1 - p0) + p1 * (1 - p1)) / (p1 - p0)^2
cat(sprintf("  n for 80%% power, ONE sample: %.1f\n", n1))
cat(sprintf("  n for 80%% power, TWO sample: %.1f  <- the lane's 435\n",
            n2))
wil <- function(k, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2); p <- k / n
  c(((p + z^2 / (2 * n)) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n),
    ((p + z^2 / (2 * n)) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n))
}
cat("\n  the rule this study ACTUALLY uses (Wilson excludes 0.95),",
    "exact binomial power:\n")
for (p in c(0.95, 0.93, 0.92, 0.90, 0.88, 0.85, 0.80)) {
  kmax <- max((0:n)[vapply(0:n, function(k) wil(k, n)[2L] < 0.95,
                           logical(1))])
  cat(sprintf("    truth %.2f: reject iff X <= %d; power %.4f\n", p,
              kmax, stats::pbinom(kmax, n, p)))
}
for (N in c(60, 100, 150, 180, 200, 205, 210, 250, 435)) {
  kmax <- max((0:N)[vapply(0:N, function(k) wil(k, N)[2L] < 0.95,
                           logical(1))])
  cat(sprintf(paste0("    n = %3d: reject iff X <= %3d; power ",
                     "at 0.90 %.4f; size at 0.95 %.4f\n"),
              N, kmax, stats::pbinom(kmax, N, 0.90),
              stats::pbinom(kmax, N, 0.95)))
}

cat("\n== 2. se/sd at 14 against its own sampling distribution ==\n")
C <- rev_read("dev/eamhier-rec/C")
sesd <- function(est, lo, hi, idx) {
  se <- (hi[idx] - lo[idx]) / (2 * stats::qnorm(0.975))
  mean(se) / stats::sd(est[idx])
}
set.seed(11)
for (v in c("lsv", "mu_cond")) {
  est <- C[[v]]; lo <- C[[paste0(v, "_lo")]]; hi <- C[[paste0(v, "_hi")]]
  full <- sesd(est, lo, hi, seq_len(60))
  first14 <- sesd(est, lo, hi, order(C$seed)[1:14])
  b <- replicate(20000, sesd(est, lo, hi, sample.int(60, 14)))
  cat(sprintf(paste0("  %-8s all 60 %.3f | first 14 %.3f | 14-subsets: ",
                     "mean %.3f sd %.3f, 2.5%% %.3f, 97.5%% %.3f\n"),
              v, full, first14, mean(b), stats::sd(b),
              stats::quantile(b, 0.025), stats::quantile(b, 0.975)))
  cat(sprintf(paste0("           P(a random 14-subset gives ",
                     "<= the first 14's value) = %.3f\n"),
              mean(b <= first14)))
}

cat("\n== 6b. is sd(log bs|s) at 5.1%% low a defect or ML shrinkage? ==\n")
A <- rev_read("dev/eamhier-rec/A")
# The ordinary ML downward bias of a variance component with q groups
# and p between-group degrees of freedom spent on fixed effects, times
# the Jensen term for reporting a standard deviation rather than a
# variance. Both random effects here sit on an intercept only: the
# condition effect is WITHIN subject and spends no between-subject df.
q <- 30
shrink <- function(p) sqrt(1 - p / q) * (1 - 1 / (4 * (q - p)))
for (p in 1:2) {
  cat(sprintf("  expected sd_hat / sd_true at q=%d, p=%d: %.5f (%.2f%% low)\n",
              q, p, shrink(p), 100 * (shrink(p) - 1)))
}
for (v in c("sd1", "sd2")) {
  truth <- if (v == "sd1") 0.35 else 0.20
  nm <- if (v == "sd1") "sd(mu|s)" else "sd(log bs|s)"
  m <- mean(A[[v]]); mc <- stats::sd(A[[v]]) / sqrt(60)
  cat(sprintf(paste0("  %-14s mean %.5f mcse %.5f | z vs truth ",
                     "%6.2f | z vs ML expectation (p=1) %6.2f\n"),
              nm, m, mc, (m - truth) / mc,
              (m - truth * shrink(1)) / mc))
}
cat("\ndone\n")
