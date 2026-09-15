# Lane eamhier, punch round 1: re-derive every number the review
# corrected, from this lane's own records, before any of them is
# written into a document again.
#
# Run: Rscript --vanilla dev/eamhier-scripts/eamhier-punch1.R
source("dev/eamhier-scripts/eamhier-summary-read.R")

cat("== B. power, one sample against a KNOWN 0.95 ==\n")
# The rule this study applies: reject when the Wilson interval on X of
# n excludes 0.95. So the rejection region is found from the rule, not
# assumed, and the power is the exact binomial probability of landing
# in it.
wilson <- function(k, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  c(((p + z^2 / (2 * n)) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n),
    ((p + z^2 / (2 * n)) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n))
}
reject_max <- function(n) {
  k <- 0:n
  hi <- vapply(k, function(i) wilson(i, n)[2L], numeric(1))
  max(k[hi < 0.95])
}
pow <- function(n, p) stats::pbinom(reject_max(n), n, p)
cat("  n = 60, rejection region X <=", reject_max(60), "\n")
for (p in c(0.95, 0.93, 0.92, 0.90, 0.88, 0.85, 0.80)) {
  cat(sprintf("    truth %.2f  power %.4f\n", p, pow(60, p)))
}
ns <- 60:400
ok <- vapply(ns, pow, numeric(1), p = 0.90) >= 0.80
# Power is NOT monotone in n here: the rejection region moves in whole
# counts, so it saws. Both numbers are reported, because quoting only
# the first crossing would promise 80 percent at a count that does not
# keep it.
cat("  first n whose power reaches 0.80 at a true 0.90: ",
    min(ns[ok]), "\n")
cat("  smallest n from which it STAYS at or above 0.80: ",
    ns[min(which(rev(cumprod(rev(ok))) == 1L))], "\n")
# the two-sample figures the findings used, reproduced so that the
# correction names what was wrong rather than only what is right
cat("  two-sample formula, for comparison: power ",
    formatC(stats::power.prop.test(n = 60, p1 = 0.95, p2 = 0.90,
                                   sig.level = 0.05)$power, digits = 4,
            format = "f"),
    ", n ",
    formatC(stats::power.prop.test(p1 = 0.95, p2 = 0.90,
                                   sig.level = 0.05,
                                   power = 0.8)$n, digits = 1,
            format = "f"), "\n", sep = "")

cat("\n== C. the sv range in the single-level control ==\n")
s1 <- read_records("dev/eamhier-rec/S1")
s1 <- s1[s1$status == "ok", , drop = FALSE]
lse <- (s1$lsv_hi - s1$lsv_lo) / (2 * stats::qnorm(0.975))
coll <- lse > 5 * stats::median(lse)
sv <- exp(s1$lsv)
cat("  all 60, six smallest: ",
    paste(formatC(sort(sv)[1:6], digits = 6, format = "f"),
          collapse = " "), "\n")
cat("  collapsed:", sum(coll), " usable range: ",
    formatC(min(sv[!coll]), digits = 6, format = "f"), " to ",
    formatC(max(sv[!coll]), digits = 6, format = "f"), "\n", sep = "")
cat("  seed of the smallest usable: ", s1$seed[!coll][
  which.min(sv[!coll])], "\n")

cat("\n== F. the within-fit signature of the collapse ==\n")
one <- s1[coll, , drop = FALSE]
se_of <- function(lo, hi) (hi - lo) / (2 * stats::qnorm(0.975))
ses <- c(mu0 = se_of(one$mu0_lo, one$mu0_hi),
         mu_cond = se_of(one$mu_cond_lo, one$mu_cond_hi),
         lbs = se_of(one$lbs_lo, one$lbs_hi),
         lsv = se_of(one$lsv_lo, one$lsv_hi))
cat("  seed", one$seed, " se(log sv)",
    formatC(ses[["lsv"]], digits = 4, format = "f"), "\n")
cat("  the others: ",
    paste(names(ses)[1:3], formatC(ses[1:3], digits = 6, format = "f"),
          collapse = "  "), "\n")
# The record does not carry the ndt coefficient's LINK-scale interval,
# and that one is the largest of the four, so it is read back from the
# refit that eamhier-sv-collapse.R prints.
ndt_se_link <- (2.4503451 - 2.2539783) / (2 * stats::qnorm(0.975))
cat("  ndt (link scale), from the refit's confint: ",
    formatC(ndt_se_link, digits = 6, format = "f"), "\n")
cat("  ratio to the largest of the other four: ",
    formatC(ses[["lsv"]] / ndt_se_link, digits = 0, format = "f"),
    "\n")
cat("  1 of 60, Wilson: ",
    paste(formatC(100 * wilson(1, 60), digits = 1, format = "f"),
          collapse = " to "), " percent\n")

cat("\n== E. arm B's ratio is an identity: the residual ==\n")
b <- read_records("dev/eamhier-rec/B")
b <- b[b$status == "ok" & b$bound == "pg", , drop = FALSE]
ratio <- b$ndt_sd_hat_ms / (b$ndt_frac * b$floor_sd_ms)
cat("  mean ", format(mean(ratio), digits = 10),
    "  sd ", format(stats::sd(ratio), digits = 4),
    "\n  min ", format(min(ratio), digits = 10),
    "  max ", format(max(ratio), digits = 10), "\n", sep = "")
cat("  identical(ratio, 1) for how many: ",
    sum(vapply(ratio, function(r) isTRUE(all.equal(r, 1,
                                                   tolerance = 0)),
               logical(1))), " of ", length(ratio), "\n", sep = "")
cat("  cor(ratio - 1, fitted ndt component): ",
    formatC(stats::cor(ratio - 1, b$sd2), digits = 4, format = "f"),
    "\n")
cat("  fitted ndt component: mean ", format(mean(b$sd2), digits = 4),
    ", under 1e-4 on ", sum(b$sd2 < 1e-4), " of ", nrow(b), "\n",
    sep = "")

cat("\n== the retraction: was 0.694 ordinary noise? ==\n")
cc <- read_records("dev/eamhier-rec/C")
cc <- cc[cc$status == "ok", , drop = FALSE]
cc <- cc[order(cc$seed), , drop = FALSE]
sesd <- function(x, est, lo, hi, truth) {
  se <- (hi - lo) / (2 * stats::qnorm(0.975))
  mean(se) / stats::sd(est)
}
first14 <- cc[1:14, , drop = FALSE]
cat("  first 14 seeds: log sv se/sd ",
    formatC(mean((first14$lsv_hi - first14$lsv_lo) /
                   (2 * stats::qnorm(0.975))) / stats::sd(first14$lsv),
            digits = 4, format = "f"),
    ", mu_cond se/sd ",
    formatC(mean((first14$mu_cond_hi - first14$mu_cond_lo) /
                   (2 * stats::qnorm(0.975))) /
              stats::sd(first14$mu_cond), digits = 4, format = "f"),
    "\n", sep = "")
set.seed(1L)
draws <- replicate(20000L, {
  j <- sample(nrow(cc), 14L)
  x <- cc[j, , drop = FALSE]
  mean((x$lsv_hi - x$lsv_lo) / (2 * stats::qnorm(0.975))) /
    stats::sd(x$lsv)
})
cat("  20,000 random 14-subsets: mean ",
    formatC(mean(draws), digits = 3, format = "f"), ", sd ",
    formatC(stats::sd(draws), digits = 3, format = "f"),
    ", 2.5 to 97.5 percent ",
    paste(formatC(stats::quantile(draws, c(0.025, 0.975)), digits = 3,
                  format = "f"), collapse = " to "), "\n", sep = "")
cat("  P(a random 14-subset is at or below 0.694) = ",
    formatC(mean(draws <= 0.694), digits = 4, format = "f"), "\n")

cat("\n== the variance components against ML shrinkage ==\n")
a <- read_records("dev/eamhier-rec/A")
a <- a[a$status == "ok", , drop = FALSE]
q <- 30
shrink <- sqrt(1 - 1 / q) * (1 - 1 / (4 * q - 4))
cat("  sqrt(1 - 1/30) * (1 - 1/116) = ", format(shrink, digits = 6),
    ", that is ", formatC(100 * (1 - shrink), digits = 2,
                          format = "f"), " percent low\n", sep = "")
for (nm in c("sd1", "sd2")) {
  tr <- if (nm == "sd1") 0.35 else 0.20
  m <- mean(a[[nm]])
  se <- stats::sd(a[[nm]]) / sqrt(nrow(a))
  cat(sprintf("  %s mean %.5f mcse %.5f  z vs truth %6.2f  z vs ML %6.2f\n",
              nm, m, se, (m - tr) / se, (m - tr * shrink) / se))
}

cat("\n== the sawtooth is paid for in SIZE ==\n")
for (n in c(60, 180, 184, 185, 199, 200, 202)) {
  k <- reject_max(n)
  cat(sprintf("  n %3d  region X <= %3d  size %.4f  power at 0.90 %.4f\n",
              n, k, stats::pbinom(k, n, 0.95), stats::pbinom(k, n, 0.90)))
}
