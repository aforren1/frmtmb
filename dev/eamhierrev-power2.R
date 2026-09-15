# REVIEW re-check: is the power of "the Wilson interval excludes 0.95"
# monotone in n? The lane says no, and reports n = 180 where power
# first reaches 0.80 and n = 202 from where it stays there. This
# review had given 205 as a single threshold.
#
# The rejection region moves in WHOLE COUNTS, so the test's size drops
# in a sawtooth as n grows and power rides the same sawtooth. Scanned
# rather than argued.
#
# Run: Rscript --vanilla dev/eamhierrev-power2.R

wil <- function(k, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  c(((p + z^2 / (2 * n)) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n),
    ((p + z^2 / (2 * n)) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n))
}

# the largest count whose Wilson interval still lies below 0.95
kmax <- function(n) {
  k <- 0:n
  hit <- vapply(k, function(kk) wil(kk, n)[2L] < 0.95, logical(1))
  if (!any(hit)) return(NA_integer_)
  max(k[hit])
}

ns <- 20:320
km <- vapply(ns, kmax, integer(1))
pw <- stats::pbinom(km, ns, 0.90)
sz <- stats::pbinom(km, ns, 0.95)

cat("== power of the Wilson-excludes-0.95 rule, truth 0.90 ==\n")
first <- ns[which(pw >= 0.80)[1L]]
runs <- rle(pw >= 0.80)
# the start of the last TRUE run that never ends inside the scan
last_false <- max(which(!(pw >= 0.80)))
sustained <- ns[last_false + 1L]
cat(sprintf("  first n with power >= 0.80 : %d (power %.4f, size %.4f)\n",
            first, pw[ns == first], sz[ns == first]))
cat(sprintf("  power dips below again at  : %d (power %.4f)\n",
            ns[last_false], pw[ns == ns[last_false]]))
cat(sprintf("  sustained from n           : %d (power %.4f, size %.4f)\n",
            sustained, pw[ns == sustained], sz[ns == sustained]))
cat(sprintf("  monotone in n?             : %s\n",
            if (all(diff(pw) >= 0)) "YES" else "NO"))
cat(sprintf("  largest drop between n and n+1: %.4f at n = %d\n",
            -min(diff(pw)), ns[which.min(diff(pw))]))
cat("\n  the sawtooth, n = 170 to 215:\n")
for (i in which(ns >= 170 & ns <= 215)) {
  cat(sprintf("   n=%3d  X<=%3d  power %.4f  size %.4f%s\n", ns[i],
              km[i], pw[i], sz[i],
              if (pw[i] >= 0.80) "  >= 0.80" else ""))
}
cat("\n  every n in 20:320 with power >= 0.80 that is followed by one\n")
cat("  below it:\n")
j <- which(pw >= 0.80 & c(pw[-1L] < 0.80, FALSE))
cat(sprintf("   n = %s\n", paste(ns[j], collapse = ", ")))

cat("\n== the figures the revised sentence rests on, n = 60 ==\n")
k60 <- kmax(60)
cat(sprintf("  rejection region X <= %d of 60 (%.1f percent)\n", k60,
            100 * k60 / 60))
cat(sprintf("  size at truth 0.95 : %.4f\n", stats::pbinom(k60, 60, 0.95)))
for (p in c(0.90, 0.88, 0.85)) {
  cat(sprintf("  power at truth %.2f: %.4f\n", p,
              stats::pbinom(k60, 60, p)))
}
cat(sprintf("  power AT the boundary rate %.4f: %.4f\n", k60 / 60,
            stats::pbinom(k60, 60, k60 / 60)))
cat("\ndone\n")
