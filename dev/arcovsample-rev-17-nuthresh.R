# REVIEW re-check, detail: the exact nu at which the as-written
# Student-t head stops being a number, against the 3.6e305 the rewritten
# core NEWS entry names. Bisected in LOG space: a geometric mean formed
# as sqrt(lo * hi) overflows to Inf at these magnitudes and collapses
# the search, which it did on the first try.
#
#   Rscript dev/arcovsample-rev-17-nuthresh.R

old_head <- function(nu, K) lgamma((nu + K) / 2) - lgamma(nu / 2)
bisect <- function(fn, llo, lhi) {
  for (it in 1:200) {
    lm <- (llo + lhi) / 2
    if (fn(exp(lm))) llo <- lm else lhi <- lm
  }
  exp(llo)
}
for (K in c(1L, 2L, 3L, 5L)) {
  th <- bisect(function(nu) is.finite(old_head(nu, K)),
               log(1e300), log(1e308))
  cat("K = ", K, ": as-written head finite to nu = ",
      format(th, digits = 8), ", NaN above it\n", sep = "")
}
cat("lgamma() finite to argument ",
    format(bisect(function(a) is.finite(lgamma(a)), log(1e300),
                  log(1e308)), digits = 8), "\n", sep = "")
cat("\nthe head at nu values the NEWS entry brackets, K = 2:\n")
for (nu in c(1e305, 2e305, 3e305, 3.6e305, 4e305, 5e305, 5.1e305,
             5.2e305)) {
  cat(sprintf("  nu %9.2e  as-written %14s   truth about %10.4f\n", nu,
              format(old_head(nu, 2L)), log(nu / 2)))
}
cat("DONE\n")
