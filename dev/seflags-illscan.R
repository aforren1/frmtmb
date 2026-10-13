# Lane seflags, punch round 1: degree-8 (and 9) raw polynomials with and
# without (1 | g) on [1, 2], seeds 1 to 10: the code, the lost set and
# the coefficients' SE against lmer() or lm(), to find a fit that
# reaches the "illcond" rule on this lane alone.
#   Rscript dev/seflags-illscan.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
for (deg in 8:9) for (re in c(TRUE, FALSE)) for (s in 1:10) {
  set.seed(s)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:deg, collapse = " + ")
  fo <- as.formula(paste("y ~", rhs, if (re) "+ (1 | g)"))
  f <- suppressWarnings(frm(fo, data = d))
  l <- suppressWarnings(ns$sdr_of(f))$se_lost
  cat(sprintf("deg %d re %s seed %2d code %d | lost %s\n", deg, re, s,
              f$opt$convergence,
              if (length(l)) paste(unique(l), table(l)[unique(l)], collapse = " ")
              else "-"))
}
