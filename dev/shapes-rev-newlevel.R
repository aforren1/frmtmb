# Reviewer, priority 2: predict(allow_new_levels = TRUE) claims to give
# a PREDICTIVE interval for a row in a group the fit never saw. brms
# draws a fresh group effect there, so its interval carries tau^2. This
# measures whether frmtmb's does: fit on ng groups, predict at m rows
# in BRAND NEW groups drawn from the same truth, and count coverage.
#
# The comparison arm is the same count at rows in KNOWN groups, which
# the lane already measured near nominal, so a gap between the two arms
# is the missing tau^2 and not a property of the harness.
#
#   Rscript dev/shapes-rev-newlevel.R <nrep>

.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
a <- commandArgs(trailingOnly = TRUE)
nrep <- if (length(a)) as.integer(a[1]) else 220L
n <- 60L; m <- 10L; ndraws <- 600L; ng <- 12L
TAU <- 0.7; SIG <- 0.8
SEED0 <- 883000L

hits <- list(new = integer(0), known = integer(0))
width <- list(new = numeric(0), known = numeric(0))
tot <- 0L; fail <- 0L
for (r in seq_len(nrep)) {
  s <- SEED0 + 100L * r
  set.seed(s)
  d <- data.frame(x = rnorm(n),
                  g = factor(rep(seq_len(ng), length.out = n)))
  u <- rnorm(ng, 0, TAU)
  d$y <- rnorm(n, -0.4 + 1.3 * d$x + u[d$g], SIG)
  fit <- tryCatch(suppressWarnings(
    frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)),
    error = function(e) NULL)
  if (is.null(fit)) { fail <- fail + 1L; next }

  # BRAND NEW groups: one fresh effect each, from the same truth
  gnew <- paste0("new", seq_len(m))
  unew <- rnorm(m, 0, TAU)
  ndn <- data.frame(x = rnorm(m),
                    g = factor(gnew, levels = c(levels(d$g), gnew)))
  yn <- rnorm(m, -0.4 + 1.3 * ndn$x + unew, SIG)
  # KNOWN groups, the control arm
  gk <- sample(seq_len(ng), m, TRUE)
  ndk <- data.frame(x = rnorm(m),
                    g = factor(gk, levels = levels(d$g)))
  yk <- rnorm(m, -0.4 + 1.3 * ndk$x + u[gk], SIG)

  set.seed(s + 1L)
  pn <- tryCatch(suppressWarnings(
    predict(fit, newdata = ndn, ndraws = ndraws,
            allow_new_levels = TRUE)), error = function(e) NULL)
  set.seed(s + 2L)
  pk <- tryCatch(suppressWarnings(
    predict(fit, newdata = ndk, ndraws = ndraws)), error = function(e) NULL)
  if (is.null(pn) || is.null(pk)) { fail <- fail + 1L; next }
  hits$new <- c(hits$new, sum(yn >= pn[, 3L] & yn <= pn[, 4L]))
  hits$known <- c(hits$known, sum(yk >= pk[, 3L] & yk <= pk[, 4L]))
  width$new <- c(width$new, mean(pn[, 4L] - pn[, 3L]))
  width$known <- c(width$known, mean(pk[, 4L] - pk[, 3L]))
  tot <- tot + m
}
cat("replicates:", length(hits$new), " failed:", fail, "\n")
cat("n:", n, " groups:", ng, " new points:", m, " ndraws:", ndraws,
    " tau:", TAU, " sigma:", SIG, " SEED0:", SEED0, "\n")
cat("script: dev/shapes-rev-newlevel.R\n\n")
for (nm in names(hits)) {
  v <- hits[[nm]]; k <- sum(v)
  ci <- stats::binom.test(k, tot, 0.95)$conf.int
  pr <- v / m; se <- stats::sd(pr) / sqrt(length(pr))
  cat(sprintf("%-6s %5d of %5d = %.4f  (%.4f, %.4f) binomial; ",
              nm, k, tot, k / tot, ci[1], ci[2]))
  cat(sprintf("replicate mean %.4f se %.4f  z=%.2f  mean width %.4f\n",
              mean(pr), se, (mean(pr) - 0.95) / se, mean(width[[nm]])))
}
cat(sprintf("\nwidth ratio new/known: %.4f\n",
            mean(width$new) / mean(width$known)))
cat(sprintf("what brms's fresh-effect draw would need: %.4f\n",
            sqrt(SIG^2 + TAU^2) / SIG))
