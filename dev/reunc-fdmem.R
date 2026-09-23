# Lane wt-reunc, punch round 1 (M1): what the ordinal fitted() route
# ALLOCATES, measured rather than asserted.
#
# The review named two dense objects: the Jacobian, at
# (n * K) x (p + levels), and the joint covariance, at
# (p + levels)^2. Batching removed the first. The second is the
# quadratic form's own matrix, because the variance of a row needs
# every pairwise covariance of the levels that row's block perturbs,
# and 0.61.0 already builds it for fitted()'s SCALAR route on any mixed
# fit. So the claim to check is which of the two this lane added.
#
# gc(reset = TRUE) then gc() gives R's own peak for the interval, in
# megabytes, which is the instrument that counts allocation rather than
# a clock.
#
#   Rscript dev/reunc-fdmem.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc2-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
lane <- !grepl("rellib-r3", find.package("frmtmb"), fixed = TRUE)

peak_mb <- function(fun) {
  invisible(gc(reset = TRUE, full = TRUE))
  invisible(fun())
  g <- gc(full = TRUE)
  # gc() names its megabyte columns "(Mb)" three times over, so the
  # max-used one is taken by position (6) rather than by name
  sum(g[, 6L])
}

make_fit <- function(ng, m = 6, seed = 5) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(ng), each = m)),
                  x = stats::rnorm(ng * m))
  lat <- 0.8 * d$x + stats::rnorm(ng, 0, 1)[d$g] + stats::rlogis(ng * m)
  d$y <- factor(cut(lat, c(-Inf, -1, 0.5, Inf), labels = FALSE),
                ordered = TRUE)
  list(fit = suppressWarnings(frm(bf(y ~ x + (1 | g)) + cumulative(),
                                  data = d)), d = d)
}

cat("\npeak Mb over the call, R's own gc() high-water mark\n")
cat(sprintf("%6s %10s %10s %10s %10s\n", "levels", "in", "in_full",
            "nd", "scalar"))
for (ng in c(200, 1000)) {
  z <- make_fit(ng)
  fit <- z$fit
  nd <- z$d[c(1, 7, 13, 19, 25), c("x", "g")]
  # the SCALAR route on a gaussian fit of the same shape, which 0.61.0
  # already runs through get_joint_cov(): the control for "this lane
  # added the dense covariance"
  dg <- z$d
  dg$yc <- as.numeric(dg$y) + stats::rnorm(nrow(dg))
  fg <- frm(bf(yc ~ x + (1 | g)), data = dg)
  m_in <- peak_mb(function() fitted(fit))
  m_nd <- peak_mb(function() fitted(fit, newdata = nd))
  m_sc <- peak_mb(function() fitted(fg))
  m_full <- if (lane) {
    gov <- frmtmb:::re_governed_b(fit)
    peak_mb(function() {
      frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x),
                         b_idx = gov)
    })
  } else NA_real_
  cat(sprintf("%6d %10.1f %10.1f %10.1f %10.1f\n", ng, m_in, m_full,
              m_nd, m_sc))
}
