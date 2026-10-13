# Reviewer of lane seflags, re-check of B1: se_design_null() against
# aliases that are exact only to rounding (centered, scaled and shifted,
# poly() against a smooth's fixed column), without random effects, a
# rank deficiency the frame already handles (double counting?), a near
# alias the QR tolerance may call exact, an alias across two nonlinear
# parameters (not a design alias), and a tensor smooth's null space.
# Per case: code, lost set, design-null columns, fitted(newdata) SEs and
# warnings against the same likelihood without the alias, and a
# hypothesis on one aliased coefficient.
#   Rscript dev/seflags-rev3-alias.R <lane|base|op2>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
mk <- function(seed = 1, m = 3, n = 200) {
  set.seed(seed)
  x <- stats::rnorm(n, m, 1)
  z <- stats::rnorm(n)
  d <- data.frame(x = x, z = z, g = factor(rep(1:10, length.out = n)))
  d$y <- 1 + 0.5 * x + sin(2 * x) + 0.3 * z + 0.2 * x * z +
    stats::rnorm(10, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.3)
  d
}
row <- function(lab, fo, ref, d, hyp = NULL, nd = d[c(1, 50, 100, 150), ]) {
  r <- cap(frm(fo, data = d))
  f <- r$v
  if (inherits(f, "rev_err")) {
    cat(sprintf("%-28s ERROR %s\n", lab, substr(f, 1, 150)))
    return(invisible())
  }
  nn <- tryCatch(NCOL(ns$se_design_null(f)), error = function(e) NA)
  ref_fit <- suppressWarnings(frm(ref, data = d))
  re_na <- if (length(f$frame$re_blocks)) NA else NULL
  p <- cap(fitted(f, newdata = nd, re_formula = re_na))
  p0 <- suppressWarnings(fitted(ref_fit, newdata = nd, re_formula = re_na))
  rel <- if (is.matrix(p$v)) max(abs(p$v[, 2] / p0[, 2] - 1)) else NA
  h <- if (!is.null(hyp)) cap(hypothesis(f, hyp, class = NULL)) else NULL
  hs <- if (is.null(h)) "-" else if (inherits(h$v, "rev_err")) "ERROR" else
    sprintf("%s (%d warn)", format(h$v$hypothesis$Est.Error, digits = 4),
            length(h$w))
  cat(sprintf(paste0("%-28s code %d | lost %s | design-null cols %s |",
                     " fit warn %d\n     fitted se %s warn %d | max",
                     " |se/ref-1| %.2g | hyp %s: %s\n"),
              lab, f$opt$convergence, lost_of(f), nn, length(r$w),
              if (is.matrix(p$v)) paste(signif(p$v[, 2], 4), collapse = " ")
              else substr(as.character(p$v), 1, 80),
              length(p$w), rel, if (is.null(hyp)) "-" else hyp, hs))
  for (w in c(r$w, p$w)) cat("     W:", substr(w, 1, 160), "\n")
}
for (seed in 1:2) {
  d <- mk(seed)
  cat("\n== seed", seed, "\n")
  d$xc <- d$x - mean(d$x)
  row("centered xc + s(x)", y ~ xc + s(x) + (1 | g), y ~ s(x) + (1 | g), d,
      "b_xc = 0")
  d$xk <- 1e4 * d$x + 5e4
  row("scaled 1e4 x + 5e4 + s(x)", y ~ xk + s(x) + (1 | g),
      y ~ s(x) + (1 | g), d, "b_xk = 0")
  row("poly(x, 2) + s(x)", y ~ poly(x, 2) + s(x) + (1 | g),
      y ~ s(x) + I(x^2) + (1 | g), d)
  row("x + s(x), no RE", y ~ x + s(x), y ~ s(x), d, "b_x = 0")
  d$x3 <- d$x + d$z
  row("frame drop x + z + x3", y ~ x + z + x3 + (1 | g),
      y ~ x + z + (1 | g), d, "b_x = 0")
  d$x2 <- 2 * d$x
  row("frame drop x + 2x + s(x)", y ~ x + x2 + s(x) + (1 | g),
      y ~ s(x) + (1 | g), d, "b_x = 0")
  for (eps in c(1e-4, 1e-6, 1e-8)) {
    set.seed(100 + seed)
    d$xe <- d$x + stats::rnorm(nrow(d), 0, eps)
    row(sprintf("near alias xe(eps %g) + s(x)", eps),
        y ~ xe + s(x) + (1 | g), y ~ s(x) + (1 | g), d, "b_xe = 0")
  }
  row("nl a ~ 1 + x, b ~ 1 + x",
      bf(y ~ a + b, a ~ 1 + x + (1 | g), b ~ 1 + x, nl = TRUE),
      y ~ x + (1 | g), d, "b_a_x = 0")
  row("x + z + t2(x, z)", y ~ x + z + t2(x, z) + (1 | g),
      y ~ t2(x, z) + (1 | g), d, "b_x = 0")
}
