# Lane wt-resmooth. An INDEPENDENT check of the finite-difference
# standard error of an ordinal category probability on a fit with a
# smooth, because adding the smooth's coefficients to the differenced set
# moves it a long way and in both directions (dev/resmooth-fdse.txt,
# dev/resmooth-popse-before.txt). The reference is Monte Carlo: draw the
# outer parameters AND the smooth's coefficients jointly from the joint
# covariance the delta method uses, recompute the probabilities, take the
# standard deviation. That is the quantity the delta method approximates,
# so the two must agree up to curvature and Monte Carlo error.
#   Rscript dev/resmooth-mcse.R > dev/resmooth-mcse.txt
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", find.package("frmtmb"), "\n")

mc_check <- function(lab, fit, nd, ndraw = 4000, seed = 4) {
  bl <- fit$frame[["re_blocks"]]
  sm <- sort(unique(unlist(lapply(bl[vapply(bl, function(b) {
    b[["covstruct"]] %in% c("smooth", "gp", "hsgp")
  }, NA)], `[[`, "b_idx"))))
  f <- function(x) frmtmb:::fitted_point(x, nd, re_formula = NA)
  ds <- frmtmb:::fit_draw_space(fit)
  map <- ds$map
  V <- frmtmb:::fd_joint_cov(fit, map, sm)
  v0 <- frmtmb:::fit_outer_vector(fit, map)
  b0 <- fit$estimates[["b"]]
  p <- length(v0)
  L <- t(chol((V + t(V)) / 2 + diag(1e-10, nrow(V))))
  set.seed(seed)
  m0 <- as.matrix(f(fit))
  out <- matrix(NA_real_, ndraw, length(m0))
  for (i in seq_len(ndraw)) {
    z <- as.vector(L %*% stats::rnorm(nrow(V)))
    g <- frmtmb:::fit_set_outer(fit, v0 + z[seq_len(p)], map)
    g$estimates[["b"]][sm] <- b0[sm] + z[p + seq_along(sm)]
    g$cache <- new.env(parent = emptyenv())
    out[i, ] <- as.vector(as.matrix(f(g)))
  }
  mc <- apply(out, 2, stats::sd)
  with_b <- as.vector(frmtmb:::fit_fd_se(fit, f, b_idx = sm, b_batch = NULL))
  no_b <- as.vector(frmtmb:::fit_fd_se(fit, f, b_idx = NULL))
  cat("==", lab, "| smooth coefficients:", length(sm), "| draws:", ndraw,
      "\n")
  cat("  Monte Carlo sd      :", sprintf("%.5f", mc), "\n")
  cat("  delta, smooth b IN  :", sprintf("%.5f", with_b), "\n")
  cat("  delta, smooth b OUT :", sprintf("%.5f", no_b), "\n")
  cat(sprintf("  median |IN/MC - 1| %.4f | median |OUT/MC - 1| %.4f\n",
              stats::median(abs(with_b / mc - 1)),
              stats::median(abs(no_b / mc - 1))))
}

# a POPULATION smooth, where the shipped value was LARGER than the
# delta method with the smooth's coefficients in
set.seed(23)
n <- 240
d1 <- data.frame(x = stats::runif(n, -2, 2))
lat <- 1.2 * sin(2 * d1$x) + stats::rlogis(n)
d1$y <- factor(cut(lat, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
f1 <- suppressWarnings(frm(bf(y ~ s(x, k = 8)) + cumulative(), data = d1))
mc_check("cumulative, s(x, k = 8)", f1, data.frame(x = c(-1.5, 0, 1.5)))

# a FACTOR smooth, where it was smaller
set.seed(21)
ng <- 6
d2 <- data.frame(g = factor(rep(seq_len(ng), each = 25)),
                 x = stats::runif(ng * 25, -2, 2))
lat2 <- stats::rnorm(ng, 0, 1)[d2$g] * sin(d2$x) + 0.5 * d2$x +
  stats::rlogis(nrow(d2))
d2$y <- factor(cut(lat2, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
f2 <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5)) + cumulative(),
                           data = d2))
mc_check("cumulative, s(x, g, bs = fs)", f2, d2[c(3, 40, 90), c("x", "g")])
