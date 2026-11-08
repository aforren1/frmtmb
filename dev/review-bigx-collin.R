# Near-aliased designs that the 1e-7 tolerance dropped and 1e-11 keeps.
# Against glm() (tol 1e-11) and lm() (1e-7). REV_ARM selects the build.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")

report <- function(lbl, fo, d, fam = gaussian(), gfam = stats::gaussian()) {
  r <- rfit(bf(fo), family = fam, data = d)
  g <- suppressWarnings(glm(fo, family = gfam, data = d))
  lmr <- qr(model.matrix(fo, d))$rank
  gse <- suppressWarnings(summary(g)$coefficients[, 2])
  if (r$err) {
    cat(sprintf("%-40s ERROR %s\n", lbl, conditionMessage(r$f)))
    return(invisible())
  }
  f <- r$f
  fe <- tryCatch(quiet(fixef(f)), error = function(e) NULL)
  se <- if (is.null(fe)) NA else fe[, 2]
  cat(sprintf(paste0("%-40s kept=%d (glm %d, lm-tol %d) conv=%s basis=%s ",
                     "ll=%.6f glm ll=%.6f dll=%+.2e\n",
                     "    max SE frm=%.3g glm=%.3g; SE ratio frm/glm ",
                     "range %.4g..%.4g; any NaN SE %s\n",
                     "    warn: %s\n    msg: %s\n"),
              lbl, length(beta_of(f)), g$rank, lmr, f$opt$convergence,
              basis_on(f), ll_of(f), as.numeric(logLik(g)),
              ll_of(f) - as.numeric(logLik(g)),
              max(se, na.rm = TRUE), max(gse),
              if (length(se) == length(gse)) min(se / gse) else NA,
              if (length(se) == length(gse)) max(se / gse) else NA,
              anyNA(se), short_w(r),
              if (length(r$msg)) substr(paste(unique(r$msg), collapse = " | "), 1, 220) else "none"))
}

for (s in 1:3) {
  set.seed(500 + s)
  n <- 200
  x1 <- rnorm(n)
  d <- data.frame(x1 = x1)
  d$x2a <- x1 + 1e-9 * rnorm(n)
  d$x2b <- x1 + 1e-8 * rnorm(n)
  d$x2c <- x1 + 1e-6 * rnorm(n)
  d$y <- 1 + 0.5 * x1 + rnorm(n)
  d$yp <- rpois(n, exp(0.3 + 0.4 * x1))
  cat("== seed", 500 + s, "\n")
  report("x2 = x1 + 1e-9 noise, gaussian", y ~ x1 + x2a, d)
  report("x2 = x1 + 1e-8 noise, gaussian", y ~ x1 + x2b, d)
  report("x2 = x1 + 1e-6 noise, gaussian", y ~ x1 + x2c, d)
  report("x2 = x1 + 1e-9 noise, poisson", yp ~ x1 + x2a, d, poisson(),
         stats::poisson())
  # proportions rounded to 9 digits: a dummy trap broken only by rounding
  p <- matrix(runif(n * 3), n); p <- p / rowSums(p)
  d$p1 <- round(p[, 1], 9); d$p2 <- round(p[, 2], 9); d$p3 <- round(p[, 3], 9)
  d$yq <- 1 + 2 * d$p1 - d$p2 + rnorm(n)
  report("p1+p2+p3 = 1 to 1e-9, gaussian", yq ~ p1 + p2 + p3, d)
  d$p1r <- round(p[, 1], 12); d$p2r <- round(p[, 2], 12)
  d$p3r <- round(p[, 3], 12)
  report("p1+p2+p3 = 1 to 1e-12, gaussian", yq ~ p1r + p2r + p3r, d)
  # raw polynomials at large x
  xx <- runif(n, 1000, 1010)
  d$xx <- xx
  d$yy <- sin((xx - 1000) / 2) + rnorm(n, 0, 0.3)
  for (deg in 5:8) {
    report(sprintf("poly(xx in [1000,1010], %d, raw)", deg),
           as.formula(sprintf("yy ~ poly(xx, %d, raw = TRUE)", deg)), d)
  }
  d$year <- sample(1990:2020, n, TRUE)
  d$yr <- 0.01 * (d$year - 2005) + 0.002 * (d$year - 2005)^2 + rnorm(n)
  report("year + I(year^2)", yr ~ year + I(year^2), d)
  d$ypy <- rpois(n, exp(0.2 + 0.03 * (d$year - 2005)))
  report("poisson year + I(year^2)", ypy ~ year + I(year^2), d, poisson(),
         stats::poisson())
}
