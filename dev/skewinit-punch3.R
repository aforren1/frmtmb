# Punch round 3, two record fixes, both verified here rather than
# copied: does autoscale rescue the no-intercept poisson at small
# covariate scales, and does the blocker's mechanism survive inside
# the one dpar that opted in (alpha ~ 0 + xs)?
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")
op <- options(digits = 12)

cat("\n=== 1. does autoscale rescue the poisson no-intercept case?\n")
set.seed(23)
n <- 250
xt <- rnorm(n)
y <- rpois(n, pmin(exp(0.5 + 0.4 * xt), 1e6))
cat(sprintf("%-8s %16s %16s %16s %12s %12s\n", "scale", "glm",
            "autoscale off", "autoscale ON", "off - glm", "ON - glm"))
for (s in c(1, 1e-2, 1e-4, 1e-6, 1e-8)) {
  dd <- data.frame(y = y, xs = xt * s)
  gl <- as.numeric(stats::logLik(
    stats::glm(y ~ 0 + xs, family = stats::poisson(), data = dd)))
  f0 <- suppressWarnings(frm(bf(y ~ 0 + xs), family = poisson(), data = dd,
                             control = frmtmb_control(autoscale = FALSE)))
  f1 <- suppressWarnings(frm(bf(y ~ 0 + xs), family = poisson(), data = dd,
                             control = frmtmb_control(autoscale = TRUE)))
  cat(sprintf("%-8g %16.6f %16.6f %16.6f %12.6f %12.6f\n", s, gl,
              as.numeric(logLik(f0)), as.numeric(logLik(f1)),
              as.numeric(logLik(f0)) - gl, as.numeric(logLik(f1)) - gl))
}

cat("\n=== 2. the same mechanism INSIDE the declaring dpar\n")
mk <- function(seed, s, n = 250) {
  set.seed(seed)
  xs <- -abs(rnorm(n)) * 3
  yy <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  data.frame(y = yy, xs = xs * s, xr = xs)
}
fo <- bf(y ~ xr, sigma ~ 1, alpha ~ 0 + xs)
cat(sprintf("%-8s %16s %16s %16s %10s\n", "scale", "default", "best sweep",
            "start coef", "conv/msg"))
for (s in c(1, 1e-3, 1e-6, 1e-9)) {
  dd <- mk(1, s)
  u <- frm(fo, family = skew_normal(), data = dd, dry_run = "objective")
  b0 <- u$estimates$betad[["alpha_xs"]]
  f <- suppressWarnings(frm(fo, family = skew_normal(), data = dd))
  # a sweep on the PREDICTOR scale: alpha = v means coef = v / mean|xs|
  sw <- vapply(c(-4, -2, -1, 1, 2, 4), function(v) {
    st <- list(betad = c(log(stats::sd(dd$y)),
                         v / sqrt(mean(dd$xs^2))))
    as.numeric(logLik(suppressWarnings(
      frm(fo, family = skew_normal(), data = dd, start = st))))
  }, 0)
  cat(sprintf("%-8g %16.7f %16.7f %16.4g %s\n", s,
              as.numeric(logLik(f)), max(sw), b0,
              substr(f$opt$message, 1, 18)))
}

cat("\n=== 3. paired sweep, lane against BASE, 12 seeds x 4 scales\n")
cat("   (run this script on both builds and pair the CSVs)\n")
rows <- list()
for (s in c(1, 1e-3, 1e-6, 1e-9)) {
  for (sd in 1:12) {
    dd <- mk(sd, s)
    f <- suppressWarnings(frm(fo, family = skew_normal(), data = dd))
    rows[[length(rows) + 1L]] <- data.frame(
      scale = s, seed = sd, ll = as.numeric(logLik(f)),
      conv = f$opt$convergence)
  }
}
r <- do.call(rbind, rows)
utils::write.csv(r, paste0(
  "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/",
  "alphascale-", Sys.getenv("SKEWINIT_TAG", "fixed"), ".csv"),
  row.names = FALSE)
cat("   wrote", nrow(r), "fits\n")
options(op)
