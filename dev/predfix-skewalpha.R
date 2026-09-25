# Item 6, the declaring-dpar hazard: bf(y ~ xr, sigma ~ 1, alpha ~ 0 + xs)
# with xs a rescaled covariate (dev/skewinit-punch3.R part 2). Does the
# default reach the best of a predictor-scale start sweep?
#   PREDFIX_ARM=base|lane Rscript dev/predfix-skewalpha.R
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
mk <- function(seed, s, n = 250) {
  set.seed(seed)
  xs <- -abs(rnorm(n)) * 3
  yy <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  data.frame(y = yy, xs = xs * s, xr = xs)
}
fo <- bf(y ~ xr, sigma ~ 1, alpha ~ 0 + xs)
cat(sprintf("%-5s %-6s %16s %16s %12s %5s %s\n", "seed", "scale", "default",
            "best sweep", "short", "conv", "autoscaled"))
for (seed in 1:12) for (s in c(1, 1e-3, 1e-6, 1e-9)) {
  dd <- mk(seed, s)
  f <- suppressWarnings(frm(fo, family = skew_normal(), data = dd))
  sw <- vapply(c(-4, -2, -1, 1, 2, 4), function(v) {
    st <- list(betad = c(log(stats::sd(dd$y)), v / sqrt(mean(dd$xs^2))))
    as.numeric(logLik(suppressWarnings(
      frm(fo, family = skew_normal(), data = dd, start = st,
          control = frmtmb_control(autoscale = FALSE)))))
  }, 0)
  cat(sprintf("%-5d %-6g %16.7f %16.7f %12.6f %5d %s\n", seed, s,
              as.numeric(logLik(f)), max(sw), max(sw) - as.numeric(logLik(f)),
              f$opt$convergence, !is.null(f$par_units)))
}
