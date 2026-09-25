# Punch round 1, M1 calibration part 3: do random-slope shortfalls
# persist at spreads 0.1 to 1, where no scaling is involved? 16 seeds.
# slope on that column start to fall short with autoscale = FALSE, and is
# autoscale = TRUE (with the Z scaling) ever short? Reference: the same
# data at spread 1, best of the two arms.
#   PREDFIX_ARM=lane Rscript dev/predfix-p1-slopecal.R > dev/predfix-log/p1-slopecal.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
ll <- function(f) if (inherits(f, "try-error")) NA else as.numeric(logLik(f))
fit <- function(fo, fam, d, a) {
  try(suppressWarnings(frm(fo, family = fam, data = d,
                           control = frmtmb_control(autoscale = a))),
      silent = TRUE)
}
mk <- function(seed, tau, fam) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x <- rnorm(300)
  eta <- 0.5 + 0.4 * x + rnorm(20, 0, 0.6)[g] + rnorm(20, 0, tau)[g] * x
  y <- if (fam == "gau") eta + rnorm(300) else rpois(300, exp(eta))
  data.frame(g = g, x0 = x, y = y)
}
fo <- bf(y ~ x + (1 + x | g))
rows <- list()
for (fam in c("gau", "pois")) for (tau in c(0.1, 0.4)) for (seed in 1:16) {
  d <- mk(seed, tau, fam)
  fm <- if (fam == "gau") gaussian() else poisson()
  d$x <- d$x0
  ref <- max(ll(fit(fo, fm, d, FALSE)), ll(fit(fo, fm, d, TRUE)), na.rm = TRUE)
  for (s in c(1, 0.5, 0.2, 0.1)) {
    d$x <- d$x0 * s
    rows[[length(rows) + 1L]] <- data.frame(
      fam = fam, tau = tau, seed = seed, s = s, sdx = sd(d$x),
      short_false = ref - ll(fit(fo, fm, d, FALSE)),
      short_true = ref - ll(fit(fo, fm, d, TRUE)))
  }
}
r <- do.call(rbind, rows)
saveRDS(r, "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/p1-slopecal3.rds")
options(width = 200)
cat("autoscale = FALSE short by > 1e-6, by spread s (64 fits per s)\n")
print(tapply(r$short_false > 1e-6, r$s, sum))
cat("max shortfall by s\n")
print(signif(tapply(r$short_false, r$s, max), 4))
cat("autoscale = TRUE short by > 1e-6, by s\n")
print(tapply(r$short_true > 1e-6, r$s, sum))
cat("NA: FALSE", sum(is.na(r$short_false)), " TRUE", sum(is.na(r$short_true)), "\n")
print(r[r$short_false > 1e-6 & r$s >= 0, ], digits = 6)
