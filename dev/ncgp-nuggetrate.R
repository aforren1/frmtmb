# How often the nugget warning fires, where the fit's sigma reads low and where it does
# not, on the warm-started build. Designs: dev/review-gpgap-ncgpfloor.R's seed 8901, n 60,
# and dev/review-ncgp-floor.R's, n 40 and 100, three length scales, three y scales. Then
# non-gaussian at large information. A row reads low when sigma is under half its true
# value; a noise-free row has no truth and is counted by itself.
#   Rscript dev/ncgp-nuggetrate.R
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
one <- function(tag, x, y, s, ysc) {
  r <- rev_fit(y ~ gp(x), data.frame(x = x, y = y))
  f <- r$fit
  if (inherits(f, "error")) return(invisible())
  g <- rev_gp_read(f)
  m <- !is.null(frmtmb::gp_nugget_message(f))
  nw <- sum(grepl("noise has reached the nugget", r$warn, fixed = TRUE))
  cat(sprintf("NR %-22s s %.0e sigma/s %s sigma/sdgp %.2e sdgp %.3g conv %d warns %d nugget-warning %d\n",
              tag, s, if (s > 0) sprintf("%.3g", g$sigma / ysc / s) else "NA",
              g$sigma / g$sd, g$sd / ysc, g$conv, length(r$warn), nw))
}
ratios <- c(1e-3, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 1e-8, 0)
for (A in c(1, 1000)) {
  set.seed(8901)
  x <- round(runif(60, 0, 6), 3)
  e <- rnorm(60)
  for (s in ratios) one(sprintf("gpgap A%g", A), x, A * (sin(2 * x) + s * e), s, A)
}
for (n in c(40L, 100L)) for (om in c(0.5, 2, 6)) for (ysc in c(1, 1000, 1e-3)) for (s in ratios) {
  set.seed(4100 + n)
  x <- sort(runif(n, 0, 10))
  e <- rnorm(n)
  one(sprintf("n%d om%g y%g", n, om, ysc), x, ysc * (sin(om * x) + s * e), s, ysc)
}
cat("== non-gaussian\n")
for (seed in 1:3) for (lm in c(4, 12, 16)) {
  set.seed(5601 + seed)
  x <- sort(round(runif(40, 0, 6), 3))
  d <- data.frame(x = x, y = rpois(40, exp(lm + sin(2 * x))))
  r <- rev_fit(y ~ gp(x), d, family = poisson())
  cat(sprintf("NRP poisson log mean %d seed %d conv %d warns %d nugget-warning %d\n", lm, seed, r$fit$opt$convergence,
              length(r$warn), sum(grepl("noise has reached the nugget", r$warn, fixed = TRUE))))
}
for (seed in 1:3) for (tr in c(50, 1e5, 1e7)) {
  set.seed(5701 + seed)
  x <- sort(round(runif(40, 0, 6), 3))
  d <- data.frame(x = x, t = tr)
  d$y <- rbinom(40, d$t, plogis(sin(2 * x)))
  r <- rev_fit(y | trials(t) ~ gp(x), d, family = binomial())
  cat(sprintf("NRB binomial trials %g seed %d conv %d warns %d nugget-warning %d\n", tr, seed, r$fit$opt$convergence,
              length(r$warn), sum(grepl("noise has reached the nugget", r$warn, fixed = TRUE))))
}
# positions repeating at a known noise: what the message says for sdgp
for (s in c(1e-3, 1e-5, 1e-6, 1e-7)) {
  set.seed(8902)
  u <- sort(runif(20, 0, 6)); x <- rep(u, each = 5)
  d <- data.frame(x = x, y = sin(2 * x) + s * rnorm(100), s = s)
  r <- rev_fit(y | se(s) ~ gp(x), d)
  g <- rev_gp_read(r$fit)
  cat(sprintf("NRS se() replicated s %.0e sdgp %.3g conv %d warns %d nugget-warning %d\n", s, g$sd, g$conv, length(r$warn),
              sum(grepl("noise has reached the nugget", r$warn, fixed = TRUE))))
}
