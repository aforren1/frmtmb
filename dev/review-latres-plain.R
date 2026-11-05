# Attack 5: non-autocor models must have the same fits and no time regression.
# Large plain models, three repeats, the frame and the objective clocked apart.
# Usage: REV_ARM=base|lane Rscript review-latres-plain.R <outfile>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
outf <- commandArgs(TRUE)[[1]]
set.seed(9)
n <- 200000L
g <- factor(sample(seq_len(2000), n, TRUE))
d <- data.frame(g = g, x = rnorm(n), z = rnorm(n), w = rnorm(n))
d$y <- 1 + .5 * d$x - .2 * d$z + rnorm(2000)[g] + rnorm(n)
d$cnt <- rpois(n, exp(.5 + .3 * d$x + rnorm(2000, 0, .3)[g]))
sp <- list(
  lm = quote(frm(y ~ x + z + w, data = d)),
  lmer = quote(frm(y ~ x + z + (1 | g), data = d)),
  pois = quote(frm(cnt ~ x + (1 | g), family = poisson(), data = d)),
  frame = quote(frm(cnt ~ x + (x | g), family = poisson(), data = d, dry_run = "frame")),
  prior = quote(default_prior(cnt ~ x + (x | g), family = poisson(), data = d)))
for (nm in names(sp)) {
  ts <- numeric(3)
  for (i in 1:3) {
    t0 <- proc.time()[[3]]
    v <- suppressWarnings(eval(sp[[nm]]))
    ts[i] <- proc.time()[[3]] - t0
  }
  ll <- if (inherits(v, "frmtmb_fit")) sprintf("%.10f", as.numeric(logLik(v))) else ""
  np <- if (inherits(v, "frmtmb_fit")) length(v$obj$env$par) else NA
  cat(rev_arm, nm, sprintf("%.2f %.2f %.2f", ts[1], ts[2], ts[3]), ll, np, "\n", file = outf, append = TRUE)
}
