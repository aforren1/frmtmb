# More than 400 outer parameters with a group term, and importance =,
# on a shifted column. REV_ARM selects the build.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
line <- function(lbl, rw, fz) {
  if (rw$err) { cat(lbl, "ERROR", conditionMessage(rw$f), "\n"); return() }
  ew <- quiet(fixef(rw$f))["w", 1:2]; ez <- quiet(fixef(fz))["z", 1:2]
  cat(sprintf("%-40s basis=%s conv=%s dll=%+.3e shift=%+.2f SE warn: %s\n",
              lbl, basis_on(rw$f), rw$f$opt$convergence,
              ll_of(rw$f) - ll_of(fz), (ew[1] - ez[1]) / ez[2], short_w(rw)))
}
for (s in 1:2) {
  set.seed(300 + s)
  n <- 2400L
  f <- factor(rep(sprintf("L%03d", 1:400), length.out = n))
  g <- factor(rep(1:30, length.out = n))
  z <- rnorm(n)
  eta <- 0.3 + 0.4 * z + rnorm(400, 0, 0.3)[f] + rnorm(30, 0, 0.6)[g]
  d <- data.frame(y = rpois(n, exp(eta)), z = z, w = z + 1e6, f = f, g = g)
  t0 <- proc.time()[[3]]
  fz <- quiet(frm(bf(y ~ z + f + (1 | g)), family = poisson(), data = d))
  rw <- rfit(bf(y ~ w + f + (1 | g)), family = poisson(), data = d)
  line(sprintf("np=%d poisson (1|g) seed %d", length(rw$f$opt$par), 300 + s),
       rw, fz)
  cat("  time", proc.time()[[3]] - t0, "\n")
}
for (s in 1:2) {
  set.seed(s)
  n <- 150L; G <- 15L
  g <- factor(rep(seq_len(G), length.out = n))
  z <- rnorm(n)
  eta <- 0.5 + 0.4 * z + rnorm(G, 0, 0.7)[g]
  d <- data.frame(y = rpois(n, exp(eta)), z = z, w = z + 1e6, g = g)
  fz <- quiet(frm(bf(y ~ z + (1 | g)), family = poisson(), data = d,
                  importance = 200L))
  rw <- rfit(bf(y ~ w + (1 | g)), family = poisson(), data = d,
             importance = 200L)
  line(sprintf("importance = 200, seed %d", s), rw, fz)
}
