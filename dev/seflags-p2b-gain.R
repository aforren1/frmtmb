# Lane seflags, punch round 2b: how far short of the maximum do the B5
# fits stop, and does the Newton model see it? The log-likelihood
# against poly(x, 8)'s (the same maximum), against the Newton gain
# g' V g / 2 over the kept covariance and over the ill-conditioned
# directions, and after Newton steps over all of them.
#   Rscript dev/seflags-p2b-gain.R <lane|m9>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
for (s in c(3, 4, 10)) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)),
                  z = stats::rnorm(200), w = factor(sample(c("a", "b"), 200, TRUE)))
  d$y <- sin(3 * d$x) + 0.3 * d$z + 0.5 * (d$w == "b") +
    stats::rnorm(20, 0, 0.3)[d$g] + stats::rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  f <- suppressWarnings(frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + z + w +
                              (1 | g), data = d))
  fr <- suppressWarnings(frm(y ~ poly(x, 8) + z + w + (1 | g), data = d))
  sdr <- suppressWarnings(ns$sdr_of(f))
  p <- f$opt$par
  g <- as.numeric(f$obj$gr(p))
  Vk <- sdr$cov_fixed_prop
  il <- sdr$se_ill
  gk <- sum(g * (Vk %*% g)) / 2
  gi <- if (is.null(il)) 0 else sum((g %*% il$D)^2 / il$lam) / 2
  short <- f$obj$fn(p) - fr$opt$objective
  q <- p
  for (it in 1:3) {
    gq <- as.numeric(f$obj$gr(q))
    step <- Vk %*% gq + if (is.null(il)) 0 else il$D %*% (t(il$D) %*% gq / il$lam)
    q <- q - as.numeric(step)
  }
  cat(sprintf("seed %d code %d | short of the max by %.4f | Newton gain kept %.4f ill %.4f | after 3 steps short by %.4f\n",
              s, f$opt$convergence, short, gk, gi, f$obj$fn(q) - fr$opt$objective))
}
