# Designs near the basis's threshold (design_whitener() returns a block
# from a scaled condition number of 1e4): both sides must give the
# centered fit. REV_ARM selects the build.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
dw <- get("design_whitener", asNamespace("frmtmb"))
gd <- function(seed, m, fam, n = 150L, G = 15L) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), length.out = n))
  z <- rnorm(n)
  eta <- 0.5 + 0.4 * z + rnorm(G, 0, 0.7)[g]
  y <- switch(fam, poisson = rpois(n, exp(eta)),
              Gamma = rgamma(n, shape = 3, rate = 3 / exp(eta)),
              bernoulli = rbinom(n, 1, plogis(eta)),
              gaussian = eta + rnorm(n, 0, 0.5))
  data.frame(y = y, z = z, w = z + m, g = g)
}
ff <- function(s) switch(s, poisson = poisson(), Gamma = Gamma(link = "log"),
                         bernoulli = bernoulli(), gaussian = gaussian())
cat(sprintf("%-10s %4s %7s %6s %9s %11s %9s %9s %s\n", "family", "seed",
            "m", "basis", "kappa", "dll", "shift", "seratio", "warn"))
worst <- list()
for (fam in c("poisson", "Gamma", "bernoulli", "gaussian")) {
  for (m in c(3000, 4000, 4500, 5000, 5250, 5500, 5750, 6000, 6500, 7000,
              8000, 1e4, 2e4)) {
    for (s in 1:5) {
      d <- gd(s, m, fam)
      X <- cbind(1, d$w)
      R <- qr.R(qr(X / rep(sqrt(colSums(X^2) / nrow(X)), each = nrow(X))))
      kap <- 1 / rcond(R, triangular = TRUE)
      fz <- quiet(frm(bf(y ~ z + (1 | g)), family = ff(fam), data = d))
      r <- rfit(bf(y ~ w + (1 | g)), family = ff(fam), data = d)
      if (r$err) { cat(fam, s, m, "ERROR", conditionMessage(r$f), "\n"); next }
      ew <- quiet(fixef(r$f))["w", 1:2]; ez <- quiet(fixef(fz))["z", 1:2]
      dll <- ll_of(r$f) - ll_of(fz)
      sh <- (ew[1] - ez[1]) / ez[2]
      cat(sprintf("%-10s %4d %7g %6s %9.3g %+11.3e %+9.2e %9.6f %s\n", fam, s,
                  m, basis_on(r$f), kap, dll, sh, ew[2] / ez[2], short_w(r)))
    }
  }
}
