# Fits whose scaled design kappa is 1e3 to 1e4 (the basis now engages
# there): log-likelihood per cell, saved per arm, plus a reference on the
# centered or orthogonal spelling. REV_ARM selects the build, REV_TAG names it.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
tag <- Sys.getenv("REV_TAG", rev_arm)
kap <- function(X) {
  X <- as.matrix(X)
  R <- qr.R(qr(X / rep(sqrt(colSums(X^2) / nrow(X)), each = nrow(X))))
  1 / rcond(R, triangular = TRUE)
}
fams <- list(poisson = poisson(), bernoulli = bernoulli(), gaussian = gaussian(),
             Gamma = Gamma(link = "log"), negbinomial = negbinomial(),
             Beta = Beta(), student = student(), cumulative = cumulative())
gen <- function(fam, eta, n) switch(fam,
  poisson = rpois(n, exp(eta)), bernoulli = rbinom(n, 1, plogis(eta)),
  gaussian = eta + rnorm(n, 0, 0.6), Gamma = rgamma(n, 3, 3 / exp(eta)),
  negbinomial = rnbinom(n, mu = exp(eta), size = 3),
  Beta = plogis(eta + rnorm(n, 0, 0.5)), student = eta + 0.6 * rt(n, 5),
  cumulative = cut(eta + rlogis(n), c(-Inf, 0, 0.8, 1.6, Inf), labels = FALSE))
out <- list()
add <- function(id, k, r, rref) {
  out[[id]] <<- list(k = k, err = r$err, ll = if (r$err) NA else ll_of(r$f),
                     basis = if (r$err) NA else basis_on(r$f),
                     warned = length(r$w) > 0,
                     llref = if (rref$err) NA else ll_of(rref$f))
}
for (s in 1:4) {
  for (fam in names(fams)) {
    for (m in c(400, 700, 1000, 1500, 2500, 4000)) {
      for (re in c(FALSE, TRUE)) {
        set.seed(1000 * s + m)
        n <- 160; G <- 16
        g <- factor(rep(seq_len(G), length.out = n)); z <- rnorm(n)
        eta <- 0.4 + 0.5 * z + if (re) rnorm(G, 0, 0.6)[g] else 0
        d <- data.frame(y = gen(fam, eta, n), z = z, w = z + m, g = g)
        fw <- if (re) bf(y ~ w + (1 | g)) else bf(y ~ w)
        fz <- if (re) bf(y ~ z + (1 | g)) else bf(y ~ z)
        r <- rfit(fw, family = fams[[fam]], data = d)
        rr <- rfit(fz, family = fams[[fam]], data = d)
        add(sprintf("%s m=%g re=%s s=%d", fam, m, re, s), kap(cbind(1, d$w)), r, rr)
      }
    }
  }
  # raw polynomials of degree 3 on [1, 2] and [0, 3], and year-like columns
  set.seed(5000 + s); n <- 200
  x <- runif(n, 1, 2); g <- factor(rep(1:20, length.out = n))
  u <- rnorm(20, 0, 0.4)[g]
  d <- data.frame(x = x, g = g, y = sin(3 * x) + u + rnorm(n, 0, 0.2),
                  yp = rpois(n, exp(sin(3 * x) + u)))
  for (deg in 2:3) for (re in c(FALSE, TRUE)) for (fam in c("gaussian", "poisson")) {
    lhs <- if (fam == "gaussian") "y" else "yp"
    fr <- as.formula(sprintf("%s ~ poly(x, %d, raw = TRUE)%s", lhs, deg, if (re) " + (1 | g)" else ""))
    fo <- as.formula(sprintf("%s ~ poly(x, %d)%s", lhs, deg, if (re) " + (1 | g)" else ""))
    r <- rfit(bf(fr), family = fams[[fam]], data = d)
    rr <- rfit(bf(fo), family = fams[[fam]], data = d)
    add(sprintf("rawpoly%d %s re=%s s=%d", deg, fam, re, s),
        kap(model.matrix(lme4::nobars(fr), d)), r, rr)
  }
  d$yr <- sample(1990:2020, n, TRUE)
  d$yy <- rpois(n, exp(0.2 + 0.03 * (d$yr - 2005) + u))
  d$yc <- d$yr - 2005
  r <- rfit(bf(yy ~ yr + I(yr^2) + (1 | g)), family = poisson(), data = d)
  rr <- rfit(bf(yy ~ yc + I(yc^2) + (1 | g)), family = poisson(), data = d)
  add(sprintf("year quad s=%d", s), kap(cbind(1, d$yr, d$yr^2)), r, rr)
}
saveRDS(out, sprintf("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-fx-band-%s.rds", tag))
cat("cells", length(out), "\n")
