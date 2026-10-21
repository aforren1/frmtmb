# Reviewer: Stan-1 reproduced independently. neg_binomial_2_log_glm_lpmf
# against neg_binomial_2_log_lpmf and an exact reference (rising
# factorial as a sum of logs) on a phi grid; then the lane's stuck draw,
# read from the lane's brms fit (read-only); then which brms programs
# write the GLM primitive.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
cat("rstan", format(packageVersion("rstan")), "StanHeaders",
    format(packageVersion("StanHeaders")), "\n")
code <- "functions {
  real g_ll(array[] int y, matrix X, real a, vector b, real phi) {
    return neg_binomial_2_log_glm_lpmf(y | X, a, b, phi);
  }
  real p_ll(array[] int y, vector eta, real phi) {
    return neg_binomial_2_log_lpmf(y | eta, phi);
  }
}
"
rstan::expose_stan_functions(rstan::stanc(model_code = code))
# exact log P(y) for NB2(mu, phi), y small: rising factorial by sums
exact <- function(y, eta, phi) {
  mu <- exp(eta)
  rf <- if (y > 0) sum(log(phi + seq_len(y) - 1)) else 0
  rf - lgamma(y + 1) + y * (eta - log(mu + phi)) -
    phi * log1p(mu / phi)
}
X <- matrix(0, 1, 1)
for (y in c(0L, 1L, 4L)) {
  for (eta in c(-1, 3)) {
    for (lp in c(2, 6, 8, 10, 12, 14, 16, 20)) {
      phi <- 10^lp
      g <- g_ll(y, X, eta, 0, phi)
      p <- p_ll(y, eta, phi)
      e <- exact(y, eta, phi)
      cat(sprintf(paste("y %d eta %4.1f phi 1e%02d glm %.12g plain %.12g",
                        "exact %.12g | rel glm %.2e plain %.2e\n"),
                  y, eta, lp, g, p, e, abs(g / e - 1), abs(p / e - 1)))
    }
  }
}
# a nonzero design: the GLM form with X b doing the work
X2 <- matrix(c(0.5, -1, 2), 3, 1)
for (lp in c(8, 12, 16)) {
  phi <- 10^lp
  y <- c(2L, 0L, 9L)
  eta <- 0.3 + X2[, 1] * 0.7
  cat(sprintf("3 rows, phi 1e%02d: glm %.12g plain %.12g exact %.12g\n", lp,
              g_ll(y, X2, 0.3, 0.7, phi), p_ll(y, eta, phi),
              sum(mapply(exact, y, eta, phi))))
}
b <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-out/brmsfit-nb2-glmbug.rds")
dr <- posterior::as_draws_df(b)
d <- b$data
s <- dr[dr$.chain == 2, ]
s <- s[nrow(s), ]
Xc <- matrix(d$x - mean(d$x), ncol = 1)
a <- s$b_Intercept + s$b_x * mean(d$x)
eta <- a + Xc[, 1] * s$b_x
cat(sprintf("stuck draw (chain 2, last): b_Intercept %.6g b_x %.6g sigma %.6g\n",
            s$b_Intercept, s$b_x, s$sigma))
cat(sprintf("  glm %.10g  plain %.10g  dnbinom %.10g  lp__ %.6g\n",
            g_ll(as.integer(d$y), Xc, a, s$b_x, 1 / s$sigma),
            p_ll(as.integer(d$y), eta, 1 / s$sigma),
            sum(dnbinom(d$y, size = 1 / s$sigma, mu = exp(eta), log = TRUE)),
            s$lp__))
cat("chain 2 distinct draws of b_Intercept:", length(unique(s$b_Intercept)),
    "of", nrow(dr[dr$.chain == 2, ]), "\n")
cat("chain 2 distinct draws (all):",
    length(unique(dr$b_Intercept[dr$.chain == 2])), "\n")
suppressMessages(library(brms))
dd <- data.frame(y = rpois(30, 3), x = rnorm(30), z = rnorm(30))
for (f in list(negbinomial(), brmsfamily("negbinomial2"))) {
  for (fo in list(y ~ x, y ~ 1, bf(y ~ x, shape ~ z), bf(y ~ x, sigma ~ z))) {
    sc <- tryCatch(stancode(fo, dd, family = f), error = function(e) "ERR")
    l <- grep("target \\+= neg_binomial", strsplit(sc, "\n")[[1]],
              value = TRUE)
    cat(f$family, deparse1(if (inherits(fo, "brmsformula")) fo$formula else
      fo), if (inherits(fo, "brmsformula")) deparse1(fo$pforms[[1]]), ":",
      if (length(l)) trimws(l) else "(refused or none)", "\n")
  }
}
