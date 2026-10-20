# Lane adhess: coefficient covariance of GLMs without random effects
# against the exact one, the inverse of the information X' W X formed
# and inverted in 256 bits at each fit's own estimate; glm()'s vcov()
# beside it, against the exact one at glm()'s estimate. Ten designs with
# covariates of very different scale, poisson and bernoulli, n = 500.
#   ADHESS_ARM=lane|base Rscript dev/adhess-glm.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
n <- 500
des <- list(
  d01_unit = function() data.frame(x = rnorm(n)),
  d02_x500 = function() data.frame(x = rnorm(n, 500, 100)),
  d03_sd100 = function() data.frame(x = rnorm(n, 0, 100)),
  d04_sd1e3 = function() data.frame(x = rnorm(n, 0, 1e-3)),
  d05_mixed2 = function() data.frame(x1 = rnorm(n), x2 = rnorm(n, 1e4, 1e3)),
  d06_shift = function() data.frame(x = runif(n, 1e3, 1e3 + 1)),
  d07_three = function() data.frame(x1 = rnorm(n), x2 = rnorm(n, 0, 1e3),
                                    x3 = rnorm(n, 0, 1e-2)),
  d08_factor = function() data.frame(f = factor(sample(letters[1:4], n, TRUE)),
                                     x = rnorm(n, 1e3, 200)),
  d09_poly3 = function() {
    x <- runif(n, 0, 10)
    data.frame(x = x, x2 = x^2, x3 = x^3)
  },
  d10_skew = function() data.frame(x = rgamma(n, shape = 2, scale = 5e3))
)
fam_frm <- list(poisson = poisson(), bernoulli = bernoulli())
fam_glm <- list(poisson = poisson(), bernoulli = binomial())
cat(sprintf("%-11s %-9s | %-28s | %-28s | %s\n", "design", "family",
            "frmtmb SE err / cor err", "glm SE err / cor err",
            "max |b_frm - b_glm| / SE, warn"))
for (k in seq_along(des)) {
  nm <- names(des)[k]
  for (fm in names(fam_frm)) {
    set.seed(100 * k + match(fm, names(fam_frm)))
    d <- des[[nm]]()
    fo <- stats::as.formula(paste("y ~", paste(setdiff(names(d), "y"),
                                               collapse = " + ")))
    X <- stats::model.matrix(fo[-2], d)
    z <- scale(X[, -1, drop = FALSE])
    b0 <- rep(0.5 / sqrt(ncol(z)), ncol(z))
    eta <- as.numeric(z %*% b0)
    d$y <- if (fm == "poisson") rpois(n, exp(0.5 + eta)) else
      rbinom(n, 1, plogis(0.2 + eta))
    r <- cap(frm(fo, family = fam_frm[[fm]], data = d))
    f <- r$v
    g <- suppressWarnings(glm(fo, family = fam_glm[[fm]], data = d,
                              control = glm.control(epsilon = 1e-12,
                                                    maxit = 100)))
    if (is_err(f)) {
      cat(sprintf("%-11s %-9s | ERROR %s\n", nm, fm, f))
      next
    }
    fe <- suppressWarnings(fixef(f))
    bf_ <- fe[, "Estimate"]
    Vf <- as.matrix(suppressWarnings(vcov(f)))
    Vxf <- glm_cov_exact(X, bf_, fm)
    Vxg <- glm_cov_exact(X, coef(g), fm)
    ef <- cov_err(Vf, Vxf)
    eg <- cov_err(unname(vcov(g)), Vxg)
    db <- max(abs(bf_ - coef(g)) / sqrt(diag(Vxf)))
    cat(sprintf("%-11s %-9s | %.2e / %.2e          | %.2e / %.2e          | %.1e, %s\n",
                nm, fm, ef[1], ef[2], eg[1], eg[2], db, wtag(r$w)))
  }
}
