# Nothing but skew_normal may move. Fits a spread of families and fit
# modes, writes every number at full precision, and a second run
# compares the two files with identical().
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
TAG <- Sys.getenv("SKEWINIT_TAG", "fixed")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

set.seed(11)
n <- 300
x <- rnorm(n)
z <- runif(n)
g <- factor(rep(1:20, length.out = n))
gg <- rnorm(20, 0, 0.7)[as.integer(g)]
eta <- 0.4 + 0.8 * x
dd <- data.frame(
  x = x, z = z, g = g,
  ygau = eta + gg + rnorm(n, 0, 0.9),
  ypos = rgamma(n, shape = 3, scale = exp(eta) / 3),
  ycnt = rpois(n, exp(0.5 + 0.4 * x)),
  ybin = rbinom(n, 1, plogis(eta)),
  ybeta = rbeta(n, 2 * plogis(eta) * 5, 2 * (1 - plogis(eta)) * 5),
  ysn0 = 2 + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5,
  ysn1 = eta + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5,
  # a FINITE alpha (3) with a symmetric covariate, so the base build
  # does not stall here and the two optima can be compared closely
  ysn2 = eta + 1.2 * (0.9486833 * abs(rnorm(n)) +
                        0.3162278 * rnorm(n))
)

jobs <- list(
  gaussian_ml   = list(bf(ygau ~ x), gaussian(), list()),
  gaussian_re   = list(bf(ygau ~ x + (1 | g)), gaussian(), list()),
  gaussian_reml = list(bf(ygau ~ x + (1 | g)), gaussian(), list(REML = TRUE)),
  gaussian_prof = list(bf(ygau ~ x + (1 | g)), gaussian(),
                       list(control = frmtmb_control(profile = TRUE))),
  gaussian_sig  = list(bf(ygau ~ x, sigma ~ x), gaussian(), list()),
  student_ml    = list(bf(ygau ~ x), student(), list()),
  gamma_log     = list(bf(ypos ~ x), Gamma(link = "log"), list()),
  lognormal_ml  = list(bf(ypos ~ x), lognormal(), list()),
  weibull_ml    = list(bf(ypos ~ x), weibull(), list()),
  invgauss_ml   = list(bf(ypos ~ x), inverse.gaussian(link = "log"), list()),
  exgaussian_ml = list(bf(ypos ~ x), exgaussian(), list()),
  asymlap_ml    = list(bf(ygau ~ x), asym_laplace(), list()),
  poisson_ml    = list(bf(ycnt ~ x), poisson(), list()),
  negbin_ml     = list(bf(ycnt ~ x), negbinomial(), list()),
  bernoulli_ml  = list(bf(ybin ~ x), bernoulli(), list()),
  beta_ml       = list(bf(ybeta ~ x), Beta(), list()),
  poisson_re    = list(bf(ycnt ~ x + (1 | g)), poisson(), list()),
  poisson_prof  = list(bf(ycnt ~ x + (1 | g)), poisson(),
                       list(control = frmtmb_control(profile = TRUE))),
  # skew_normal WITHOUT a mu covariate: the residual is the centred
  # response, so nothing about its start may move either
  skew_icpt     = list(bf(ysn0 ~ 1, sigma ~ 1, alpha ~ 1), skew_normal(),
                       list()),
  skew_icpt_re  = list(bf(ysn0 ~ 1 + (1 | g), sigma ~ 1, alpha ~ 1),
                       skew_normal(), list()),
  # skew_normal WITH one: this one is expected to move
  skew_cov      = list(bf(ysn1 ~ x, sigma ~ 1, alpha ~ 1), skew_normal(),
                       list()),
  skew_cov_reml = list(bf(ysn1 ~ x, sigma ~ 1, alpha ~ 1), skew_normal(),
                       list(REML = TRUE)),
  skew_cov_prof = list(bf(ysn1 ~ x, sigma ~ 1, alpha ~ 1), skew_normal(),
                       list(control = frmtmb_control(profile = TRUE))),
  # finite alpha, no stall on either build
  skew_fin      = list(bf(ysn2 ~ x, sigma ~ 1, alpha ~ 1), skew_normal(),
                       list()),
  skew_fin_reml = list(bf(ysn2 ~ x, sigma ~ 1, alpha ~ 1), skew_normal(),
                       list(REML = TRUE)),
  skew_fin_prof = list(bf(ysn2 ~ x, sigma ~ 1, alpha ~ 1), skew_normal(),
                       list(control = frmtmb_control(profile = TRUE))),
  skew_fin_sig  = list(bf(ysn2 ~ x, sigma ~ x, alpha ~ 1), skew_normal(),
                       list())
)

out <- list()
for (nm in names(jobs)) {
  j <- jobs[[nm]]
  r <- tryCatch({
    f <- do.call(frm, c(list(j[[1]], family = j[[2]], data = dd), j[[3]]))
    list(ll = as.numeric(logLik(f)), obj = f$opt$objective,
         par = unname(as.numeric(f$opt$par)),
         se = unname(sqrt(diag(vcov(f)))),
         conv = f$opt$convergence)
  }, error = function(e) list(error = conditionMessage(e)))
  out[[nm]] <- r
  cat(sprintf("%-14s %s\n", nm,
              if (is.null(r$error)) format(r$ll, digits = 15) else
                paste("ERROR", r$error)))
}
p <- paste0("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/",
            "nomove-", TAG, ".rds")
saveRDS(out, p)
cat("wrote", p, "\n")
