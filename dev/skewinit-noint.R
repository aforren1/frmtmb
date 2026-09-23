# The blocker battery for punch round 2: no-intercept designs across
# families, plus the reviewer's ill-scaled poisson at four covariate
# scales. Every non-skew_normal job must be BITWISE identical to base,
# because only a dpar that declares a stationary point may have its
# start spread over an intercept-less design.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
TAG <- Sys.getenv("SKEWINIT_TAG", "fixed")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

set.seed(23)
n <- 250
g <- factor(rep(1:4, length.out = n))
x <- rnorm(n)
eta <- 0.4 + 0.8 * x
dd <- data.frame(
  x = x, g = g,
  ygau = eta + rnorm(n, 0, 0.9),
  ypos = rgamma(n, shape = 3, scale = exp(eta) / 3),
  ycnt = rpois(n, exp(0.5 + 0.4 * x)),
  ybin = rbinom(n, 1, plogis(eta)),
  ybeta = rbeta(n, 2 * plogis(eta) * 5, 2 * (1 - plogis(eta)) * 5),
  yzi = ifelse(rbinom(n, 1, 0.25) == 1, 0, rpois(n, exp(0.5 + 0.4 * x))),
  ysn = eta + 1.2 * (0.9486833 * abs(rnorm(n)) + 0.3162278 * rnorm(n))
)
# the reviewer's construction: the same response, the covariate rescaled
set.seed(23)
xt_base <- rnorm(n)
dd$ypois_t <- rpois(n, pmin(exp(0.5 + 0.4 * xt_base), 1e6))
scales <- c(s0 = 1, s2 = 1e-2, s4 = 1e-4, s6 = 1e-6)
for (nm in names(scales)) dd[[paste0("xt_", nm)]] <- xt_base * scales[[nm]]

jobs <- list(
  gaussian_noint   = list(bf(ygau ~ 0 + g), gaussian()),
  gaussian_sig_noi = list(bf(ygau ~ x, sigma ~ 0 + g), gaussian()),
  student_noint    = list(bf(ygau ~ 0 + g), student()),
  student_nu_noint = list(bf(ygau ~ x, nu ~ 0 + g), student()),
  gamma_noint      = list(bf(ypos ~ 0 + g), Gamma(link = "log")),
  gamma_shape_noi  = list(bf(ypos ~ x, shape ~ 0 + g), Gamma(link = "log")),
  exgaussian_noint = list(bf(ypos ~ 0 + g), exgaussian()),
  lognormal_noint  = list(bf(ypos ~ 0 + g), lognormal()),
  poisson_noint    = list(bf(ycnt ~ 0 + g), poisson()),
  negbin_noint     = list(bf(ycnt ~ 0 + g), negbinomial()),
  bernoulli_noint  = list(bf(ybin ~ 0 + g), bernoulli()),
  beta_noint       = list(bf(ybeta ~ 0 + g), Beta()),
  beta_phi_noint   = list(bf(ybeta ~ x, phi ~ 0 + g), Beta()),
  zi_noint         = list(bf(yzi ~ 0 + g), zero_inflated_poisson()),
  zi_zi_noint      = list(bf(yzi ~ x, zi ~ 0 + g), zero_inflated_poisson()),
  weibull_noint    = list(bf(ypos ~ 0 + g), weibull()),
  # skew_normal: alpha DECLARES, so alpha ~ 0 + g is expected to move;
  # sigma does not, so sigma ~ 0 + g must be identical to base
  sn_alpha_noint   = list(bf(ysn ~ x, sigma ~ 1, alpha ~ 0 + g),
                          skew_normal()),
  sn_sigma_noint   = list(bf(ysn ~ x, sigma ~ 0 + g, alpha ~ 1),
                          skew_normal()),
  sn_mu_noint      = list(bf(ysn ~ 0 + g, sigma ~ 1, alpha ~ 1),
                          skew_normal())
)
for (nm in names(scales)) {
  v <- paste0("xt_", nm)
  jobs[[paste0("pois_", nm)]] <-
    list(bf(stats::as.formula(paste("ypois_t ~ 0 +", v))), poisson())
}

out <- list()
for (nm in names(jobs)) {
  j <- jobs[[nm]]
  r <- tryCatch({
    nw <- 0L
    f <- withCallingHandlers(
      frm(j[[1]], family = j[[2]], data = dd),
      warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") })
    list(ll = as.numeric(logLik(f)), par = unname(as.numeric(f$opt$par)),
         conv = f$opt$convergence, msg = f$opt$message, warns = nw,
         grad = max(abs(f$obj$gr(f$opt$par))))
  }, error = function(e) list(error = conditionMessage(e)))
  out[[nm]] <- r
  cat(sprintf("%-18s %s\n", nm,
              if (is.null(r$error)) sprintf("%s conv %d warn %d",
                                            format(r$ll, digits = 15),
                                            r$conv, r$warns)
              else paste("ERROR", r$error)))
}
# the external oracle for the poisson cases
for (nm in names(scales)) {
  v <- paste0("xt_", nm)
  gl <- stats::glm(stats::as.formula(paste("ypois_t ~ 0 +", v)),
                   family = stats::poisson(), data = dd)
  out[[paste0("glm_", nm)]] <-
    list(ll = as.numeric(stats::logLik(gl)), par = unname(coef(gl)),
         conv = 0L, msg = "glm", warns = 0L, grad = NA_real_)
  cat(sprintf("%-18s %s  (glm oracle, scale %g)
", paste0("glm_", nm),
              format(as.numeric(stats::logLik(gl)), digits = 15),
              scales[[nm]]))
}
p <- paste0("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/",
            "noint-", TAG, ".rds")
saveRDS(out, p)
cat("wrote", p, "\n")
