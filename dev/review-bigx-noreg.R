# No-regression: well-conditioned designs must be bitwise unchanged, and
# ill-conditioned ones base fitted right must agree. Saves an RDS per
# arm; review-bigx-noreg-cmp.R compares. REV_ARM selects the build,
# REV_TAG names the output.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
tag <- Sys.getenv("REV_TAG", rev_arm)
set.seed(7001)
n <- 240L
g <- factor(rep(1:16, length.out = n)); h <- factor(rep(1:6, each = 40))
x <- rnorm(n); x2 <- rnorm(n); f <- factor(sample(c("a", "b", "c"), n, TRUE))
t <- runif(n); o <- sample(1:4, n, TRUE); wt <- runif(n, 0.5, 2)
u <- rnorm(16, 0, 0.6)[g]; us <- rnorm(16, 0, 0.3)[g]
eta <- 0.3 + 0.5 * x - 0.3 * x2 + u + us * x
d <- data.frame(g, h, x, x2, f, t, o, wt,
                yg = eta + rnorm(n, 0, 0.7),
                yp = rpois(n, exp(eta)),
                yb = rbinom(n, 1, plogis(eta)),
                yG = rgamma(n, 3, 3 / exp(eta)),
                ynb = rnbinom(n, mu = exp(eta), size = 2),
                ybeta = plogis(eta + rnorm(n, 0, 0.5)),
                yo = cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                         labels = FALSE),
                yt = eta + rt(n, 4) * 0.7,
                yz = rpois(n, exp(eta)) * rbinom(n, 1, 0.7),
                ycat = sample(1:3, n, TRUE),
                year = sample(1990:2020, n, TRUE))
d$yyr <- rpois(n, exp(0.2 + 0.03 * (d$year - 2005) + u))
d$xr <- runif(n, 1, 20)
d$ypr <- 1 + 0.2 * d$xr - 0.01 * d$xr^2 + rnorm(n)
d$xs <- 100 * x + 3e4  # kappa about 6e2 * 5 ... ill but base-fittable
d$ys <- rpois(n, exp(0.3 + 0.005 * (d$xs - 3e4) + u))
specs <- list(
  lm = list(bf(yg ~ x + x2), gaussian()),
  lmf = list(bf(yg ~ x * f), gaussian()),
  lmm = list(bf(yg ~ x + (1 | g)), gaussian()),
  lmm_rs = list(bf(yg ~ x + (1 + x | g)), gaussian()),
  lmm_reml = list(bf(yg ~ x + (1 | g)), gaussian(), TRUE),
  lmm_2re = list(bf(yg ~ x + (1 | g) + (1 | h)), gaussian()),
  pois = list(bf(yp ~ x + x2), poisson()),
  glmm_p = list(bf(yp ~ x + (1 | g)), poisson()),
  glmm_prs = list(bf(yp ~ x + (1 + x || g)), poisson()),
  glmm_pf = list(bf(yp ~ x + f + (1 | g)), poisson()),
  glmm_pw = list(bf(yp | weights(wt) ~ x + (1 | g)), poisson()),
  bern = list(bf(yb ~ x + x2), bernoulli()),
  glmm_b = list(bf(yb ~ x + (1 | g)), bernoulli()),
  gamma = list(bf(yG ~ x + (1 | g)), Gamma(link = "log")),
  nb = list(bf(ynb ~ x + (1 | g)), negbinomial()),
  beta = list(bf(ybeta ~ x + (1 | g)), Beta()),
  student = list(bf(yt ~ x + (1 | g)), student()),
  cum = list(bf(yo ~ x + (1 | g)), cumulative()),
  acat = list(bf(yo ~ x), acat()),
  sratio_cs = list(bf(yo ~ cs(x)), sratio()),
  zip = list(bf(yz ~ x + (1 | g), zi ~ x2), zero_inflated_poisson()),
  hurdle = list(bf(yz ~ x, hu ~ 1), hurdle_poisson()),
  sigma = list(bf(yg ~ x, sigma ~ x2), gaussian()),
  smooth = list(bf(yg ~ s(t) + x), gaussian()),
  smooth_p = list(bf(yp ~ s(t) + (1 | g)), poisson()),
  mo = list(bf(yp ~ mo(o) + x), poisson()),
  cat = list(bf(ycat ~ x), categorical()),
  offset = list(bf(yp ~ x + offset(log(wt))), poisson()),
  nl = list(bf(yg ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
            gaussian()),
  lognormal = list(bf(yG ~ x + (1 | g)), lognormal()),
  interact = list(bf(yg ~ x * x2 + (1 | g)), gaussian()),
  weibull = list(bf(yG ~ x), weibull()),
  # ill-conditioned designs base fitted right
  ill_year = list(bf(yyr ~ year + (1 | g)), poisson()),
  ill_rawpoly = list(bf(ypr ~ poly(xr, 2, raw = TRUE)), gaussian()),
  ill_shift = list(bf(ys ~ xs + (1 | g)), poisson()),
  ill_year_g = list(bf(yg ~ year), gaussian())
)
out <- list()
for (nm in names(specs)) {
  sp <- specs[[nm]]
  r <- rfit(sp[[1]], family = sp[[2]], data = d,
            REML = if (length(sp) > 2) sp[[3]] else FALSE)
  if (r$err) {
    out[[nm]] <- list(err = conditionMessage(r$f)); next
  }
  f <- r$f
  se <- tryCatch(sqrt(diag(quiet(frmtmb:::sdr_of(f)$cov.fixed))),
                 error = function(e) NA)
  out[[nm]] <- list(par = f$opt$par, obj = f$opt$objective,
                    conv = f$opt$convergence, se = se, w = r$w,
                    basis = basis_on(f), ll = ll_of(f))
  cat(sprintf("%-12s conv=%s basis=%s obj=%.10f\n", nm, f$opt$convergence,
              basis_on(f), f$opt$objective))
}
saveRDS(out, sprintf("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-noreg-%s.rds", tag))
