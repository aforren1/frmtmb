# Lane adhess: how far the change moves the covariance of ordinary,
# well-scaled fits. Per fit: the estimates, vcov(full = TRUE), the
# warnings, and a reference: the inverse of the AD Hessian without
# random effects, and a Richardson Jacobian of the AD gradient in
# standard-error units (numDeriv, six halvings, base step 0.01) with
# them. Compare the two arms with dev/adhess-movecmp.R.
#   Rscript dev/adhess-move.R <lane|base> <out.rds>
a <- commandArgs(TRUE)
Sys.setenv(ADHESS_ARM = a[1])
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
out <- a[2]
mk <- function(s, n = 300) {
  set.seed(s)
  d <- data.frame(x = stats::rnorm(n), z = stats::runif(n),
                  g = factor(rep(1:20, length.out = n)))
  u <- stats::rnorm(20, 0, 0.5)[d$g]
  us <- stats::rnorm(20, 0, 0.3)[d$g]
  d$y <- 1 + 0.5 * d$x + sin(3 * d$z) + u + us * d$x + stats::rnorm(n, 0, 0.4)
  d$yb <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * d$x + u))
  d$yc <- stats::rpois(n, exp(0.2 + 0.4 * d$x + u))
  d$yg <- stats::rgamma(n, 4, 4 / exp(0.5 + 0.3 * d$x))
  d$yn <- stats::rnbinom(n, mu = exp(1 + 0.3 * d$x + u), size = 2)
  d$yt <- 1 + 0.5 * d$x + stats::rt(n, 5)
  d$yo <- cut(d$y, c(-Inf, quantile(d$y, c(0.25, 0.5, 0.75)), Inf),
              ordered_result = TRUE)
  mu <- stats::plogis(-0.5 + 0.5 * d$x)
  d$ybeta <- stats::rbeta(n, mu * 20, (1 - mu) * 20)
  d$y2 <- 0.5 * d$y + stats::rnorm(n)
  d
}
specs <- list(
  lm = list(y ~ x + z, gaussian()),
  gamma = list(yg ~ x + z, Gamma(link = "log")),
  pois_glm = list(yc ~ x + z, poisson()),
  bern_glm = list(yb ~ x + z, bernoulli()),
  student = list(yt ~ x, student()),
  cumul = list(yo ~ x + z, cumulative()),
  beta = list(ybeta ~ x, Beta()),
  sigma_glm = list(bf(y ~ x, sigma ~ z), gaussian()),
  mv = list(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = TRUE), gaussian()),
  ri = list(y ~ x + z + (1 | g), gaussian()),
  rs = list(y ~ x + z + (1 + x | g), gaussian()),
  smooth = list(y ~ x + s(z) + (1 | g), gaussian()),
  sigma = list(bf(y ~ x + (1 | g), sigma ~ z), gaussian()),
  bern = list(yb ~ x + z + (1 | g), bernoulli()),
  pois = list(yc ~ x + z + (1 | g), poisson()),
  negbin = list(yn ~ x + (1 | g), negbinomial()),
  poly3 = list(y ~ poly(x, 3) + z + (1 | g), gaussian()))
ref_cov <- function(f) {
  obj <- f$obj
  p <- ns$sdr_outer_point(obj)
  st <- ns$obj_state_save(obj)
  on.exit(ns$obj_state_restore(obj, st))
  if (!length(obj$env$random)) {
    H <- tryCatch(obj$he(p), error = function(e) NULL)
    if (is.null(H)) return(NULL)
    return(tryCatch(solve((H + t(H)) / 2), error = function(e) NULL))
  }
  H0 <- ns$fd_hessian(p, obj$gr)$H
  D <- sqrt(abs(diag(H0)))
  g0 <- function(dz) as.numeric(obj$gr(p + dz / D)) / D
  J <- numDeriv::jacobian(g0, rep(0, length(p)), method = "Richardson",
                          method.args = list(eps = 0.01, d = 0.01,
                                             zero.tol = 1, r = 6, v = 2))
  J <- (J + t(J)) / 2 * outer(D, D)
  tryCatch(solve(J), error = function(e) NULL)
}
res <- list()
for (s in 1:3) {
  d <- mk(s)
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    r <- cap(frm(sp[[1]], family = sp[[2]], data = d))
    f <- r$v
    key <- paste(nm, s)
    if (is_err(f)) {
      res[[key]] <- list(error = as.character(f))
      cat(sprintf("%-10s seed %d ERROR %s\n", nm, s, f))
      next
    }
    V <- tryCatch(unname(suppressWarnings(vcov(f, full = TRUE))),
                  error = function(e) NULL)
    res[[key]] <- list(par = f$opt$par, V = V, warn = r$w,
                       lost = ns$sdr_of(f)$se_lost,
                       re = length(f$obj$env$random) > 0,
                       ref = if (arm == "lane") ref_cov(f))
    cat(sprintf("%-10s seed %d np %d warn %s lost %s\n", nm, s,
                length(f$opt$par), wtag(r$w), lost_of(f)))
  }
}
saveRDS(res, out)
