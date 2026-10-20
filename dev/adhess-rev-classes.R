# Review of lane adhess, item 1: every fit class on the no-RE path.
# Per case: the fit's vcov() against the inverse of a Richardson
# Jacobian of the AD gradient at the covariance's own point, and the AD
# Hessian against that Jacobian.
#   REV_ARM=lane|base|merged Rscript dev/adhess-rev-classes.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-common.R")
suppressMessages({
  library(frmtmb.eam)
  library(frmtmb.coupling)
})
cases <- list()
cases$map_prior <- function() {
  set.seed(11)
  d <- data.frame(x = rnorm(150, 500, 100))
  d$y <- rpois(150, exp(0.3 + 0.004 * (d$x - 500)))
  frm(y ~ x, family = poisson(), data = d,
      prior = set_prior("normal(0, 0.001)", class = "b"))
}
cases$map_prior_unit <- function() {
  set.seed(12)
  d <- data.frame(x = rnorm(60))
  d$y <- rbinom(60, 1, plogis(0.3 + 1.5 * d$x))
  frm(y ~ x, family = bernoulli(), data = d,
      prior = set_prior("normal(0, 0.5)", class = "b"))
}
cases$bound_lb <- function() {
  set.seed(13)
  d <- data.frame(x = rnorm(100, 50, 10))
  d$y <- 1 + 0.05 * d$x + rnorm(100)
  frm(y ~ x, data = d, family = gaussian(),
      prior = set_prior("", class = "b", coef = "x", lb = 0.1))
}
cases$mo <- function() {
  set.seed(14)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 300, TRUE), levels = lev, ordered = TRUE)
  x <- rnorm(300, 200, 40)
  y <- rpois(300, exp(0.5 + c(0, 0.3, 0.5, 0.6)[income] + 0.002 * (x - 200)))
  frm(y ~ mo(income) + x, family = poisson(), data = data.frame(income, x, y))
}
cases$gev <- function() {
  set.seed(15)
  x <- rnorm(300, 30, 5)
  y <- extraDistr::rgev(300, 2 + 0.1 * (x - 30), 1, 0.15)
  frm(y ~ x, family = gen_extreme_value(), data = data.frame(x, y))
}
cases$mixture <- function() {
  set.seed(16)
  z <- rbinom(400, 1, 0.4)
  y <- ifelse(z == 1, rnorm(400, -1.5, 0.7), rnorm(400, 1.5, 1))
  frm(bf(y ~ 1, mu1 ~ 1, mu2 ~ 1), family = mixture(gaussian(), gaussian()),
      data = data.frame(y), start = list(beta = c(-1, 1)))
}
cases$cumulative <- function() {
  set.seed(17)
  x <- rnorm(400, 50, 10)
  yl <- 0.08 * (x - 50) + rlogis(400)
  y <- factor(cut(yl, c(-Inf, -1, 0, 1.2, Inf), labels = FALSE),
              ordered = TRUE)
  frm(y ~ x, family = cumulative(), data = data.frame(x, y))
}
cases$cens <- function() {
  set.seed(18)
  x <- rnorm(300, 100, 20)
  ys <- 1 + 0.02 * x + rnorm(300, 0, 0.5)
  cc <- ifelse(ys > 3.3, "right", "none")
  y <- pmin(ys, 3.3)
  frm(y | cens(cc) ~ x, data = data.frame(x, y, cc))
}
cases$trunc_pois <- function() {
  set.seed(19)
  x <- rnorm(400, 300, 60)
  y <- rpois(400, exp(0.2 + 0.004 * (x - 300)))
  d <- data.frame(x, y)[y >= 1, ]
  frm(y | trunc(lb = 1) ~ x, family = poisson(), data = d)
}
cases$ar <- function() {
  set.seed(20)
  n <- 200
  e <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
  x <- rnorm(n, 40, 8)
  d <- data.frame(t = seq_len(n), x, y = 1 + 0.05 * x + e)
  frm(y ~ x + ar(t), data = d)
}
cases$se_aterm <- function() {
  set.seed(21)
  x <- rnorm(80, 500, 100)
  s <- runif(80, 0.5, 2)
  y <- 0.1 + 0.002 * x + rnorm(80, 0, s)
  frm(y | se(s) ~ x, data = data.frame(x, y, s))
}
cases$asym_laplace <- function() {
  set.seed(22)
  x <- rnorm(300)
  y <- 1 + 0.5 * x + rt(300, 3)
  frm(y ~ x, family = asym_laplace(), data = data.frame(x, y))
}
cases$const_sigma <- function() {
  set.seed(23)
  x <- rnorm(200, 500, 100)
  y <- rgamma(200, shape = 3, rate = 3 / exp(0.5 + 0.003 * (x - 500)))
  frm(bf(y ~ x, shape = 3), family = Gamma(link = "log"),
      data = data.frame(x, y))
}
cases$wiener_ndt <- function() {
  set.seed(24)
  d <- ddm_simulate(300, mu = 1.2, bs = 1.4, ndt = 0.25, bias = 0.5)
  frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
      family = wiener(), data = d)
}
cases$coupling <- function() {
  tf <- paste0("C:/Users/adf44/source/r/frmtmb-wt-adhess/extensions/",
               "frmtmb.coupling/tests/testthat/test-coherence.R")
  ex <- parse(tf)
  for (e in ex[1:3]) eval(e, globalenv())
  set.seed(25)
  d <- cp_units(rep(2, 15L), rep(0.4, 15L), n = 6L)
  frm(cp_bf(), family = cross_wishart(), data = d)
}
cases$spline <- function() {
  set.seed(26)
  x <- runif(200, 0, 100)
  y <- sin(x / 15) + rnorm(200, 0, 0.3)
  frm(y ~ s(x, k = 8), data = data.frame(x, y))
}
cases$skew_normal <- function() {
  set.seed(27)
  x <- rnorm(300, 50, 10)
  y <- 1 + 0.05 * x + rnorm(300) + abs(rnorm(300))
  frm(y ~ x, family = skew_normal(), data = data.frame(x, y))
}

only <- commandArgs(TRUE)
if (length(only)) cases <- cases[only]
for (nm in names(cases)) {
  set.seed(1)
  r <- cap(cases[[nm]]())
  fit <- r$v
  if (inherits(fit, "rev_err")) {
    cat(sprintf("%-15s ERROR %s\n", nm, fit))
    next
  }
  nre <- length(fit$obj$env$random)
  at <- the_point(fit)
  V <- tryCatch(suppressWarnings(vcov(fit, full = TRUE)),
                error = function(e) NULL)
  lost <- names(sdr_of(fit)$se_lost %||% character(0))
  st <- ns$obj_state_save(fit$obj)
  Hx <- tryCatch(fit$obj$he(at), error = function(e) NULL)
  ns$obj_state_restore(fit$obj, st)
  Hr <- tryCatch(rich_hess(fit$obj, at), error = function(e) NULL)
  pen <- tryCatch(fit$obj$fn(at) + as.numeric(logLik(fit)),
                  error = function(e) NA)
  hd <- if (!is.null(Hx) && !is.null(Hr)) {
    D <- sqrt(abs(diag(Hr)))
    max(abs((Hx + t(Hx)) / 2 - Hr) / outer(D, D))
  } else NA
  Vr <- tryCatch(solve(Hr), error = function(e) NULL)
  e <- if (!is.null(V) && !is.null(Vr)) cov_err(V, Vr) else c(se = NA, cor = NA)
  bnd <- any(grepl("bound", c(r$m, r$w), ignore.case = TRUE))
  hc <- fit$cache$hessian_fixed
  cat(sprintf(paste0("%-15s np %2d re %d units %s | fn+logLik %.3g | ",
                     "ADH-Rich %.2e | SE err %.2e cor err %.2e | exact %s ",
                     "| lost %s | warn %d msg %d bnd %s\n"),
              nm, length(at), nre, !is.null(fit$par_units), pen, hd,
              e[["se"]], e[["cor"]], format(isTRUE(hc$exact)),
              if (length(lost)) paste(lost, collapse = ",") else "-",
              length(r$w), length(r$m), bnd))
  if (length(r$w)) cat("   warn:", substr(r$w, 1, 150), sep = "\n   ")
  if (length(r$m)) cat("   msg:", substr(r$m, 1, 150), sep = "\n   ")
  if (!is.null(V)) cat("   SE:", format(sqrt(diag(V)), digits = 8), "\n")
}
