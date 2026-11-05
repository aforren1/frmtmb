# Item 4: beta_binomial() with no overdispersion. Rounding or a real limit? Seeds as in
# dev/nlfix-bbscan2.R. FRM_ARM selects the build.
source("dev/review-nlfix-common.R")
suppressPackageStartupMessages(library(testthat))
ns <- asNamespace("frmtmb")
E <- new.env(parent = ns); E$`%||%` <- function(a, b) if (is.null(a)) b else a
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) sys.source(f, envir = E)
exact_ll <- function(y, size, p, phi) {
  # the finite product: no lgamma of a large argument anywhere
  a <- p * phi; b <- (1 - p) * phi
  vapply(seq_along(y), function(i) {
    s <- 0
    if (y[i] > 0) s <- s + sum(log(a[i] + 0:(y[i] - 1)))
    if (size[i] - y[i] > 0) s <- s + sum(log(b[i] + 0:(size[i] - y[i] - 1)))
    s <- s - sum(log(phi[i] + 0:(size[i] - 1)))
    lchoose(size[i], y[i]) + s
  }, 0)
}
dat <- function(s) { set.seed(s); n <- 300; d <- data.frame(x = round(rnorm(n), 3), ntr = 10L); d$y <- rbinom(n, 10, plogis(0.3 + 0.5 * d$x)); d }
cat("--- the fits\n")
for (s in c(2, 4, 6, 9, 12)) {
  d <- dat(s)
  a <- suppressWarnings(frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = d))
  b <- suppressWarnings(frm(bf(y | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = beta_binomial(), data = d))
  bi <- as.numeric(logLik(glm(cbind(y, ntr - y) ~ x, binomial, d)))
  p <- plogis(a$estimates$beta[1] + a$estimates$beta[2] * d$x)
  cat(sprintf("seed %2d linear %.4f code %d log phi %.2f | body %.4f code %d log phi %.2f | binomial limit %.4f | exact density at the linear estimate %.4f\n", s,
              as.numeric(logLik(a)), a$opt$convergence, a$estimates$betad, as.numeric(logLik(b)), b$opt$convergence, b$estimates$betad, bi,
              sum(exact_ll(d$y, d$ntr, p, rep(exp(a$estimates$betad), nrow(d))))))
}
cat("\n--- frm()'s objective against the exact density and brms's Stan program, one data set, a grid of log phi\n")
d <- dat(9)
fit0 <- suppressWarnings(frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = d))
tp <- ns$make_start(fit0$frame, NULL)
beta <- coef(glm(cbind(y, ntr - y) ~ x, binomial, d))
tp$beta[] <- beta
nll <- ns$build_objective(fit0$frame)
fb <- brms::bf(y | trials(ntr) ~ x, center = FALSE)
fam <- brms::brmsfamily("beta_binomial")
prior <- E$brms_flat_prior(fb, data = d, family = fam)
code <- brms::make_stancode(fb, data = d, family = fam, prior = prior)
sdat <- brms::make_standata(fb, data = d, family = fam, prior = prior)
sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 0))
cat(grep("target \\+=|phi", strsplit(as.character(code), "\n")[[1]], value = TRUE)[1:5], sep = "\n")
p <- plogis(beta[1] + beta[2] * d$x)
for (lp in c(2, 10, 15, 20, 25, 30, 35, 39, 45, 60)) {
  tp$betad[] <- lp
  ours <- -as.numeric(nll(tp))
  ex <- sum(exact_ll(d$y, d$ntr, p, rep(exp(lp), nrow(d))))
  st <- tryCatch(rstan::log_prob(sf, rstan::unconstrain_pars(sf, list(b = as.array(unname(beta)), phi = exp(lp))), adjust_transform = FALSE), error = function(e) NA)
  cat(sprintf("log phi %3d: frm %.6g | exact %.6f | brms Stan %.6f | frm - exact %.3g\n", lp, ours, ex, st, ours - ex))
}
cat("\n--- the change: lgamma(a + k) - lgamma(a) by lgamma_shift_diff(), as a family\n")
print(args(ns$lgamma_shift_diff))
bb2 <- frmtmb_family("bb2", dpars = c("mu", "phi"), links = list(mu = "logit", phi = "log"), type = "discrete",
  required_aterms = "trials",
  lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
    size <- ns$aterm_trials(aterms, "bb2")
    a <- dpars$mu * dpars$phi; b <- (1 - dpars$mu) * dpars$phi
    "c" <- RTMB::ADoverload("c")
    lgamma(size + 1) - lgamma(y + 1) - lgamma(size - y + 1) +
      ns$lgamma_shift_diff(a, y) + ns$lgamma_shift_diff(b, size - y) - ns$lgamma_shift_diff(dpars$phi, size)
  },
  init_dpars = list(mu = function(y, aterms) 0.5, phi = function(y, aterms) 5))
for (s in c(2, 4, 6, 9, 12)) {
  d <- dat(s)
  r <- tryCatch(withCallingHandlers(frm(bf(y | trials(ntr) ~ x), family = bb2, data = d), warning = function(w) invokeRestart("muffleWarning")), error = function(e) e)
  bi <- as.numeric(logLik(glm(cbind(y, ntr - y) ~ x, binomial, d)))
  r2 <- tryCatch(suppressWarnings(frm(bf(y | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = bb2, data = d)), error = function(e) e)
  cat(sprintf("seed %2d patched linear %s | body %s | binomial limit %.5f\n", s,
              if (inherits(r, "error")) substr(conditionMessage(r), 1, 140) else sprintf("%.5f code %d log phi %.2f", as.numeric(logLik(r)), r$opt$convergence, r$estimates$betad),
              if (inherits(r2, "error")) substr(conditionMessage(r2), 1, 80) else sprintf("%.5f code %d", as.numeric(logLik(r2)), r2$opt$convergence), bi))
}
cat("\n--- patched against the exact density, moderate phi, overdispersed data, and the same fit on the current density\n")
set.seed(21); d <- data.frame(x = rnorm(300), ntr = 10L); p <- plogis(0.3 + 0.5 * d$x); d$y <- rbinom(300, 10, rbeta(300, p * 8, (1 - p) * 8))
a <- suppressWarnings(frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = d)); b <- suppressWarnings(frm(bf(y | trials(ntr) ~ x), family = bb2, data = d))
cat("current", format(as.numeric(logLik(a)), digits = 12), a$estimates$betad, " patched", format(as.numeric(logLik(b)), digits = 12), b$estimates$betad, "\n")
cat("\n--- other families with the same form\n")
for (fm in c("zero_inflated_beta_binomial", "dirichlet_multinomial")) {
  res <- character(0)
  for (s in c(2, 4, 6, 9, 12)) {
    set.seed(s); n <- 300; dd <- data.frame(x = round(rnorm(n), 3), ntr = 10L)
    if (fm == "dirichlet_multinomial") {
      m <- t(sapply(seq_len(n), function(i) rmultinom(1, 10, c(1, exp(0.5 * dd$x[i]), 1)))); colnames(m) <- c("a", "b", "c"); dd$y <- m
    } else dd$y <- rbinom(n, 10, plogis(0.3 + 0.5 * dd$x))
    r <- tryCatch(suppressWarnings(frm(bf(y | trials(ntr) ~ x), family = get(fm)(), data = dd)), error = function(e) e)
    res <- c(res, if (inherits(r, "error")) "ERR" else sprintf("%.2f", as.numeric(logLik(r))))
  }
  cat(fm, ":", res, "\n")
}
cat("\n--- a positive logLik on a discrete family anywhere else\n")
set.seed(5); dp <- data.frame(x = rnorm(100)); dp$y <- rbinom(100, 1, plogis(40 * dp$x))
cat("bernoulli separation:", as.numeric(logLik(suppressWarnings(frm(y ~ x, family = bernoulli(), data = dp)))), "\n")
