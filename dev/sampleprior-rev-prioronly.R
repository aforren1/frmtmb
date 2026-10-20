# Reviewer of lane sampleprior: a PRIOR-ONLY tape built from exactly the
# entries frm_sample() resolves for its defaults (sample_resolve_priors()
# then neg_log_prior_fn(), as prior_augmented_obj() composes them, with
# no likelihood), every parameter but the tested dpars fixed by `map`,
# sampled by tmbstan. The draws of each natural parameter must follow
# brms's prior on it: one Jacobian too few or too many shows as a
# distribution off the target. Sampler seed 11, 4 chains x 4000.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
ns <- asNamespace("frmtmb.sample")
prior_only <- function(f, keep) {
  fit <- suppressMessages(suppressWarnings(ns$sample_assemble(
    f[[1]], f[[2]], f[[3]], data2 = list(), start = NULL,
    fit_control = frmtmb_control(), na.action = stats::na.omit,
    REML = FALSE)))
  rp <- suppressMessages(ns$sample_resolve_priors(fit, NULL, base = NULL,
                                                  defaults = TRUE))
  nlp <- asNamespace("frmtmb")$neg_log_prior_fn(rp$ri$entries)
  pars <- ns$mo_sample_pars(fit, fit$estimates)
  map <- lapply(pars, function(v) factor(rep(NA, length(v))))
  for (cp in names(keep)) {
    m <- rep(NA_integer_, length(pars[[cp]]))
    m[keep[[cp]]] <- seq_along(keep[[cp]])
    map[[cp]] <- factor(m)
  }
  obj <- RTMB::MakeADFun(function(p) nlp(p), pars, map = map, silent = TRUE)
  sf <- suppressWarnings(real_tmbstan(obj, chains = 4, iter = 4000,
                                      seed = 11, refresh = 0))
  a <- rstan::extract(sf, permuted = FALSE)
  m <- do.call(rbind, lapply(seq_len(dim(a)[2]), function(ch) a[, ch, ]))
  div <- sum(rstan::get_divergent_iterations(sf))
  list(m = m[, seq_len(ncol(m) - 1L), drop = FALSE], div = div)
}
report <- function(lab, v, cdf, q) {
  ks <- suppressWarnings(stats::ks.test(v, cdf))
  # the KS p-value overstates on autocorrelated draws, so the quantiles
  # are printed beside it, against the target's
  cat(sprintf("%-34s n %d KS D %.4f p %.3f | q05/50/95 draws %s target %s\n",
              lab, length(v), ks$statistic, ks$p.value,
              paste(sprintf("%.3f", stats::quantile(v, c(.05, .5, .95))),
                    collapse = " "),
              paste(sprintf("%.3f", q(c(.05, .5, .95))), collapse = " ")))
}
tg <- function(a, b, lo) {
  list(cdf = function(v) (pgamma(v, a, b) - pgamma(lo, a, b)) /
         (1 - pgamma(lo, a, b)),
       q = function(p) qgamma(pgamma(lo, a, b) + p * (1 - pgamma(lo, a, b)),
                              a, b))
}

# 1. student nu: gamma(2, 0.1) on nu > 1, through the logm1 link
r <- prior_only(eval(CASES$student$frm), list(betad = 2))
cat("divergences", r$div, "\n")
t1 <- tg(2, 0.1, 1)
report("student nu ~ gamma(2,0.1)[1,]", 1 + exp(r$m[, 1]), t1$cdf, t1$q)

# 2. zero_inflated_poisson zi without a formula: beta(1, 1), logit link
f <- list(y ~ x, D$zcnt, frmtmb::brmsfamily("zero_inflated_poisson"))
r <- prior_only(f, list(betad = 1))
cat("divergences", r$div, "\n")
report("zip zi ~ beta(1,1)", plogis(r$m[, 1]), punif, qunif)

# 3. negbinomial shape: inv_gamma(0.4, 0.3), log link
r <- prior_only(eval(CASES$negbin$frm), list(betad = 1))
cat("divergences", r$div, "\n")
report("negbin 1/shape ~ gamma(0.4, 0.3)", exp(-r$m[, 1]),
       function(v) pgamma(v, 0.4, 0.3), function(p) qgamma(p, 0.4, 0.3))

# 4. shifted_lognormal ndt: uniform(0, min_Y), log link, NaN past min_Y
f <- eval(CASES$shifted_ln$frm)
my <- min(f[[2]]$y)
r <- prior_only(f, list(betad = 2))
cat("divergences", r$div, " max ndt / min_Y", max(exp(r$m[, 1])) / my, "\n")
report("shifted_ln ndt ~ U(0, min_Y)", exp(r$m[, 1]) / my, punif, qunif)

# 5. wiener bs, ndt, bias: gamma(1, 1); U(0, min_Y) on scaled_logit;
#    beta(1, 1)
f <- eval(CASES$wiener$frm)
my <- min(f[[2]]$rt)
r <- prior_only(f, list(betad = 1:3))
cat("divergences", r$div, "\n")
report("wiener bs ~ gamma(1,1)", exp(r$m[, 1]),
       function(v) pgamma(v, 1, 1), function(p) qgamma(p, 1, 1))
report("wiener ndt ~ U(0, min_Y)", plogis(r$m[, 2]), punif, qunif)
report("wiener bias ~ beta(1,1)", plogis(r$m[, 3]), punif, qunif)

# 6. a dpar WITH a formula, phi ~ 1 of beta: student_t(3, 0, 2.5) on the
#    log-link intercept itself, no Jacobian
f <- list(bf(y ~ x, phi ~ 1), D$unit, frmtmb::brmsfamily("beta"))
r <- prior_only(f, list(betad = 1))
cat("divergences", r$div, "\n")
report("beta Intercept_phi ~ t(3,0,2.5)", r$m[, 1] / 2.5,
       function(v) pt(v, 3), function(p) qt(p, 3) * 2.5 / 2.5)

# 7. a mixture's weights: dirichlet(1), flat on the simplex
f <- eval(CASES$mixture$frm)
fit <- suppressMessages(suppressWarnings(ns$sample_assemble(
  f[[1]], f[[2]], f[[3]], data2 = list(), start = NULL,
  fit_control = frmtmb_control(), na.action = stats::na.omit, REML = FALSE)))
rp <- suppressMessages(ns$sample_resolve_priors(fit, NULL, base = NULL,
                                                defaults = TRUE))
fit$prior <- rp$effective
snlp <- ns$sample_extra_nlp(fit, rp$ri$entries)
pars <- ns$mo_sample_pars(fit, fit$estimates)
map <- lapply(pars, function(v) factor(rep(NA, length(v))))
map$betad <- factor(c(NA, NA, 1))
obj <- RTMB::MakeADFun(function(p) snlp(p), pars, map = map, silent = TRUE)
sf <- suppressWarnings(real_tmbstan(obj, chains = 4, iter = 4000, seed = 11,
                                    refresh = 0))
v <- plogis(as.numeric(rstan::extract(sf, permuted = FALSE)[, , 1]))
report("mixture theta2 ~ dirichlet(1)", v, punif, qunif)
