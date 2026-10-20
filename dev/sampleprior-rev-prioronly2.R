# Reviewer of lane sampleprior: the negbinomial shape's prior-only run
# again, longer and thinned so the KS test sees nearly independent draws
# (sampler seed 12, 4 chains x 22000, thin 10, warmup 2000)
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
ns <- asNamespace("frmtmb.sample")
f <- eval(CASES$negbin$frm)
fit <- suppressMessages(suppressWarnings(ns$sample_assemble(
  f[[1]], f[[2]], f[[3]], data2 = list(), start = NULL,
  fit_control = frmtmb_control(), na.action = stats::na.omit, REML = FALSE)))
rp <- suppressMessages(ns$sample_resolve_priors(fit, NULL, base = NULL,
                                                defaults = TRUE))
nlp <- asNamespace("frmtmb")$neg_log_prior_fn(rp$ri$entries)
pars <- ns$mo_sample_pars(fit, fit$estimates)
map <- lapply(pars, function(v) factor(rep(NA, length(v))))
map$betad <- factor(1)
obj <- RTMB::MakeADFun(function(p) nlp(p), pars, map = map, silent = TRUE)
# the density on log(shape) against the reviewer's own, on a grid
u <- seq(-4, 8, by = 0.5)
fr <- vapply(u, function(z) -obj$fn(z), 1)
tg <- 0.4 * log(0.3) - lgamma(0.4) - 1.4 * u - 0.3 * exp(-u) + u
cat("grid: max |(frm - target) - mean| =",
    max(abs((fr - tg) - mean(fr - tg))), "\n")
sf <- suppressWarnings(real_tmbstan(obj, chains = 4, iter = 22000,
                                    warmup = 2000, thin = 10, seed = 12,
                                    refresh = 0,
                                    control = list(adapt_delta = 0.95)))
a <- rstan::extract(sf, permuted = FALSE)[, , 1]
cat("divergences", sum(rstan::get_divergent_iterations(sf)),
    " ess_bulk", round(rstan::ess_bulk(a)), " rhat", round(rstan::Rhat(a), 4),
    "\n")
v <- exp(-as.numeric(a))
ks <- suppressWarnings(stats::ks.test(v, function(z) pgamma(z, 0.4, 0.3)))
cat(sprintf("negbin 1/shape n %d KS D %.4f p %.3f | q05/50/95 %s target %s\n",
            length(v), ks$statistic, ks$p.value,
            paste(sprintf("%.3f", quantile(v, c(.05, .5, .95))), collapse = " "),
            paste(sprintf("%.3f", qgamma(c(.05, .5, .95), 0.4, 0.3)),
                  collapse = " ")))
