# The two ways to get prior draws, measured on model B of
# dev/priordraws-brms-names.R (a cauchy(0, 1) sd among its priors):
# direct draws from each prior (what sample_prior = "yes" does) against
# a run of the sampler on the prior-only objective (what "only" does).
# Effective sample size per draw, the exact-prior KS on a thinned run,
# divergences, and time. Seed 21, 4 chains of 2000.
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb.sample"), "\n")
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)))
d$y <- 1 + 0.5 * d$x + rnorm(60)
pr <- set_prior("normal(0, 1)", class = "b") +
  set_prior("normal(0, 2)", class = "b", coef = "z") +
  set_prior("normal(1, 3)", class = "Intercept") +
  set_prior("exponential(2)", class = "sigma") +
  set_prior("lkj(2)", class = "cor") +
  set_prior("cauchy(0, 1)", class = "sd", group = "g", coef = "x")
f <- bf(y ~ x + z + (1 + x | g))
q <- function(e) suppressWarnings(suppressMessages(e))
t_only <- system.time(do <- q(frm_sample(f, data = d, family = gaussian(),
                                         prior = pr, sample_prior = "only",
                                         chains = 4, iter = 2000, seed = 21,
                                         refresh = 0)))
t_no <- system.time(dn <- q(frm_sample(f, data = d, family = gaussian(),
                                       prior = pr, chains = 4, iter = 2000,
                                       seed = 21, refresh = 0)))
t_yes <- system.time(dy <- q(frm_sample(f, data = d, family = gaussian(),
                                        prior = pr, sample_prior = "yes",
                                        chains = 4, iter = 2000, seed = 21,
                                        refresh = 0)))
n <- ndraws(do)
div <- sum(subset(nuts_params(do), Parameter == "divergent__")$Value)
cat(sprintf("only run: %d draws, %.1f s elapsed, %d divergent\n", n,
            t_only[["elapsed"]], div))
cat(sprintf("posterior run without prior draws %.1f s, with %.1f s\n",
            t_no[["elapsed"]], t_yes[["elapsed"]]))
# the direct draws alone, timed over 20 repetitions
plan <- frmtmb.sample:::sample_prior_plan(
  dy$fit, list(effective = dy$fit$prior), "yes")
tm <- system.time(for (i in 1:20) {
  frmtmb.sample:::prior_draws_make(plan, n, dy$fit$prior)
})
cat(sprintf("direct draws of %d rows: %.4f s each (20 repetitions)\n", n,
            tm[["elapsed"]] / 20))
ref <- list(b_x = pnorm, b_z = function(v) pnorm(v, 0, 2),
            sigma = function(v) pexp(v, 2),
            sd_g__x = function(v) pmax(2 * pcauchy(v) - 1, 0),
            sd_g__Intercept = function(v) pmax(2 * pt(v / 2.5, 3) - 1, 0))
mo <- as.matrix(do)
pd <- prior_draws(dy)
for (v in names(ref)) {
  a <- matrix(mo[, v], ncol = 4)
  eb <- posterior::ess_bulk(a) / n
  et <- posterior::ess_tail(a) / n
  # KS on every 10th draw of the run, to thin the autocorrelation away
  thin <- mo[seq(1, n, by = 10), v]
  ksr <- ks.test(thin, ref[[v]])$p.value
  ksd <- ks.test(pd[[v]] %||% pd[[sub("_[xz]$", "", v)]], ref[[v]])$p.value
  cat(sprintf("  %-16s run: ess_bulk/n %.2f ess_tail/n %.2f KS(thin 10) p=%.3f | direct: ess/n 1 by construction, KS p=%.3f\n",
              v, eb, et, ksr, ksd))
}
qs <- c(0.99, 0.999)
cat("  sd_g__x upper quantiles, exact:", round(tan(pi * qs / 2), 1),
    "| run:", round(quantile(mo[, "sd_g__x"], qs), 1),
    "| direct:", round(quantile(pd$sd_g__x, qs), 1), "\n")
