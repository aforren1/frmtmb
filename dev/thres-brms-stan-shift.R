# Lane thres: brms's compiled cumulative log density against frmtmb's
# objective at the ML estimate and at a point moved off it (every
# parameter shifted by 0.03 times its position), to show that the exact
# agreement in dev/thres-brms-stan-log.txt is not a property of the
# optimum alone. Data seed 11, as there.
#
#   Rscript dev/thres-brms-stan-shift.R > dev/thres-brms-stan-shift-log.txt
.libPaths(c("/opt/rlib/stan", "/opt/rlib/lane-thres", "/opt/rlib/base",
            "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages({library(frmtmb); library(brms); library(rstan)})
set.seed(11)
n <- 240
d <- data.frame(x = rnorm(n), g = sample(c("a", "b", "c"), n, TRUE))
u <- rlogis(n, 0.7 * d$x)
tau <- list(a = c(-1, 0.2, 1.3), b = c(-0.5, 0.6),
            c = c(-1.5, -0.4, 0.5, 1.6))
d$y <- vapply(seq_len(n), function(i) 1L + sum(u[i] > tau[[d$g[i]]]), 1L)
fit <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative())
bf0 <- bf(y | thres(gr = g) ~ x)
pr <- get_prior(bf0, data = d, family = brms::cumulative())
pr$prior <- ""
code <- make_stancode(bf0, data = d, family = brms::cumulative(),
                      prior = pr)
sdat <- make_standata(bf0, data = d, family = brms::cumulative(),
                      prior = pr)
sf <- suppressMessages(sampling(stan_model(model_code = code),
                                data = sdat, chains = 0))
lay <- frmtmb:::thres_layout(c(3, 2, 4))
for (shift in c(0, 0.3)) {
  p <- fit$opt$par + shift * seq_along(fit$opt$par) / 10
  th <- fit$spec$responses$y$family$post$ord_thresholds(p[-1])
  pars <- list(b = array(p[1], 1))
  for (k in 1:3) {
    pars[[paste0("Intercept_", k)]] <- as.array(th[lay$start[k]:lay$end[k]])
  }
  lp <- log_prob(sf, unconstrain_pars(sf, pars), adjust_transform = FALSE)
  ours <- -fit$obj$fn(p)
  cat(sprintf("shift %.1f stan %.17g frmtmb %.17g diff %.3g\n", shift,
              lp, ours, lp - ours))
}
cat(sprintf("logLik %.17g  -opt$objective %.17g\n",
            as.numeric(logLik(fit)), -fit$opt$objective))
