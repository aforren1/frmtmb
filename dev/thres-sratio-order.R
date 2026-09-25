# Lane thres: frmtmb's sratio() orders its thresholds and brms 2.23.0's
# does not (brms:::has_ordered_thres() is FALSE for sratio: its Stan
# program declares `vector[nthres] Intercept`, no `ordered`). At
# frmtmb's grouped sratio estimate (data seed 11, as in
# dev/thres-brms-stan.R) level a's second and third thresholds sit on
# frmtmb's ordering boundary, and brms's gradient there is +-0.544 on
# exactly those two: brms's likelihood rises when they cross.
#
#   Rscript dev/thres-sratio-order.R > dev/thres-sratio-order-log.txt 2>&1
.libPaths(c("/opt/rlib/stan", "/opt/rlib/lane-thres", "/opt/rlib/base",
            "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages({library(frmtmb); library(brms); library(rstan)})
set.seed(11)
n <- 240
d <- data.frame(x = rnorm(n), g = sample(c("a", "b", "c"), n, TRUE))
u <- rlogis(n, 0.7 * d$x)
tau <- list(a = c(-1, 0.2, 1.3), b = c(-0.5, 0.6), c = c(-1.5, -0.4, 0.5, 1.6))
d$y <- vapply(seq_len(n), function(i) 1L + sum(u[i] > tau[[d$g[i]]]), 1L)
fit <- frm(y | thres(gr = g) ~ x, data = d, family = sratio())
bf0 <- bf(y | thres(gr = g) ~ x)
pr <- get_prior(bf0, data = d, family = brms::sratio()); pr$prior <- ""
code <- make_stancode(bf0, data = d, family = brms::sratio(), prior = pr)
sdat <- make_standata(bf0, data = d, family = brms::sratio(), prior = pr)
sf <- suppressMessages(sampling(stan_model(model_code = code),
                                data = sdat, chains = 0))
th <- fit$spec$responses$y$family$post$ord_thresholds(fit$estimates$tau_raw)
lay <- frmtmb:::thres_layout(c(3, 2, 4))
pars <- list(b = array(fit$estimates$beta, 1))
for (k in 1:3) {
  pars[[paste0("Intercept_", k)]] <- as.array(th[lay$start[k]:lay$end[k]])
}
up <- unconstrain_pars(sf, pars)
print(up); print(fit$opt$par)
print(grad_log_prob(sf, up, adjust_transform = FALSE))
print(fit$obj$gr(fit$opt$par))
# move level a's third threshold 0.05 BELOW its second, which brms's
# unordered vector allows and frmtmb's ordered one does not
up2 <- up
up2[4] <- up[4] - 0.05
cat("brms log_prob gain from crossing:",
    log_prob(sf, up2, adjust_transform = FALSE) -
      log_prob(sf, up, adjust_transform = FALSE), "\n")
