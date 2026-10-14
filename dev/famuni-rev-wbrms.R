# Reviewer probe: wiener(link_ndt = "log") with a covariate on ndt,
# against brms 2.23.0 on the same data. Seed 20261113, n = 300.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("frmtmb from", find.package("frmtmb"), "; eam from",
    find.package("frmtmb.eam"), "\n")
seed <- 20261113L
set.seed(seed)
n <- 300
x <- rnorm(n)
mu <- 0.7 + 0.5 * x
ndt <- exp(log(0.3) + 0.15 * x)
sims <- lapply(seq_len(n), function(i) {
  RWiener::rwiener(1, alpha = 1.2, tau = ndt[i], beta = 0.5, delta = mu[i])
})
dat <- data.frame(x = x, rt = vapply(sims, function(s) s$q, 0),
                  d = vapply(sims, function(s) as.integer(s$resp == "upper"),
                             0L))
fit <- frm(frmtmb::bf(rt | dec(d) ~ x, ndt ~ x),
           family = wiener(link_ndt = "log"), data = dat)
cat("eam ML conv", fit$opt$convergence, "logLik",
    format(as.numeric(logLik(fit)), digits = 10), "\n")
fe <- fixef(fit)
print(fe)
sp <- summary(fit)$spec_pars
print(sp)
# the log density at the ML against RWiener at the same parameters
lp <- function(nm) frm_linpred(fit, dpar = nm, type = "response")
dr <- mapply(function(t, r, a, tau, w, v) {
  RWiener::dwiener(t, alpha = a, tau = tau, beta = w, delta = v,
                   resp = if (r == 1) "upper" else "lower", give_log = TRUE)
}, dat$rt, dat$d, lp("bs"), lp("ndt"), lp("bias"), lp("mu"))
cat(sprintf("sum RWiener log density at the ML %.10f; logLik %.10f; diff %.2e\n",
            sum(dr), as.numeric(logLik(fit)), sum(dr) - as.numeric(logLik(fit))))

suppressMessages(library(brms))
t_b <- system.time(bfit <- brm(brms::bf(rt | dec(d) ~ x, ndt ~ x),
                               family = brms::wiener(), data = dat,
                               chains = 4, iter = 2000, warmup = 1000,
                               seed = seed, backend = "rstan", refresh = 0,
                               silent = 2))
cat("brms fit in", t_b[["elapsed"]], "s\n")
print(prior_summary(bfit))
drw <- as_draws_df(bfit)
cat("max Rhat", max(rhat(bfit), na.rm = TRUE), "; divergences",
    sum(nuts_params(bfit)$Value[nuts_params(bfit)$Parameter == "divergent__"]),
    "\n")
ix <- c(Intercept = "b_Intercept", x = "b_x", ndt_Intercept = "b_ndt_Intercept",
        ndt_x = "b_ndt_x")
tab <- data.frame(par = c(names(ix), "bs", "bias"),
  ml = c(fe[names(ix), "Estimate"], sp["bs", "Estimate"], sp["bias", "Estimate"]),
  se = c(fe[names(ix), "Est.Error"], sp["bs", "Est.Error"], sp["bias", "Est.Error"]),
  med = c(vapply(ix, function(p) median(drw[[p]]), 0), median(drw$bs),
          median(drw$bias)),
  psd = c(vapply(ix, function(p) sd(drw[[p]]), 0), sd(drw$bs), sd(drw$bias)))
tab$z <- (tab$ml - tab$med) / tab$psd
tab$se_sd <- tab$se / tab$psd
print(format(tab, digits = 5), row.names = FALSE)
cat("max |z|", max(abs(tab$z)), " SE/SD range", range(tab$se_sd), "\n")
