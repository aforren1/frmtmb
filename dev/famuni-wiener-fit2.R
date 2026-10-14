# Fit agreement: frmtmb.eam's wiener (ML) against brms 2.23.0 (rstan,
# 4 chains x 1000 post-warmup draws) on one simulated data set with a
# covariate on mu and on bs. Also the post-fit semantics: log_lik(),
# posterior_predict() and posterior_epred().
# Seed 20261013 for the data, 20261013 for brms.
# Output: dev/famuni-out/wiener-fit2.txt; the brmsfit is cached in
# dev/famuni-out/wiener-brmsfit.rds.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam",
  quiet = TRUE))
suppressMessages(library(brms))
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out"
seed <- 20261013L

## data: RWiener draws, one per row
set.seed(seed)
n <- 400
x <- rnorm(n)
truth <- c(b_Intercept = 0.6, b_x = 0.8, b_bs_Intercept = 0.25,
           b_bs_x = 0.2, ndt = 0.3, bias = 0.5)
mu <- truth[["b_Intercept"]] + truth[["b_x"]] * x
bs <- exp(truth[["b_bs_Intercept"]] + truth[["b_bs_x"]] * x)
sims <- lapply(seq_len(n), function(i) {
  RWiener::rwiener(1, alpha = bs[i], tau = truth[["ndt"]],
                   beta = truth[["bias"]], delta = mu[i])
})
dat <- data.frame(x = x, rt = vapply(sims, function(s) s$q, 0),
                  d = vapply(sims, function(s) as.integer(s$resp == "upper"),
                             0L))
cat("n", n, "min rt", min(dat$rt), "share upper", mean(dat$d), "\n")

## frmtmb.eam, ML
form <- frmtmb::bf(rt | dec(d) ~ x, bs ~ x)
bform <- brms::bf(rt | dec(d) ~ x, bs ~ x)
t_ml <- system.time(fit <- frm(form, family = frmtmb.eam::wiener(), data = dat))
cat("\nfrmtmb.eam ML in", t_ml[["elapsed"]], "s; convergence",
    fit$opt$convergence, "logLik", format(as.numeric(logLik(fit)),
                                         digits = 10), "\n")
print(summary(fit))
fe <- frmtmb::fixef(fit)
print(fe)
ub <- family(fit)$ndt_bound$ub
cat("ndt bound", ub, "\n")

## brms, cached
bpath <- file.path(out_dir, "wiener-brmsfit.rds")
if (file.exists(bpath)) {
  bfit <- readRDS(bpath)
  cat("\nbrms fit read from cache\n")
} else {
  t_b <- system.time(bfit <- brm(bform, family = brms::wiener(), data = dat,
                                 chains = 4, iter = 2000, warmup = 1000,
                                 seed = seed, backend = "rstan",
                                 refresh = 0, silent = 2))
  cat("\nbrms fit in", t_b[["elapsed"]], "s\n")
  saveRDS(bfit, bpath)
}
print(prior_summary(bfit))
print(summary(bfit))
dr <- as_draws_df(bfit)
cat("max Rhat", max(rhat(bfit), na.rm = TRUE), " min bulk ESS",
    min(posterior::summarise_draws(dr)$ess_bulk, na.rm = TRUE), "\n")
cat("divergences", sum(nuts_params(bfit)$Value[
  nuts_params(bfit)$Parameter == "divergent__"]), "\n")

## the comparison. The scalar ndt and bias are reported by summary() on
## their natural scale, as brms reports them, with a delta-method SE.
bd <- list(
  mu_Intercept = dr$b_Intercept, mu_x = dr$b_x,
  bs_Intercept = dr$b_bs_Intercept, bs_x = dr$b_bs_x,
  ndt = dr$ndt, bias = dr$bias)
ix <- c(mu_Intercept = "Intercept", mu_x = "x",
        bs_Intercept = "bs_Intercept", bs_x = "bs_x")
stopifnot(all(ix %in% rownames(fe)))
sp <- summary(fit)$spec_pars
print(sp)
est <- c(setNames(fe[ix, "Estimate"], names(ix)),
         ndt = sp["ndt", "Estimate"], bias = sp["bias", "Estimate"])
ses <- c(setNames(fe[ix, "Est.Error"], names(ix)),
         ndt = sp["ndt", "Est.Error"], bias = sp["bias", "Est.Error"])
tab <- data.frame(
  par = names(est), truth = truth,
  ml = est, se = ses,
  post_median = vapply(bd[names(est)], stats::median, 0),
  post_sd = vapply(bd[names(est)], stats::sd, 0))
tab$z_shift <- (tab$ml - tab$post_median) / tab$post_sd
tab$se_over_sd <- tab$se / tab$post_sd
cat("\n-- ML against posterior median, in posterior SD\n")
print(tab, digits = 4, row.names = FALSE)
cat("max |z_shift|", max(abs(tab$z_shift)), "\n")

## log_lik(): both are the pointwise log density at the observed boundary.
## Evaluate eam's density at brms's draws and compare to brms's log_lik.
set.seed(seed)
# sorted: brms returns the draws of draw_ids in ascending order
ids <- sort(sample(nrow(dr), 50))
ll_b <- brms::log_lik(bfit, draw_ids = ids)
mu_d <- outer(dr$b_Intercept[ids], rep(1, n)) + outer(dr$b_x[ids], dat$x)
bs_d <- exp(outer(dr$b_bs_Intercept[ids], rep(1, n)) +
              outer(dr$b_bs_x[ids], dat$x))
ll_e <- matrix(NA_real_, length(ids), n)
for (k in seq_along(ids)) {
  ll_e[k, ] <- ddm_lpdf_both(dat$rt - dr$ndt[ids[k]], mu_d[k, ], bs_d[k, ],
                             dr$bias[ids[k]], dat$d)
}
cat("\n-- log_lik: eam density at 50 brms draws against brms log_lik\n")
cat(sprintf("max abs diff %.3e, max rel %.3e\n", max(abs(ll_e - ll_b)),
            max(abs(ll_e - ll_b) / pmax(1, abs(ll_b)))))
# frmtmb's log_lik() on an ML fit refuses by design (it is per draw), so
# the ML comparison is the family density at the fitted dpars, summed
# against logLik(), and against brms's R dwiener (its log_lik backend)
r <- tryCatch(frmtmb::log_lik(fit), error = function(e) conditionMessage(e))
cat("frmtmb::log_lik(ML fit):", substr(paste(r), 1, 90), "\n")
dp <- lapply(c("mu", "bs", "ndt", "bias"), function(p) {
  frm_linpred(fit, dpar = p, type = "response")
})
names(dp) <- c("mu", "bs", "ndt", "bias")
ll_ml <- ddm_lpdf_both(dat$rt - dp$ndt, dp$mu, dp$bs, dp$bias, dat$d)
cat(sprintf("sum of the density at the ML %.10f, logLik %.10f\n",
            sum(ll_ml), as.numeric(logLik(fit))))
ll_bml <- brms::dwiener(dat$rt, alpha = dp$bs, tau = dp$ndt, beta = dp$bias,
                        delta = dp$mu, resp = dat$d, log = TRUE)
cat(sprintf("eam density against brms::dwiener at the ML: max abs %.3e\n",
            max(abs(ll_ml - ll_bml))))

## posterior_predict(): brms draws the boundary as well and returns the
## time only, or the time signed by the DRAWN boundary.
set.seed(seed)
pp <- brms::posterior_predict(bfit, ndraws = 400)
set.seed(seed)
ppn <- brms::posterior_predict(bfit, ndraws = 400, negative_rt = TRUE)
cat("\n-- posterior_predict\n")
cat("negative values, default:", sum(pp < 0), " negative_rt = TRUE share",
    format(mean(ppn < 0), digits = 4), " observed share lower",
    format(mean(dat$d == 0), digits = 4), "\n")
neg_by_row <- colMeans(ppn < 0)
cat("share of lower-boundary draws on rows observed UPPER",
    format(mean(neg_by_row[dat$d == 1]), digits = 4), "; on rows observed",
    "LOWER", format(mean(neg_by_row[dat$d == 0]), digits = 4), "\n")
set.seed(seed)
sm <- simulate(fit, nsim = 400)
smm <- as.matrix(sm)
cat("eam simulate() dims", paste(dim(smm), collapse = " x "),
    " any negative", any(smm < 0), "\n")
# eam conditions on the observed boundary, so its mean draw per row is
# the conditional mean; brms's is the marginal mean
up_rows <- dat$d == 1
cat("mean simulated rt, rows observed upper: eam",
    format(mean(smm[up_rows, ]), digits = 5), " brms",
    format(mean(pp[, up_rows]), digits = 5), "\n")
cat("mean simulated rt, rows observed lower: eam",
    format(mean(smm[!up_rows, ]), digits = 5), " brms",
    format(mean(pp[, !up_rows]), digits = 5), "\n")

## posterior_epred(): brms's formula, called as brms calls it
cat("\n-- posterior_epred\n")
ep_b <- brms::posterior_epred(bfit, ndraws = 400)
cat("brms epred dims", paste(dim(ep_b), collapse = " x "),
    " rows observed upper", format(mean(ep_b[, up_rows]), digits = 5),
    " lower", format(mean(ep_b[, !up_rows]), digits = 5), "\n")
f_eam <- fitted(fit)
cat("eam fitted() rows observed upper",
    format(mean(f_eam[up_rows, 1]), digits = 5), " lower",
    format(mean(f_eam[!up_rows, 1]), digits = 5), "\n")
epf <- get("posterior_epred_wiener", asNamespace("brms"))
chk <- expand.grid(mu = c(-1.5, -0.5, 0.5, 1.5), bs = c(0.8, 1, 1.6, 2.5),
                   bias = c(0.3, 0.5, 0.7), ndt = 0.3)
chk$brms <- epf(list(dpars = as.list(chk[c("mu", "bs", "ndt", "bias")])))
pu <- ddm_p_upper(chk$mu, chk$bs, chk$bias)
chk$exact <- chk$ndt +
  pu * ddm_cond_mean_dt(chk$mu, chk$bs, chk$bias, rep(1, nrow(chk))) +
  (1 - pu) * ddm_cond_mean_dt(chk$mu, chk$bs, chk$bias, rep(0, nrow(chk)))
# the same Grasman, Wagenmakers and van der Maas (2009) formula with the
# ABSOLUTE start point z = bias * bs, which is what it takes
z <- chk$bias * chk$bs
chk$brms_z <- chk$ndt - z / chk$mu + chk$bs / chk$mu *
  (exp(-2 * chk$mu * z) - 1) / (exp(-2 * chk$mu * chk$bs) - 1)
# Monte Carlo from RWiener at four of the cells
set.seed(seed)
mc_rows <- c(which(chk$bs == 1)[1], which(chk$bs == 2.5)[c(1, 6, 12)])
chk$mc <- NA_real_
for (r in mc_rows) {
  s <- RWiener::rwiener(20000, chk$bs[r], chk$ndt[r], chk$bias[r], chk$mu[r])
  chk$mc[r] <- mean(s$q)
}
chk$rel_err_brms <- (chk$brms - chk$exact) / chk$exact
cat("brms epred formula against the exact unconditional mean\n")
print(chk, digits = 5, row.names = FALSE)
cat("max |rel err| brms formula, bs == 1:",
    format(max(abs(chk$rel_err_brms[chk$bs == 1])), digits = 3),
    "; bs != 1:", format(max(abs(chk$rel_err_brms[chk$bs != 1])), digits = 3),
    "; with z = bias * bs:",
    format(max(abs((chk$brms_z - chk$exact) / chk$exact)), digits = 3), "\n")
cat("brms epred at mu = 0:", epf(list(dpars = list(mu = 0, bs = 1.4,
                                                   ndt = 0.3, bias = 0.5))),
    "\n")
