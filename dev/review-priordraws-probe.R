# Reviewer probes on the lane build: leaks, warmup readers, thinning,
# evidence ratio against brms's own functions on the same draws,
# update(), memory.
#   Rscript dev/review-priordraws-probe.R
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-stancache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
options(mc.cores = 1)
q <- function(expr) suppressWarnings(suppressMessages(expr))
tryv <- function(expr) tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
set.seed(1)
d <- data.frame(x = rnorm(80), z = rnorm(80), g = factor(rep(1:8, 10)))
d$y <- 1 + 0.5 * d$x + 0.2 * d$z + rnorm(8, 0, 0.7)[d$g] + rnorm(80)
pr <- set_prior("normal(0, 1)", class = "b")
f <- bf(y ~ x + z + (1 | g))
run <- function(...) q(frm_sample(f, data = d, family = gaussian(), prior = pr,
                                  chains = 2, iter = 1000, seed = 21,
                                  refresh = 0, ...))

cat("\n#### P1 leaks of prior_only\n")
only <- run(sample_prior = "only")
cat("only$fit$frame$prior_only:", deparse(only$fit$frame[["prior_only"]]), "\n")
ll <- q(log_lik(only))
cat("log_lik(only) dim", dim(ll), "finite", all(is.finite(ll)), "\n")
# log_lik of one draw against an independent normal density
dr1 <- as.data.frame(only)[1, ]
mu1 <- dr1[["b_Intercept"]] + dr1[["b_x"]] * d$x + dr1[["b_z"]] * d$z +
  unlist(dr1[paste0("r_g[", 1:8, ",Intercept]")])[as.integer(d$g)]
ref <- dnorm(d$y, mu1, dr1[["sigma"]], log = TRUE)
cat("log_lik(only) draw 1 vs reference, max rel diff",
    max(abs(ll[1, ] - ref) / abs(ref)), "\n")
tf <- tempfile(fileext = ".rds")
saveRDS(only, tf)
only2 <- readRDS(tf)
cat("readRDS prior_only:", deparse(only2$fit$frame[["prior_only"]]), "\n")
fm <- q(frm(f, data = d, family = gaussian()))
cat("frm after only: logLik", format(as.numeric(logLik(fm)), digits = 10),
    " frame prior_only", deparse(fm$frame[["prior_only"]]), "\n")
nobj <- frmtmb::build_objective(fm$frame)
cat("build_objective(fit frame) at par:", nobj(fm$obj$env$parList()), "\n")

cat("\n#### P2 prior columns and the methods\n")
yes <- run(sample_prior = "yes")
no <- run()
cat("variables(yes) tail:", tail(variables(yes), 5), "\n")
cat("identical summary:", identical(summary(yes), summary(no)), "\n")
cat("identical fixef:", identical(fixef(yes), fixef(no)), "\n")
cat("identical ranef:", identical(ranef(yes), ranef(no)), "\n")
ps <- posterior_summary(yes)
cat("posterior_summary rows with prior_:", sum(grepl("^prior_", rownames(ps))), "\n")
set.seed(4); p1 <- posterior_predict(yes)
set.seed(4); p2 <- posterior_predict(no)
cat("identical posterior_predict:", identical(p1, p2), "\n")
cat("identical log_lik:", identical(q(log_lik(yes)), q(log_lik(no))), "\n")
l1 <- q(loo(yes)); l2 <- q(loo(no))
cat("identical loo elpd:", identical(l1$estimates, l2$estimates), "\n")
set.seed(4); b1 <- q(bayes_R2(yes)); set.seed(4); b2 <- q(bayes_R2(no))
cat("identical bayes_R2:", identical(b1, b2), "\n")
ce1 <- q(conditional_effects(yes, "x")); ce2 <- q(conditional_effects(no, "x"))
cat("identical conditional_effects:", identical(ce1[[1]]$estimate__, ce2[[1]]$estimate__), "\n")
mp <- tryv(q(mcmc_plot(yes)))
cat("mcmc_plot default variables:", if (is.character(mp)) mp else
  paste(unique(as.character(mp$data$parameter)), collapse = " "), "\n")
pl <- tryv(q(plot(yes, plot = FALSE)))
cat("plot() default pages:", if (is.character(pl)) pl else length(pl), "\n")
pv <- tryv(q(pairs(yes)))
cat("pairs default ok:", !is.character(pv), "\n")
rh <- rhat(yes)
cat("rhat names with prior_:", sum(grepl("^prior_", names(rh))), "\n")
print(tryv(q(prior_summary(yes))))

cat("\n#### P3 warmup never read by other methods\n")
bare <- yes
bare$warmup <- NULL
chk <- function(nm, f) {
  set.seed(7); a <- tryv(q(f(yes)))
  set.seed(7); b <- tryv(q(f(bare)))
  cat(sprintf("  %-22s identical %s\n", nm, identical(a, b)))
}
chk("ranef", ranef)
chk("log_lik", log_lik)
chk("posterior_predict", posterior_predict)
chk("posterior_epred", posterior_epred)
chk("loo", function(x) loo(x)$estimates)
chk("bayes_R2", bayes_R2)
chk("pp_check data", function(x) pp_check(x, ndraws = 5)$data)
chk("conditional_effects", function(x) conditional_effects(x, "x")[[1]])
chk("hypothesis", function(x) hypothesis(x, "x = 0", seed = 1)$hypothesis)
chk("nuts_params", nuts_params)
chk("log_posterior", log_posterior)
chk("ndraws", ndraws)
chk("variables", variables)
chk("print", function(x) capture.output(print(x)))
chk("as_draws_df", as_draws_df)
chk("fitted", fitted)
chk("predict", predict)
chk("coef", coef)
chk("waic", function(x) waic(x)$estimates)
chk("kfold-free: loo_R2", loo_R2)

cat("\n#### P4 thinning and inc_warmup\n")
th <- run(thin = 3, warmup = 400)
sf <- th$stanfit
cat("rstan dims: inc_warmup", dim(rstan::extract(sf, permuted = FALSE, inc_warmup = TRUE))[1:2],
    " post", dim(rstan::extract(sf, permuted = FALSE))[1:2], " warmup2", sf@sim$warmup2, "\n")
cat("as_draws_array(inc_warmup) dims", dim(as_draws_array(th, inc_warmup = TRUE))[1:2], "\n")
cat("as.mcmc(inc_warmup, thin 3):", tryv(dim(as.mcmc(th, inc_warmup = TRUE)[[1]])), "\n")
cat("nsamples(incl_warmup):", tryv(q(nsamples(th, incl_warmup = TRUE))), "\n")

cat("\n#### P5 evidence ratio against brms's own functions, same draws\n")
er_brms <- function(post, prior, sign) {
  brms:::evidence_ratio(post, cut = 0, wsign = sign, prior_samples = prior)
}
dr <- as.data.frame(yes)
hs <- c("x = 0", "x > 0", "x < 0.3", "x - z = 0", "x + z = 0.5")
h <- hypothesis(yes, hs, seed = 3)
set.seed(3)
pd <- list(prior_draws(yes, "b_x")$b_x)
pd[[4]] <- with(prior_draws(yes, c("b_x", "b_z")), b_x - b_z)
pd[[5]] <- with(prior_draws(yes, c("b_x", "b_z")), b_x + b_z - 0.5)
post <- list(dr$b_x, dr$b_x, dr$b_x - 0.3, dr$b_x - dr$b_z, dr$b_x + dr$b_z - 0.5)
sg <- c("equal", "greater", "less", "equal", "equal")
for (k in seq_along(hs)) {
  bv <- er_brms(post[[k]], if (sg[k] == "equal") pd[[k]], sg[k])
  cat(sprintf("  %-12s frmtmb %.10g  brms-fn %.10g\n", hs[k],
              h$hypothesis$Evid.Ratio[k], bv))
}
# the prior draws of x - z: two independent columns of prior_b, sd sqrt(2)
pk <- h$prior_samples[[4]]
cat("prior of x - z: sd", sd(pk), " (normal(0,1) difference: 1.414)\n")
hs2 <- c("sd_g__Intercept = 0.5", "sigma = 1")
h2 <- hypothesis(yes, hs2, class = "", seed = 3)
set.seed(3)
p_sd <- prior_draws(yes, "sd_g__Intercept")[[1]] - 0.5
p_sg <- prior_draws(yes, "sigma")[[1]] - 1
cat(sprintf("  sd: frmtmb %.10g brms-fn %.10g\n", h2$hypothesis$Evid.Ratio[1],
            er_brms(dr$sd_g__Intercept - 0.5, p_sd, "equal")))
cat(sprintf("  sigma: frmtmb %.10g brms-fn %.10g\n", h2$hypothesis$Evid.Ratio[2],
            er_brms(dr$sigma - 1, p_sg, "equal")))
wi <- NULL
hi <- withCallingHandlers(hypothesis(yes, "Intercept = 0"),
                          warning = function(w) {
                            wi <<- c(wi, conditionMessage(w)); invokeRestart("muffleWarning")})
cat("Intercept = 0: ER", hi$hypothesis$Evid.Ratio, " warnings:", wi, "\n")
hn <- hypothesis(no, "x = 0")
cat("no prior draws: ER", hn$hypothesis$Evid.Ratio, "\n")
hp <- tryv(hypothesis(yes, "x = 0")); grDevices::pdf(NULL)
pl2 <- tryv(plot(hp, plot = TRUE)); grDevices::dev.off()
cat("plot hypothesis ok:", !is.character(pl2), "\n")

cat("\n#### P6 update()\n")
u1 <- q(update(no, sample_prior = "yes"))
cat("update(no, yes) == yes draws:", identical(u1$draws, yes$draws), "\n")
u2 <- q(update(yes, seed = 22))
cat("update(seed=22) differs:", !identical(u2$draws, yes$draws),
    " ncol same:", ncol(u2$draws) == ncol(yes$draws), "\n")
u3 <- q(update(yes, chains = 1, iter = 400))
cat("update(chains=1, iter=400): ndraws", ndraws(u3), " nchains", nchains(u3), "\n")
u4 <- q(update(yes, sample_prior = "only"))
cat("update(yes, only): sample_prior", u4$sample_prior, " prior cols",
    sum(grepl("^prior_", colnames(u4$draws))), "\n")
u5 <- q(update(only, sample_prior = "no"))
cat("update(only, no): sd(b_x)", sd(as.data.frame(u5)$b_x), " vs only",
    sd(as.data.frame(only)$b_x), "\n")
cat("update(formula):", tryv(update(no, y ~ x)), "\n")
cat("update(. ~ . + 0):", tryv(update(no, formula. = . ~ . - z)), "\n")
cat("update(newdata):", tryv(update(no, newdata = d)), "\n")
cat("update(unnamed):", tryv(update(no, 5)), "\n")
cat("update(prior = flat):", tryv({u6 <- q(update(no, prior = "flat")); "OK"}), "\n")
fitroute <- q(frm(f, data = d, family = gaussian()))
dsf <- q(frm_sample(fitroute, prior = pr, chains = 2, iter = 1000, seed = 21, refresh = 0))
uf <- q(update(dsf, sample_prior = "yes"))
dsfy <- q(frm_sample(fitroute, prior = pr, chains = 2, iter = 1000, seed = 21,
                     refresh = 0, sample_prior = "yes"))
cat("fit route update(yes) == direct:", identical(uf$draws, dsfy$draws), "\n")

cat("\n#### P7 memory\n")
big <- data.frame(x = rnorm(2000), g = factor(rep(1:400, 5)))
big$y <- 1 + 0.5 * big$x + rnorm(400)[big$g] + rnorm(2000)
for (sw in c(TRUE, FALSE)) {
  t0 <- Sys.time()
  bd <- q(frm_sample(bf(y ~ x + (1 | g)), data = big, family = gaussian(),
                     prior = pr, chains = 2, iter = 1000, seed = 2,
                     refresh = 0, save_warmup = sw))
  cat(sprintf("save_warmup=%s: draws %.1f MB, warmup %.1f MB, stanfit %.1f MB, total %.1f MB, %.1f s\n",
              sw, object.size(bd$draws) / 2^20, object.size(bd$warmup) / 2^20,
              object.size(bd$stanfit) / 2^20, object.size(bd) / 2^20,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
