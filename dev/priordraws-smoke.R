# First run of sample_prior, prior_draws(), inc_warmup and update() on
# the lane's build.
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)),
                m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$g] + rnorm(60)
q <- function(lab, expr) {
  cat("\n=====", lab, "=====\n")
  withCallingHandlers(
    tryCatch(expr, error = function(e) cat("ERROR:", conditionMessage(e), "\n")),
    warning = function(w) {
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
}
pr <- set_prior("normal(0, 1)", class = "b")
ds <- q("yes", frm_sample(bf(y ~ x + z + (1 + x | g)), data = d,
                          family = gaussian(), prior = pr,
                          sample_prior = "yes", chains = 2, iter = 600,
                          seed = 1, refresh = 0))
q("names", print(names(prior_draws(ds))))
q("tail cols", print(tail(colnames(ds$draws), 8)))
q("dims", print(c(nrow(ds$draws), nrow(ds$warmup), ndraws(ds))))
q("prior means", print(round(colMeans(prior_draws(ds)), 3)))
q("prior_draws b_x", print(head(prior_draws(ds, "b_x"), 3)))
q("prior_draws mixed", print(names(prior_draws(ds, c("b_x", "sigma",
                                                     "sd_g__x", "zzz")))))
q("prior_draws b_Intercept", print(dim(prior_draws(ds, "b_Intercept"))))
q("hyp", print(hypothesis(ds, c("x = 0", "x > 0", "x - z = 0",
                                "Intercept = 0"))$hypothesis))
q("hyp sigma", print(hypothesis(ds, "sigma = 1", class = NULL)$hypothesis))
q("summary rows", print(rownames(summary(ds))))
q("inc_warmup df", print(nrow(as_draws_df(ds, inc_warmup = TRUE))))
q("inc_warmup array", print(dim(as_draws_array(ds, inc_warmup = TRUE))))
q("as.mcmc warmup", print(dim(as.mcmc(ds, inc_warmup = TRUE)[[1]])))
q("as.mcmc combine warmup", as.mcmc(ds, inc_warmup = TRUE,
                                    combine_chains = TRUE))
q("rhat warmup", print(head(rhat(ds, inc_warmup = TRUE), 3)))
q("rhat", print(head(rhat(ds), 3)))
q("nuts", print(dim(nuts_params(ds))))
q("nuts warm", print(dim(nuts_params(ds, inc_warmup = TRUE))))
q("only flat", frm_sample(bf(y ~ x + (1 | g)), data = d, family = gaussian(),
                          sample_prior = "only", chains = 1, iter = 200,
                          refresh = 0))
do <- q("only", frm_sample(bf(y ~ x + z + (1 + x | g)), data = d,
                           family = gaussian(), prior = pr,
                           sample_prior = "only", chains = 2, iter = 1000,
                           seed = 2, refresh = 0))
q("only summary", print(summary(do)[, 1:4]))
q("only prior_draws", print(dim(prior_draws(do))))
q("only pp", print(dim(posterior_predict(do))))
d0 <- q("no", frm_sample(bf(y ~ x + z + (1 + x | g)), data = d,
                         family = gaussian(), prior = pr, chains = 2,
                         iter = 600, seed = 1, refresh = 0))
q("no vs yes same post", print(identical(d0$draws,
                                         ds$draws[, colnames(d0$draws)])))
q("no hyp", print(hypothesis(d0, "x = 0")$hypothesis))
du <- q("update", update(d0, sample_prior = "yes"))
q("update same", print(identical(du$draws[, colnames(d0$draws)],
                                 ds$draws[, colnames(d0$draws)])))
q("update prior", print(identical(prior_draws(du), prior_draws(ds))))
