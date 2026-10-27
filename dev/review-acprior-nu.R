# Reviewer: the bounded-default analogue the lane left out. brms keeps
# nu's declared lower bound 1 under a user prior (real<lower=1> nu;
# dev/review-acprior-log/brms.txt). Does frm_sample()?
# Data: 120 draws of a t with 1.2 degrees of freedom, seed 12; prior
# normal(2, 2) on nu; 4 chains of 2000, seed 3. brms's identity nu link.
# Log: dev/review-acprior-log/nu.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
.libPaths(c("C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(12)
d <- data.frame(x = rnorm(120))
d$y <- 1 + d$x + rt(120, df = 1.2)
fit <- frm(y ~ x, data = d, family = brmsfamily("student", link_nu = "identity"))
cat("nu link:", fit$spec$family$link_nu %||% "?", "
")
pr <- set_prior("normal(2, 2)", class = "nu")
r0 <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))$ri
r1 <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pr))$ri
cat("default lower:", paste(names(r0$lower), if (length(r0$lower)) signif(unlist(r0$lower), 6)),
    "\nuser-prior lower:", if (length(r1$lower)) paste(names(r1$lower),
    signif(unlist(r1$lower), 6)) else "(none)", "\n")
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, prior = pr, chains = 4, iter = 2000, seed = 3,
             refresh = 0)))
nu <- as.matrix(ds)[, "nu"]
cat(sprintf("nu: mean %.3f q05 %.3f share below 1 %.4f (brms: 0 by declaration)\n",
            mean(nu), quantile(nu, 0.05), mean(nu < 1)))
