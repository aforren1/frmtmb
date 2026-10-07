# Lane surface, item 7: what frmtmb does with class "Intercept" on a
# sum-to-zero threshold vector, alone and as a mixture component, beside
# brms's default_prior() rows (dev/surface-brms-stz.R has brms's Stan
# code for the same models).
#
#   Rscript dev/surface-stz.R lane|base > dev/surface-out/stz-<arm>.txt
arm <- commandArgs(TRUE)[1]
source("dev/surface-env.R")
surface_env(arm)
suppressPackageStartupMessages(library(frmtmb))
show <- surface_show
set.seed(3)
dp <- data.frame(x = rnorm(200))
dp$y <- sample(1:5, 200, TRUE)
fam1 <- cumulative(threshold = "sum_to_zero")
famm <- mixture(cumulative(threshold = "sum_to_zero"), cumulative())
show("one family, set_prior(class = 'Intercept')",
     frm(y ~ x, data = dp, family = fam1,
         prior = set_prior("normal(0, 3)", class = "Intercept")))
show("mixture, set_prior(class = 'Intercept', dpar = 'mu1')",
     frm(y ~ x, data = dp, family = famm,
         prior = set_prior("normal(0, 3)", class = "Intercept",
                           dpar = "mu1")))
show("mixture, set_prior(class = 'Intercept', dpar = 'mu2') (control)",
     frm(y ~ x, data = dp, family = famm,
         prior = set_prior("normal(0, 3)", class = "Intercept",
                           dpar = "mu2")))
f <- frm(y ~ x, data = dp, family = famm)
cat("thresholds of the stz component, sum:",
    format(sum(fixef(f)[grep("^mu1_Intercept", rownames(fixef(f))),
                        "Estimate"])), "\n")
cat("free threshold parameters of the stz component:",
    length(f$estimates[["tau_raw_mu1"]] %||% NA), "\n")
print(names(f$estimates))
