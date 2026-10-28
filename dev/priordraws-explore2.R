# What frmtmb's prior tables hold for a brms-example-1-like model: the
# slot table default_prior() builds and the effective prior list
# frm_sample() resolves (rellib-r14).
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
d <- brms::epilepsy
set.seed(1)
d$Exp <- sample(1:5, nrow(d), TRUE)
d$volume <- rnorm(nrow(d))
d$visit <- factor(d$visit)
f <- bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + (1 + Trt | visit),
        sigma ~ Trt)
fit <- frm(f, data = d, family = student())
tab <- default_prior(fit, route = "sample")
print(as.data.frame(tab))
cat("class of tab:", class(tab), "\n")
pr <- set_prior("normal(0,2)", class = "b") +
  set_prior("normal(0,3)", class = "b", dpar = "sigma") +
  set_prior("cauchy(0,2)", class = "sd")
rp <- frmtmb.sample:::sample_resolve_priors(fit, pr)
print(as.data.frame(rp$effective))
str(lapply(rp$ri$entries, function(e) e[c("comp", "idx", "scale")]))
print(frmtmb::brms_coef_table(fit))
print(frmtmb::outer_par_names(fit))
