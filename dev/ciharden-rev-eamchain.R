# Reviewer: how well does test-sampling.R's one wiener chain mix, and
# does it sit at the mode? The test allows its ESS and R-hat warnings.
.libPaths(c("C:/Users/adf44/source/r/wt-ciharden-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(frmtmb.eam))
ns <- asNamespace("frmtmb.eam")
set.seed(77)
dat <- ns$ddm_simulate(250, mu = 0.9, bs = 1.4, ndt = 0.25)
fit <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5), family = wiener(),
           data = dat)
w <- character()
smp <- withCallingHandlers(
  frmtmb.sample::frm_sample(fit, chains = 1, iter = 400, warmup = 200,
                            seed = 3, refresh = 0),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
cat("R", R.home(), "threads", Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
cat("warnings:", length(w), "\n")
for (x in w) cat("  -", substr(gsub("\n", " ", x), 1, 160), "\n")
sf <- smp$stanfit
s <- rstan::summary(sf)$summary
print(round(s[, c("mean", "sd", "n_eff", "Rhat")], 4))
sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)[[1]]
cat("divergent:", sum(sp[, "divergent__"]), " max treedepth hits:",
    sum(sp[, "treedepth__"] >= 10), "\n")
# the mode on the sampler's scale, against the posterior mean and sd
est <- fit$opt$par
pm <- s[seq_along(est), "mean"]
psd <- s[seq_along(est), "sd"]
cat("z = (posterior mean - mode) / posterior sd:\n")
print(round((pm - est) / psd, 3))
