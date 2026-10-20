# Lane sampleprior, fix round 2: why rdm() under ndt_group() still
# fails to start with the closed guards (data seed 20261121): the ndt
# coordinate at the ML mode and the sampling objective there.
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
  library(frmtmb.eam)
})
set.seed(20261121)
d <- rdm_simulate(200, v = c(2.4, 1.4), A = 0.5, k = 0.5, ndt = 0.2)
d$subj <- factor(rep(1:4, length.out = nrow(d)))
fit <- suppressWarnings(frm(bf(rt | vint(choice) + ndt_group(subj) ~ 1),
                            family = rdm(2), data = d))
p <- fit$obj$env$last.par.best
print(p)
lp <- Filter(function(l) identical(l$dpar, "ndt"), fit$frame$linpreds)[[1]]
cat("ndt link", lp$link$name, "\n")
j <- grep("ndt", names(fit$frame$par_template$betad))
eta <- fit$estimates$betad[j]
cat("eta", eta, " linkinv", sprintf("%.17g", lp$link$linkinv(eta)),
    " mu_eta", lp$link$mu_eta(eta), "\n")
rp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))
fit2 <- fit
fit2$prior <- rp$effective
obj <- frmtmb.sample:::prior_augmented_obj(fit2, rp$ri$entries)
cat("likelihood nll at mode", fit$obj$fn(p), "\n")
cat("sampling nll at mode", obj$fn(p), "\n")
nlp <- frmtmb::neg_log_prior_fn(rp$ri$entries)
for (e in rp$ri$entries) {
  one <- frmtmb::neg_log_prior_fn(list(e))
  cat(sprintf("  entry %s[%s] %s: %s\n", e$comp,
              paste(e$idx, collapse = ","), e$dist$kind,
              format(one(fit$estimates))))
}
