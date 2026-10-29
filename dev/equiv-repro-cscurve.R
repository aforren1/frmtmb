# Repro for frm_curve() on an ordinal fit with a cs() term (pair p013,
# seed 13). frm_curve() and frm_linpred(type = "link") return the shared
# predictor eta alone, while brms 2.23.0's posterior_linpred(), at the
# same parameter values (a Fixed_param brmsfit), returns eta + cs_k, one
# column per threshold. fitted(scale = "linear") already returns brms's
# layout; its own comment calls the shared part "the predictor at no
# threshold of the model" (R/predict.R, ord_cs_offsets()).
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
spec <- E$eq_specs[["cs_pred() x frm_curve"]]()
fit <- E$eq_frm_fit(spec)
nd <- spec$margs$newdata
cat("frm_curve .estimate:\n")
print(frm_curve(fit, newdata = nd, simultaneous = FALSE)$.estimate)
cat("frm_linpred(type = 'link'):\n")
print(unname(frm_linpred(fit, newdata = nd, type = "link")))
cat("fitted(scale = 'linear') Estimate:\n")
print(fitted(fit, newdata = nd, scale = "linear")[, "Estimate", ])
bfit <- E$eq_fixed(spec, brms::bf(y ~ x + cs(z)), brms::sratio(),
                   spec$data, fit)
cat("brms posterior_linpred at the same values (draw 1):\n")
print(brms::posterior_linpred(bfit, newdata = nd, re_formula = NA)[1, , ])
