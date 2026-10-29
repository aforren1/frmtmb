# Repro for thres() x nl (pair p049, seed 49): an ordinal model whose
# predictor is written as a nonlinear body b1 * x must have the same
# maximum as the linear y ~ x, and brms 2.23.0's maximum (flat priors,
# rstan::optimizing) is the reference.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)
d <- E$eq_dat("ordinal", seed = 49)
f_lin <- frm(bf(y | thres(3) ~ x), family = cumulative(), data = d)
f_nl <- frm(bf(y | thres(3) ~ b1 * x, b1 ~ 1, nl = TRUE),
            family = cumulative(), data = d)
f_nl0 <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(),
             data = d)
cat("logLik linear:", logLik(f_lin), "\n")
cat("logLik nl with thres(3):", logLik(f_nl), "\n")
cat("logLik nl without thres():", logLik(f_nl0), "\n")
print(fixef(f_lin)); print(fixef(f_nl))
f_0 <- frm(bf(y ~ 1), family = cumulative(), data = d)
cat("logLik intercept only:", logLik(f_0), "\n")
print(frm_compat("nl", "cumulative"))
print(names(f_nl$estimates)); print(f_nl$estimates)
for (fam in c("sratio", "cratio", "acat", "cumulative")) {
  fl <- frm(bf(y ~ x), family = get(fam)(), data = d)
  fn <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = get(fam)(), data = d)
  cat(sprintf("%-10s linear %.4f  nl %.4f  nl coefficients: %d\n", fam,
              logLik(fl), logLik(fn), length(unlist(fn$estimates$beta))))
}
