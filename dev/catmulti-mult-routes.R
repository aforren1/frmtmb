# Routes that read a family's mean, on multinomial() now that it has
# one (post$mean_fn): what each does on the lane build and on 0.71.0.
# Arg 1: "base" or "lane". Output: dev/catmulti-out/mult-routes-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- "C:/Users/adf44/source/r/wt-catmulti-lib"
base <- "C:/Users/adf44/source/r/rellib-r9"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(lib, base, user) else c(base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
source("dev/catmulti-data.R")
dm <- catmulti_mult_data(12)
fm <- frm(bf(Yc | trials(N) ~ x), family = multinomial(K = 3), data = dm)
tr <- function(nm, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(nm, ": ")
  if (is.character(r)) cat(substr(r, 1, 300), "\n") else {
    cat(class(r)[1], paste(dim(r) %||% length(r), collapse = " x "), "\n")
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a
tr("residuals", residuals(fm))
tr("residuals pearson", residuals(fm, type = "pearson"))
tr("conditional_effects", conditional_effects(fm))
tr("conditional_effects categorical", conditional_effects(fm,
                                                           categorical = TRUE))
tr("emmeans epred", emmeans::emmeans(fm, ~ x, epred = TRUE))
tr("emmeans link", emmeans::emmeans(fm, ~ x, dpar = "mu2"))
tr("marginaleffects predictions",
   marginaleffects::predictions(fm))
tr("bayes_R2", bayes_R2(fm))
tr("loo", loo(fm))
tr("pp_check", pp_check(fm))
tr("insight get_predicted", insight::get_predicted(fm))
tr("frm_linpred response", frm_linpred(fm, type = "response"))
tr("frm_linpred response se", frm_linpred(fm, type = "response",
                                          se.fit = TRUE))
tr("predict", predict(fm))
tr("dharma", dharma_residuals(fm))
