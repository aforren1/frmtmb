# Reviewer: other post-fit generics on simplex fits, to find internal
# errors (not a frm_stop() sentence). Output: dev/famsimplex-rev-generics.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
show <- function(lab, expr) {
  r <- tryCatch(suppressMessages(suppressWarnings({
    v <- expr
    paste("OK:", paste(class(v), collapse = "/"),
          paste(dim(v) %||% length(v), collapse = "x"))
  })), error = function(e) {
    cl <- conditionCall(e)
    paste0("ERROR", if (inherits(e, "frmtmb_error") || is.null(cl)) "" else
      paste0(" [in ", deparse1(cl)[1], "]"), ": ",
      substr(conditionMessage(e), 1, 160))
  })
  cat(sprintf("%-28s %s\n", lab, r))
}
dd <- fs_data_dirichlet(seed = 12, n = 300, ng = 15, phi_formula = FALSE)
dl <- fs_data_logistic_normal(seed = 14, n = 300, ng = 15,
                              sigma_formula = FALSE)
fits <- list(
  dir = frm(bf(Y ~ x + (1 | g)), family = dirichlet(), data = dd),
  ln = frm(bf(Y ~ x + (1 | g)), family = logistic_normal(), data = dl))
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat("\n==", nm, "\n")
  show("ranef", ranef(f))
  show("coef", coef(f))
  show("VarCorr", VarCorr(f))
  show("confint_varcorr", confint_varcorr(f))
  show("predict type link", predict(f, type = "link"))
  show("fitted re_formula NA", fitted(f, re_formula = NA))
  show("fitted newlevel", fitted(f, newdata = data.frame(x = 0, g = "new"),
                                 allow_new_levels = TRUE))
  show("predict newlevel", predict(f, newdata = data.frame(x = 0, g = "new"),
                                   allow_new_levels = TRUE))
  show("frm_bootstrap nsim 3", frm_bootstrap(f, nsim = 3, seed = 1))
  show("confint boot", confint(f, parm = "mub_x", method = "boot", nsim = 3,
                               seed = 1))
  show("profile confint", confint(f, parm = "mub_x", method = "profile"))
  show("dharma_residuals", dharma_residuals(f, nsim = 20, seed = 1))
  show("cooks.distance", cooks.distance(f))
  show("influence", influence(f))
  show("model.frame", model.frame(f))
  show("nobs", nobs(f))
  show("AIC", AIC(f))
  show("vcov", vcov(f))
  show("simulate re_formula NA", simulate(f, nsim = 2, re_formula = NA))
  show("posterior_linpred", frm_linpred(f))
  show("mfx avg_slopes", marginaleffects::avg_slopes(f))
  show("mfx predictions", marginaleffects::predictions(f))
  show("emmeans dpar mub", emmeans::emmeans(f, ~ x, dpar = "mub"))
  show("tidy", generics::tidy(f))
  show("glance", generics::glance(f))
  show("augment", generics::augment(f))
  show("summary print", capture.output(print(summary(f))))
  show("drop1", drop1(f))
  show("update REML", update(f, REML = TRUE))
}
