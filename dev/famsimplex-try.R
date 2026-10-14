# First fits of the simplex families on the lane build, every post-fit
# route. Output: dev/famsimplex-try.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
show <- function(lab, expr) {
  cat("\n---", lab, "\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r) && length(r) == 1) cat(r, "\n") else {
    cat("class:", class(r), " dim:", dim(r), "\n")
    if (!is.null(dimnames(r))) str(dimnames(r))
  }
  invisible(r)
}
dd <- fs_data_dirichlet(seed = 11)
fit <- frm(bf(Y ~ x, phi ~ z), family = dirichlet(), data = dd)
print(fit)
print(summary(fit))
show("fitted", f <- fitted(fit)); print(head(f[, , 1]))
show("predict", p <- predict(fit)); print(head(p[, , 1]))
show("simulate", s <- simulate(fit, nsim = 2))
show("frm_linpred dpar", frm_linpred(fit, dpar = "mub"))
show("frm_linpred response", frm_linpred(fit, type = "response"))
show("residuals", residuals(fit))
show("ce", conditional_effects(fit, "x"))
show("default_prior", print(default_prior(fit)))
dl <- fs_data_logistic_normal(seed = 13)
fit2 <- frm(bf(Y ~ x, sigmab ~ z), family = logistic_normal(), data = dl)
print(fit2)
show("ln fitted", fitted(fit2))
show("ln predict", predict(fit2))
d2 <- fs_data_dirichlet2(seed = 15)
fit3 <- frm(bf(Y ~ x), family = brmsfamily("dirichlet2"), data = d2)
print(fit3)
