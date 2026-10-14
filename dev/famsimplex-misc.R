# Remaining routes on the lane build: fitted(scale = "linear"), emmeans,
# frm_linpred(type = "link"), mvbf() with a dirichlet response, update().
# Output: dev/famsimplex-misc.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
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
show("fitted linear", fitted(fit, scale = "linear"))
show("fitted dpar mub", fitted(fit, dpar = "mub"))
show("fitted newdata", fitted(fit, newdata = dd[1:5, ]))
show("frm_linpred link", frm_linpred(fit, type = "link"))
show("update", fixef(update(fit, . ~ . + z)))
dd$g <- factor(rep(1:10, 30))
show("emmeans", if (requireNamespace("emmeans", quietly = TRUE)) {
  print(emmeans::emmeans(fit, ~ 1, dpar = "mub")); "ok" })
d2 <- dd
d2$y <- rnorm(nrow(d2))
show("mvbf dirichlet + gaussian", {
  f2 <- frm(mvbf(bf(Y ~ x, family = dirichlet()),
                 bf(y ~ x, family = gaussian())), data = d2)
  print(f2); fitted(f2, resp = "Y")
})
dl <- fs_data_logistic_normal(seed = 13)
dl$y <- rnorm(nrow(dl))
show("mvbf logistic_normal", frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                                      bf(y ~ x, family = gaussian())),
                                 data = dl))
show("confint ln", print(confint(frm(bf(Y ~ x), family = logistic_normal(),
                                     data = dl))))
show("diagnose", print(diagnose(fit)))
show("vcov", vcov(fit))
show("nobs", nobs(fit))
show("ranef-free coef", coef(fit))
