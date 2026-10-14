# Reviewer: routes the lane did not test: mvbf(), update(), newdata
# without the response, emmeans without dpar, weights.
# Output: dev/famsimplex-rev-routes.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
show <- function(lab, expr) {
  cat("\n---", lab, "\n")
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  if (length(w)) cat("WARNINGS:", unique(substr(w, 1, 300)), sep = "\n  ")
  if (is.character(r) && length(r) == 1) cat(r, "\n") else print(r)
  invisible(r)
}
dd <- fs_data_dirichlet(seed = 11)
dd$y2 <- 1 + dd$x + rnorm(nrow(dd))
f1 <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
show("mvbf dirichlet + gaussian", {
  f <- frm(mvbf(bf(Y ~ x, family = dirichlet()), bf(y2 ~ x)), data = dd)
  list(ll = c(logLik(f), logLik(f1) +
                logLik(frm(bf(y2 ~ x), family = gaussian(), data = dd))),
       coef = rownames(fixef(f)), fitted_dim = dim(fitted(f, resp = "Y")))
})
show("mvbf logistic_normal", frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                                      bf(y2 ~ x)), data = dd))
show("update data with renamed columns", {
  d2 <- dd
  colnames(d2$Y) <- c("p", "q", "r")
  rownames(fixef(update(f1, data = d2)))
})
show("update formula", rownames(fixef(update(f1, . ~ . + z))))
show("fitted newdata without Y", dim(fitted(f1, newdata = data.frame(x = 0:2))))
show("predict newdata without Y", dim(predict(f1, newdata = data.frame(x = 0:2))))
show("simulate newdata", {
  s <- simulate(f1, nsim = 1, seed = 1, newdata = data.frame(x = 0:2))
  s[[1]]
})
show("emmeans no dpar", emmeans::emmeans(f1, ~ x, at = list(x = 0)))
show("frm_linpred dpar phi", head(frm_linpred(f1, dpar = "phi")))
dd$w <- rep(c(1, 3), 150)
show("weights double the rows", {
  fw <- frm(bf(Y | weights(w) ~ x), family = dirichlet(), data = dd)
  d3 <- rbind(dd[dd$w == 1, ], dd[dd$w == 3, ], dd[dd$w == 3, ],
              dd[dd$w == 3, ])
  fr <- frm(bf(Y ~ x), family = dirichlet(), data = d3)
  c(weighted = as.numeric(logLik(fw)), replicated = as.numeric(logLik(fr)))
})
show("anova LRT phi ~ z", {
  f2 <- frm(bf(Y ~ x, phi ~ z), family = dirichlet(), data = dd)
  anova(f1, f2)
})
show("loo on ML fit", tryCatch(loo(f1), error = function(e) conditionMessage(e)))
show("log_lik on ML fit dims", dim(log_lik(f1)))
show("print fit with refcat c", frm(bf(Y ~ x), family = dirichlet(refcat = "c"),
                                    data = dd))
show("default_prior dirichlet", default_prior(bf(Y ~ x), data = dd,
                                              family = dirichlet()))
show("brms default_prior dirichlet", brms::default_prior(brms::bf(Y ~ x),
     data = dd, family = brms::dirichlet()))
dl <- fs_data_logistic_normal(seed = 13)
show("default_prior ln", default_prior(bf(Y ~ x), data = dl,
                                       family = logistic_normal()))
show("brms default_prior ln", brms::default_prior(brms::bf(Y ~ x), data = dl,
     family = brms::logistic_normal()))
show("brms set_prior class lncor accepted", {
  brms::validate_prior(brms::set_prior("lkj(2)", class = "lncor"),
                       brms::bf(Y ~ x), data = dl,
                       family = brms::logistic_normal())
})
show("brms set_prior class Llncor", {
  brms::validate_prior(brms::set_prior("lkj(2)", class = "Llncor"),
                       brms::bf(Y ~ x), data = dl,
                       family = brms::logistic_normal())
})
show("frm prior class Llncor", frm(bf(Y ~ x), family = logistic_normal(),
     data = dl, prior = set_prior("lkj(2)", class = "Llncor")))
