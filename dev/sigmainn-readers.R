# Every reader of sigma on an ar(cov = TRUE) fit, before (rellib-r14)
# and after (the lane library): what each returns for sigma on the AR(1)
# 0.9 data of dev/acprior-agree.R (seed 20261027), and under sigma ~ x2
# on the sigx data of dev/sigmainn-route.R.
# Usage: Rscript dev/sigmainn-readers.R before|after
# Log: dev/sigmainn-out/readers-<arm>.txt
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-sigmainn-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
src <- readLines(file.path(wt, "dev", "acprior-agree.R"))
i0 <- grep("^make_data <- function", src)
i1 <- i0 + which(src[(i0 + 1):length(src)] == "}")[1]
eval(parse(text = src[i0:i1]))
d <- make_data()$ar
fit <- frm(y ~ x + ar(time = t, gr = g, cov = TRUE), data = d)
show <- function(lab, expr) {
  v <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.numeric(v)) v <- paste(format(signif(v, 6)), collapse = " ")
  cat(sprintf("%-46s %s\n", lab, paste(v, collapse = " ")))
}
nat <- frmtmb:::autocor_natural(fit$estimates$thetaac, fit$frame$autocor[[1]])
show("ar[1]", unname(nat))
show("exp(sigma_(Intercept))", exp(fit$estimates$betad[[1]]))
show("sigma(fit)", sigma(fit))
show("summary(fit)$spec_pars", unlist(summary(fit)$spec_pars["sigma", ]))
show("VarCorr(fit)$residual__$sd", VarCorr(fit)$residual__$sd)
show("confint(fit)[sigma_(Intercept), ]",
     exp(confint(fit)["sigma_(Intercept)", ]))
show("confint_varcorr(fit) sigma row", {
  cv <- confint_varcorr(fit)
  unlist(cv[grepl("sigma", rownames(cv)) | grepl("sigma", cv[[1]]), ][1, ])
})
show("vcov(full = TRUE) sigma row sd", sqrt(vcov(fit, full = TRUE)[3, 3]))
show("frm_linpred(dpar = 'sigma')[1:2]",
     frm_linpred(fit, dpar = "sigma", type = "response")[1:2])
show("fitted(dpar = 'sigma')[1, 1]", fitted(fit, dpar = "sigma")[1, 1])
show("predict(dpar = 'sigma')[1:2]",
     tryCatch(predict(fit, dpar = "sigma")[1:2], error = function(e) stop(e)))
show("posterior_epred(dpar = 'sigma')", tryCatch(
  mean(posterior_epred(fit, dpar = "sigma", ndraws = 20)),
  error = function(e) stop(e)))
show("residuals pearson sd", sd(residuals(fit, type = "pearson")[, 1]))
show("simulate() residual sd", {
  s <- as.matrix(simulate(fit, nsim = 400, seed = 1))
  sd(as.vector(s - fitted(fit)[, 1]))
})
show("posterior_predict() sd", {
  pp <- posterior_predict(fit, ndraws = 400)
  sd(as.vector(sweep(pp, 2, fitted(fit)[, 1])))
})
show("hypothesis(fit, 'sigma > 0')",
     tryCatch(hypothesis(fit, "sigma > 0")$hypothesis$Estimate,
              error = function(e) stop(e)))
show("variables(fit) with sigma", grep("sigma", variables(fit), value = TRUE))
show("insight::get_sigma", if (requireNamespace("insight", quietly = TRUE))
  as.numeric(insight::get_sigma(fit)) else NA)
show("logLik", as.numeric(logLik(fit)))
show("prior_summary rows", nrow(prior_summary(fit)))
show("frmtmb_autocor_pars(fit)", tryCatch(
  paste(capture.output(print(frmtmb_autocor_pars(fit))), collapse = " | "),
  error = function(e) stop(e)))
fx <- {
  set.seed(20261027)
  ng <- 20; k <- 6
  dx <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                   x = rnorm(ng * k), x2 = rnorm(ng * k))
  dx$y <- 1 + 0.5 * dx$x + exp(-0.8 + 0.3 * dx$x2) * rnorm(nrow(dx))
  frm(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x2), data = dx)
}
show("sigma ~ x2: fixef sigma rows", {
  fe <- fixef(fx)
  fe[grepl("sigma", rownames(fe)), "Estimate"]
})
show("sigma ~ x2: frm_linpred(dpar = sigma)[1:2]",
     frm_linpred(fx, dpar = "sigma", type = "response")[1:2])
show("sigma ~ x2: conditional_effects sigma", {
  ce <- conditional_effects(fx, effects = "x2", dpar = "sigma")[[1]]
  range(ce$estimate__)
})
if (identical(arm, "after")) {
  old <- fit
  old$frame$autocor[[1]]$sigma_is <- NULL
  for (m in c("summary", "sigma", "VarCorr", "fixef", "vcov", "confint",
              "fitted", "residuals", "simulate", "print", "logLik",
              "predict", "coef", "posterior_predict", "update")) {
    show(paste0("saved 0.76.0 fit: ", m),
         if (m == "update") as.numeric(logLik(update(old)))
         else {
           r <- do.call(m, list(old))
           if (is.numeric(r)) r[1] else class(r)[1]
         })
  }
}
