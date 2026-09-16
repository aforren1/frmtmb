# Spot check, part 2: the RETURN SHAPES brms's own suite asserts, run
# against one frmtmb fit. Every assertion below is copied from
# tests.brmsfit-methods.R with nothing changed but the fit it runs on
# and the numbers that depend on this design (nobs, ngrps).
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
load("dev/brms-suite/brms/data/epilepsy.rda")

fit <- frm(count ~ zBase * Trt + (1 | patient), epilepsy, family = poisson())
nobs_ <- nrow(epilepsy)
npatients <- length(unique(epilepsy$patient))
cat("fit:", nobs_, "rows,", npatients, "patients\n\n")

res <- list()
chk <- function(label, expr) {
  out <- tryCatch({
    withCallingHandlers(expr,
                        warning = function(w) invokeRestart("muffleWarning"),
                        message = function(m) invokeRestart("muffleMessage"))
    "PASS"
  }, error = function(e) {
    msg <<- conditionMessage(e)
    if (inherits(e, "expectation_failure")) "FAIL" else "ERROR"
  })
  res[[length(res) + 1L]] <<- data.frame(
    what = label, outcome = out,
    message = if (identical(out, "PASS")) "" else substr(gsub("[\r\n]+", " ", msg), 1, 160),
    stringsAsFactors = FALSE)
}
msg <- ""

# --- brms: "fitted has reasonable outputs" ----------------------------
chk("dim(fitted(fit)) == c(nobs, 4)",
    expect_equal(dim(fitted(fit)), c(nobs_, 4)))
chk("colnames(fitted(fit)) == Estimate/Est.Error/Q2.5/Q97.5",
    expect_equal(colnames(fitted(fit)),
                 c("Estimate", "Est.Error", "Q2.5", "Q97.5")))
chk("fitted(fit, newdata = d[1:2, ]) is 2 x 4",
    expect_equal(dim(fitted(fit, newdata = epilepsy[1:2, ])), c(2, 4)))
chk("fitted(fit, re_formula = NA) runs",
    expect_true(is.numeric(fitted(fit, re_formula = NA))))
chk("fitted(fit, dpar = 'inv') errors with \"Invalid argument 'dpar'\"",
    expect_error(fitted(fit, dpar = "inv"), "Invalid argument 'dpar'"))
chk("fitted(fit, scale = 'linear') differs from response scale",
    expect_true(!isTRUE(all.equal(fitted(fit),
                                  fitted(fit, scale = "linear")))))

# --- brms: "predict has reasonable outputs" ---------------------------
chk("dim(predict(fit)) == c(nobs, 4)",
    expect_equal(dim(predict(fit)), c(nobs_, 4)))
chk("colnames(predict(fit)) == Estimate/Est.Error/Q2.5/Q97.5",
    expect_equal(colnames(predict(fit)),
                 c("Estimate", "Est.Error", "Q2.5", "Q97.5")))
chk("predict(fit, probs = c(.2,.5,.8)) has 5 columns",
    expect_equal(ncol(predict(fit, probs = c(0.2, 0.5, 0.8))), 5))
chk("predict(fit, newdata, allow_new_levels = TRUE) is 2 x 4", {
  nd <- epilepsy[1:2, ]
  nd$patient <- factor(c("new1", "new2"))
  expect_equal(dim(predict(fit, newdata = nd, allow_new_levels = TRUE)),
               c(2, 4))
})

# --- brms: "residuals has reasonable outputs" -------------------------
chk("dim(residuals(fit)) == c(nobs, 4)",
    expect_equal(dim(residuals(fit)), c(nobs_, 4)))
chk("residuals(fit, type = 'pearson', probs = 0.65) is nobs x 3",
    expect_equal(dim(residuals(fit, type = "pearson", probs = c(0.65))),
                 c(nobs_, 3)))

# --- brms: ngrps, nobs, model.frame, formula, family ------------------
chk("ngrps(fit) == list(patient = npatients)",
    expect_equal(ngrps(fit), list(patient = npatients)))
chk("nobs(fit) == nobs", expect_equal(nobs(fit), nobs_))
chk("model.frame(fit) equals fit$data",
    expect_equal(model.frame(fit), fit$data))
chk("family(fit) inherits from brmsfamily",
    expect_true(is(family(fit), "brmsfamily")))
chk("print(family(fit), links = TRUE) shows the link",
    expect_output(print(family(fit), links = TRUE), "poisson.*log"))

# --- brms: "fixef/ranef/coef/VarCorr have reasonable outputs" ---------
chk("rownames(fixef(fit))",
    expect_equal(rownames(fixef(fit)),
                 c("Intercept", "zBase", "Trt1", "zBase:Trt1")))
chk("dim(ranef(fit)$patient) == c(npatients, 4, 1)",
    expect_equal(dim(ranef(fit)$patient), c(npatients, 4, 1)))
chk("dim(coef(fit)$patient) == c(npatients, 4, 4)",
    expect_equal(dim(coef(fit)$patient), c(npatients, 4, 4)))
chk("names(VarCorr(fit)) == 'patient'",
    expect_equal(names(VarCorr(fit)), "patient"))
chk("dimnames(VarCorr(fit)$patient$cov)[c(1,3)]",
    expect_equal(dimnames(VarCorr(fit)$patient$cov)[c(1, 3)],
                 list("Intercept", "Intercept")))

# --- brms: "summary / variables / vcov / log_lik" ---------------------
chk("colnames(summary(fit)$fixed)",
    expect_equal(colnames(summary(fit)$fixed),
                 c("Estimate", "Est.Error", "l-95% CI", "u-95% CI",
                   "Rhat", "Bulk_ESS", "Tail_ESS")))
chk("print(summary(fit)) says 'Regression Coefficients:'",
    expect_output(print(summary(fit)), "Regression Coefficients:"))
chk("print(summary(fit)) says 'Multilevel Hyperparameters:'",
    expect_output(print(summary(fit)), "Multilevel Hyperparameters:"))
chk("variables(fit) has b_Intercept and sd_patient__Intercept",
    expect_true(all(c("b_Intercept", "sd_patient__Intercept") %in%
                      variables(fit))))
chk("dim(vcov(fit)) == c(4, 4)", expect_equal(dim(vcov(fit)), c(4, 4)))
chk("dim(vcov(fit, cor = TRUE)) == c(4, 4)",
    expect_equal(dim(vcov(fit, cor = TRUE)), c(4, 4)))
chk("posterior_summary(fit) has 4 columns",
    expect_equal(ncol(posterior_summary(fit)), 4))

# --- brms: "hypothesis has reasonable outputs" ------------------------
chk("hypothesis(fit, 'zBase > Trt1')$hypothesis is 1 x 8",
    expect_equal(dim(hypothesis(fit, "zBase > Trt1")$hypothesis), c(1, 8)))
chk("hypothesis(fit, 'Intercept = 0', class = 'sd', group = 'patient')",
    expect_equal(dim(hypothesis(fit, "Intercept = 0", class = "sd",
                                group = "patient")$hypothesis), c(1, 8)))
chk("hypothesis(fit, 1) errors on a non-character hypothesis",
    expect_error(hypothesis(fit, 1),
                 "Argument 'hypothesis' must be a character vector"))
chk("hypothesis(fit, 'Intercept > x') names the missing term",
    expect_error(hypothesis(fit, "Intercept > x"), "b_x", fixed = TRUE))

# --- brms: log_lik, loo, waic, pp_check, conditional_effects ----------
chk("log_lik(fit) has nobs columns",
    expect_equal(ncol(log_lik(fit)), nobs_))
chk("is.numeric(loo(fit)$estimates)",
    expect_true(is.numeric(loo(fit)$estimates)))
chk("is.numeric(waic(fit)$estimates)",
    expect_true(is.numeric(waic(fit)$estimates)))
chk("pp_check(fit) returns a ggplot",
    expect_true(is(pp_check(fit), "ggplot")))
chk("nrow(conditional_effects(fit)[[1]]) == 100",
    expect_equal(nrow(conditional_effects(fit)[[1]]), 100))
chk("conditional_effects(fit, 'nope') errors",
    expect_error(conditional_effects(fit, effects = "nope"),
                 "All specified effects are invalid for this model"))

tab <- do.call(rbind, res)
write.table(tab, "dev/brmssuite-spotcheck2.tsv", sep = "\t",
            row.names = FALSE, quote = TRUE)
print(table(tab$outcome))
cat("\n")
for (i in seq_len(nrow(tab))) {
  cat(sprintf("[%-5s] %s\n", tab$outcome[i], tab$what[i]))
  if (nzchar(tab$message[i])) cat("          ", tab$message[i], "\n")
}
