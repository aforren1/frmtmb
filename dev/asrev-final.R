## REVIEW, last probes: the predict() scale on a brms-shaped call, and
## predictive_interval()'s formals here against what brms ACCEPTS (not
## just what it declares).
.libPaths(c("C:/Users/adf44/source/r/asrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))

cat("==== predictive_interval / pp_check formals, ours ====\n")
for (nm in c("predictive_interval.frmtmb_draws", "pp_check.frmtmb_draws",
             "pp_check.frmtmb_fit")) {
  cat(sprintf("%-34s %s\n", nm,
              paste(names(formals(getFromNamespace(nm, if (grepl("draws", nm))
                "frmtmb.sample" else "frmtmb"))), collapse = ", ")))
}
cat("\nbrms predictive_interval.brmsfit forwards its dots to:\n")
print(body(getFromNamespace("predictive_interval.brmsfit", "brms")))
cat("posterior_predict.brmsfit has re.form: ",
    "re.form" %in% names(formals(getFromNamespace("posterior_predict.brmsfit",
                                                  "brms"))), "\n")
cat("=> brms ACCEPTS predictive_interval(x, re.form = ) through dots.\n")

cat("\n==== the ported brms call predict(fit, re_formula = NA) ====\n")
set.seed(11)
n <- 90
dd <- data.frame(x = rnorm(n), g = factor(rep(1:9, each = 10)))
eta <- 0.3 + 0.5 * dd$x + rnorm(9, 0, .6)[dd$g]
dd$y <- rpois(n, exp(eta))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
a <- predict(fit, re_formula = NA)
b <- predict(fit, re_formula = NA, type = "response")
cat(sprintf("predict(fit, re_formula = NA)[1]                  = %.6f\n", a[1]))
cat(sprintf("predict(fit, re_formula = NA, type='response')[1] = %.6f\n", b[1]))
cat(sprintf("ratio b/a                                        = %.4f\n",
            b[1] / a[1]))
cat("brms's posterior_epred(re_formula = NA) is the RESPONSE scale.\n")
cat("At 0.57.0 this call warned 'ignoring unknown arguments to",
    "predict(): re_formula'; now it is silent.\n")
