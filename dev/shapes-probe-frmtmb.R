# What frmtmb returns today for the shapes item 2.6f moves, on the same
# fixtures the 2.6b port uses. Paired with dev/shapes-probe-brms.R.
#   Rscript dev/shapes-probe-frmtmb.R > dev/shapes-log/probe-frmtmb.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat); library(frmtmb)
})
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())

show <- function(label, expr) {
  cat("\n== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(expr,
         warning = function(w) {
           cat("  [warning] ", conditionMessage(w), "\n", sep = "")
           invokeRestart("muffleWarning") },
         message = function(m) invokeRestart("muffleMessage")),
    error = function(e) {
      cat("  ERROR: ", conditionMessage(e), "\n", sep = ""); NULL })
  invisible(r)
}
dd <- function(x) cat("dim:", paste(dim(x), collapse = " x "),
                      " len:", length(x),
                      " class:", paste(class(x), collapse = ","), "\n")

fit1 <- brms_fixture(1); fit2 <- brms_fixture(2); fit4 <- brms_fixture(4)
fit5 <- brms_fixture(5); fit6 <- brms_fixture(6)

show("fitted(fit1)", { x <- fitted(fit1); dd(x); print(utils::head(x, 3)) })
show("fitted(fit4)", { x <- fitted(fit4); dd(x); print(utils::head(x, 2)) })
show("residuals(fit2)", { x <- residuals(fit2); dd(x) })
show("predict(fit1)", { x <- predict(fit1); dd(x) })
show("predict(fit4)", { x <- predict(fit4); dd(x); print(colnames(x)) })
show("fixef(fit1)", { x <- fixef(fit1); str(x) })
show("fixef(fit1, flatten = TRUE)", print(names(fixef(fit1, flatten = TRUE))))
show("variables(fit1)", print(variables(fit1)))
show("variables(fit4)", print(variables(fit4)))
show("ngrps(fit1)", str(ngrps(fit1)))
show("ngrps(fit6)", str(ngrps(fit6)))
show("vcov(fit1)", { x <- vcov(fit1); dd(x); print(rownames(x)) })
show("vcov(fit1, full = TRUE) rows", print(rownames(vcov(fit1, full = TRUE))))
show("confint(fit1) rows", print(rownames(confint(fit1))))
show("summary(fit1) names", { s <- summary(fit1); print(names(s)) })
show("nobs", print(c(nobs(fit1), nobs(fit2), nobs(fit4), nobs(fit5))))
show("brms_coef_table(fit1)",
     print(frmtmb:::brms_coef_table(fit1)[, c("internal", "brms", "natural",
                                              "dpar")]))
show("brms_coef_table(fit4)",
     print(frmtmb:::brms_coef_table(fit4)[, c("internal", "brms", "natural",
                                              "dpar")]))
show("estimated_coef_names(fit1)",
     print(frmtmb:::estimated_coef_names(fit1)))
show("re_blocks fit1", {
  for (bk in fit1$frame[["re_blocks"]]) {
    cat(" block:", bk[["term_label"]], " group:", bk[["group_name"]],
        " dim:", bk[["dim"]], " levels:", bk[["n_levels"]],
        " covstruct:", bk[["covstruct"]], "\n")
  } })
show("varcorr_matrices(fit1)", print(frmtmb:::varcorr_matrices(fit1)))
show("varcorr_matrices(fit6)", print(frmtmb:::varcorr_matrices(fit6)))
show("print(summary(fit6))", print(summary(fit6)))
show("print(fit1)", print(fit1))
