# Lane remeasure, bin 2T: the constructions behind the verdicts of the
# coef, loo_compare and logLik rows (brmsfit-methods:136, 140, 142,
# 434, 469, 471). Seeds are the fixtures' own (helper-brms-suite-draws.R).
# brms's side reads its stored fixtures through brms:::rename_pars(), as
# brms's own test file does.
#
#   Rscript dev/remeasure-bin2-probe1.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample); library(testthat)})
`%||%` <- function(a, b) if (is.null(a)) b else a
source("tests/testthat/helper-brms-suite.R")
source("extensions/frmtmb.sample/tests/testthat/helper-brms-suite-draws.R")
cat("frmtmb", format(packageVersion("frmtmb")), "frmtmb.sample",
    format(packageVersion("frmtmb.sample")), "loo",
    format(packageVersion("loo")), "\n")
bx <- function(k) {
  brms:::rename_pars(get(paste0("brmsfit_example", k), asNamespace("brms")))
}
f1 <- brms_fixture(1); f2 <- brms_fixture(2); f4 <- brms_fixture(4)

cat("\n## coef(fit1) on the ML fit\n")
c1 <- coef(f1)
print(class(c1$visit)); print(dim(c1$visit)); print(colnames(c1$visit))
cat("brms coef(fit1)$visit:", dim(suppressMessages(coef(bx(1)))$visit),
    "\n")
print(dimnames(suppressMessages(coef(bx(1)))$visit)[[3]])

cat("\n## coef(fit2) on the ML fit (nonlinear a, b)\n")
str(coef(f2), max.level = 2)
b2 <- suppressMessages(coef(bx(2)))
cat("brms coef(fit2)$patient:", dim(b2$patient), "\n")
print(dimnames(b2$patient)[[3]])

cat("\n## coef(fit4): brms, the frmtmb draws, the frmtmb ML fit\n")
b4 <- suppressMessages(coef(bx(4)))
cat("brms:", dim(b4$subject), "\n"); print(dimnames(b4$subject)[[3]])
cd4 <- coef(brms_fixture_draws(4))
cat("frmtmb draws:", dim(cd4$subject), "\n"); print(dimnames(cd4$subject)[[3]])
cat("frmtmb fit:", dim(coef(f4)$subject), "\n")
print(colnames(coef(f4)$subject))

cat("\n## loo(fit6) against loo(fit6, newdata = fit6$data)\n")
d6 <- brms_fixture_draws(6)
cat("is d6$data NULL:", is.null(d6$data), "\n")
l1 <- suppressWarnings(loo(d6, cores = 1))
l2 <- suppressWarnings(loo(d6, cores = 1, newdata = d6$data))
cat("frmtmb elpd_loo, without and with newdata = d6$data:",
    l1$estimates[1, 1], l2$estimates[1, 1], "\n")
lc <- loo::loo_compare(l1, l2)
print(class(lc)); print(lc)
cat("frmtmb loo_compare(...)[2, 1]:", format(lc[2, 1]), "\n")
bf6 <- bx(6)
bl1 <- suppressWarnings(loo(bf6, cores = 1))
bl2 <- suppressWarnings(loo(bf6, cores = 1, newdata = bf6$data))
blc <- loo::loo_compare(bl1, bl2)
print(class(blc)); print(blc)
cat("brms loo_compare(...)[2, 1]:", format(blc[2, 1]), "\n")
cat("\n## loo(fit6 draws, newdata = the fixture's data)\n")
print(tryCatch(loo(d6, cores = 1, newdata = brms_fixture(6)$data)$estimates,
               error = function(e) paste("ERROR:", conditionMessage(e))))

cat("\n## logLik on draws, and brms's\n")
print(tryCatch(dim(logLik(brms_fixture_draws(1))),
               error = function(e) paste("ERROR:", conditionMessage(e))))
cat("brms dim(logLik(fit1)):", dim(logLik(bx(1))), "\n")

cat("\n## does loo() on draws read newdata? half the rows\n")
nd <- brms_fixture(6)$data[1:20, ]
r <- tryCatch(suppressWarnings(loo(d6, cores = 1, newdata = nd)),
              error = function(e) e)
if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n") else
  cat("pointwise rows:", nrow(r$pointwise), " elpd_loo:",
      r$estimates[1, 1], "\n")
print(names(formals(getS3method("loo", "frmtmb_draws"))))
