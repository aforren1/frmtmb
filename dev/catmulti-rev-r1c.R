# Reviewer, fix round 1: conditional_effects(categorical = TRUE) on a
# categorical(refcat = "b") fit with (1 | p | g) against brms 2.23.0
# (dev/catmulti-rev-out/brms-cat3-refb.rds), and the Y_a lp check
# perturbed. Log dev/catmulti-rev-out/r1c.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib", "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win", NOT_CRAN = "true",
           FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-catmulti/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
source("dev/catmulti-rev-data.R")
d3 <- rev_cat_data(21, 3)
fm <- frm(bf(y ~ x + (1 | p | g)), family = categorical(refcat = "b"), data = d3)
cf <- suppressMessages(conditional_effects(fm, "x", categorical = TRUE))[[1]]
bm <- readRDS("dev/catmulti-rev-out/brms-cat3-refb.rds")
cb <- suppressMessages(brms::conditional_effects(bm, "x", categorical = TRUE))[[1]]
cat("rows frmtmb", nrow(cf), "brms", nrow(cb), "cats", levels(cf$cats__), "|", levels(cb$cats__), "\n")
k <- match(paste(signif(cb$x, 10), cb$cats__), paste(signif(cf$x, 10), cf$cats__))
cat("matched rows:", sum(!is.na(k)), "\n")
d <- abs(cf$estimate__[k] - cb$estimate__) / cb$se__
r <- cf$se__[k] / cb$se__
cat(sprintf("|frmtmb - brms| / brms se: median %.3f max %.3f; se ratio median %.3f range %.3f to %.3f\n",
            median(d, na.rm = TRUE), max(d, na.rm = TRUE), median(r, na.rm = TRUE),
            min(r, na.rm = TRUE), max(r, na.rm = TRUE)))
nd <- cf[cf$cats__ == "c", ]
ff <- fitted(fm, newdata = nd, re_formula = NA)
cat("ce vs fitted(newdata, re_formula = NA): est rel", max(abs(nd$estimate__ - ff[, "Estimate", "P(Y = c)"]) /
    ff[, "Estimate", "P(Y = c)"]), "se rel", max(abs(nd$se__ - ff[, "Est.Error", "P(Y = c)"]) /
    ff[, "Est.Error", "P(Y = c)"]), "\n")

env <- testthat::test_env("frmtmb")
for (f in list.files("tests/testthat", "^helper-.*[.]R$", full.names = TRUE)) sys.source(f, env)
set.seed(47); n <- 300
dd <- data.frame(x = stats::rnorm(n))
W <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(1, .4, .4, 1), 2))
E <- cbind(0, 0.4 + 0.5 * dd$x + 0.6 * W[, 1], -0.3 - 0.4 * dd$x + 0.9 * W[, 2])
Y <- exp(E) / rowSums(exp(E)); colnames(Y) <- c("a", "b", "c")
dd$Y_a <- Y; dd$w <- 1 + 0.5 * dd$x + stats::rnorm(n)
fit <- frm(mvbf(bf(Y_a ~ x, family = logistic_normal()), bf(w ~ x, family = gaussian())), data = dd,
           control = frmtmb_control(grad_tol = 1e-6, restarts = 3))
fix <- function(code) {
  out <- gsub("multiply_lower_tri_self_transpose(Llncor);", "multiply_lower_tri_self_transpose(Llncor_Ya);",
              code, fixed = TRUE)
  gsub("[nlncor] lncor;", "[nlncor_Ya] lncor;", out, fixed = TRUE)
}
bform <- brms::bf(Y_a ~ x, family = brms::logistic_normal()) + brms::bf(w ~ x, family = gaussian()) +
  brms::set_rescor(FALSE)
run <- function(f, nm) {
  r <- tryCatch({ x <- env$brms_lp_check(bform, list(brms::logistic_normal(), gaussian()), dd, f,
                                         code_fix = fix)
                  sprintf("PASSED const %.3e grad %.3e", x$measured_const, x$max_grad) },
                error = function(e) paste("FAILED:", gsub("\n", " ", conditionMessage(e))))
  cat(nm, ":", substr(r, 1, 200), "\n")
}
run(fit, "unperturbed")
cat("estimate blocks:", names(fit$estimates), "\n")
for (blk in c("Y_a_lncor_raw", "betad", "beta")) {
  if (is.null(fit$estimates[[blk]])) next
  for (j in seq_along(fit$estimates[[blk]])) {
    fp <- fit; fp$estimates[[blk]][j] <- fp$estimates[[blk]][j] + 0.05
    run(fp, paste0("perturb ", blk, "[", j, "] (", names(fit$estimates[[blk]])[j] %||% "", ")"))
  }
}
cat("DONE\n")
