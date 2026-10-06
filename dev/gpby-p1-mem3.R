# Punch round 1, m5: where the 2000-row extra_cov peak sits, step by step.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
src <- parse(file.path(wt, "tests/testthat/test-gp-by.R"))
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e)
d <- gpby_data()
fit <- frm(bf(y ~ gp(x)), data = d)
big <- data.frame(x = seq(-1, 9, length.out = 2000))
lp <- fit$frame$linpreds[["y.mu"]]
own <- function(lab, expr) {
  invisible(gc())
  g0 <- gc(reset = TRUE)
  r <- eval.parent(substitute(expr))
  g1 <- gc()
  cat(sprintf("%-40s own peak %6.1f MB (one n x n is %.1f MB)\n", lab,
              g1[2, 6] - g0[2, 2], 2000^2 * 8 / 2^20))
  invisible(r)
}
ed <- own("lp_eta_design", ns$lp_eta_design(fit, lp, big, TRUE, FALSE))
krig <- Filter(Negate(is.null), lapply(ed$sm_parts, `[[`, "krig"))[[1]]
cat("unseen rows", length(krig$rows), "\n")
S <- own("gp_krig_cov", ns$gp_krig_cov(krig))
rm(S)
E <- own("lp_extra_cov", ns$lp_extra_cov(fit, ed, TRUE))
rm(E)
own("frm_lp_basis(extra_cov = TRUE)",
    frm_lp_basis(fit, newdata = big, extra_cov = TRUE))
own("frm_lp_basis()", frm_lp_basis(fit, newdata = big))
