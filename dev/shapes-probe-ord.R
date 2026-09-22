# The ordinal fixture's parameter components: what variables() misses
# (defect S6 of dev/brmsport-findings.md) and where those values live.
#   Rscript dev/shapes-probe-ord.R > dev/shapes-log/probe-ord.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(testthat); library(frmtmb) })
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())

f4 <- brms_fixture(4)
cat("== names(par_template)\n"); print(names(f4$frame[["par_template"]]))
cat("== lengths\n"); print(lengths(f4$frame[["par_template"]]))
cat("== names(estimates)\n"); print(names(f4$estimates))
cat("== estimates\n"); print(f4$estimates[!names(f4$estimates) %in% "b"])
cat("== extra_names\n"); print(f4$frame[["extra_names"]])
cat("== confint rows\n"); print(rownames(confint(f4)))
cat("== outer_par_names\n"); print(frmtmb:::outer_par_names(f4))
cat("== rownames(cov.fixed)\n")
print(rownames(frmtmb:::sdr_of(f4)$cov.fixed))
cat("== linpreds\n")
for (lp in f4$frame[["linpreds"]]) {
  cat(" key:", lp[["resp"]], lp[["dpar"]], " par:", lp[["par"]],
      " cols:", paste(colnames(lp[["X"]]), collapse = ","),
      " cs:", length(lp[["cs"]] %||% list()), "\n")
  for (ct in lp[["cs"]] %||% list()) {
    cat("   cs term: label=", ct$label, " par=", ct$par, "\n", sep = "")
  }
}
cat("== fitted(f4) dim\n"); print(dim(fitted(f4)))
cat("== ordinal_ncat\n"); print(frmtmb:::ordinal_ncat(f4))
cat("== a cumulative fit with no cs, for the simple case\n")
set.seed(1)
d <- data.frame(x = rnorm(120))
d$y <- factor(cut(1 + 0.8 * d$x + rnorm(120), c(-Inf, -0.5, 0.5, Inf),
                  labels = 1:3), ordered = TRUE)
fo <- frm(bf(y ~ x) + cumulative(), data = d)
print(names(fo$estimates))
print(fo$estimates[["tau_raw"]])
print(variables(fo))
print(rownames(confint(fo)))
print(frmtmb:::outer_par_names(fo))
