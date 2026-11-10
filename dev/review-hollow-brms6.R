# Reviewer: is brms's loo(fit6, newdata = fit6$data) difference all in the
# gp response? brms 2.23.0 from the user library, no compile needed.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
fit6 <- brms:::rename_pars(brms:::brmsfit_example6)
set.seed(469)
for (r in c("volume", "count")) {
  a <- log_lik(fit6, resp = r)
  b <- log_lik(fit6, resp = r, newdata = fit6$data)
  cat(sprintf("resp %-6s identical %s, max abs diff %.3g\n", r,
              identical(a, b), max(abs(a - b))))
}
# repeat the volume call: is brms's new-row gp prediction random?
b1 <- log_lik(fit6, resp = "volume", newdata = fit6$data)
b2 <- log_lik(fit6, resp = "volume", newdata = fit6$data)
cat(sprintf("volume newdata twice: max abs diff %.3g\n", max(abs(b1 - b2))))
