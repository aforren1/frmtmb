# Reviewer: p167 poisson x arma(), brms spelling with cov = TRUE.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/review-equiv-cache",
                                             winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
s <- E$eq_specs[["poisson x arma()"]]()
d <- s$data
tr <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) e)
b <- tr(brms::make_stancode(brms::bf(y ~ x + arma(time = t, gr = g,
                                                         cov = TRUE)),
                            family = stats::poisson(), data = d))
cat("brms poisson arma(cov = TRUE):", if (inherits(b, "error"))
  conditionMessage(b) else "builds", "\n")
if (!inherits(b, "error")) {
  cc <- tr({E$brms_stan_model(b); TRUE})
  cat("compile:", if (inherits(cc, "error")) conditionMessage(cc) else "ok",
      "\n")
}
b0 <- tr(brms::make_stancode(brms::bf(y ~ x + arma(time = t, gr = g)),
                             family = stats::poisson(), data = d))
cat("brms poisson arma() (lane spelling, cov = FALSE):",
    if (inherits(b0, "error")) conditionMessage(b0) else "builds", "\n")
