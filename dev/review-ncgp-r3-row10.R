# Short check, item 4: does row 10a of test-brms-likelihood.R, at its loosened bound, still catch a real constant error? One mutant per
# process, in the namespace, then the row alone. Args: mutant
m <- commandArgs(TRUE)[1]
LIB <- "C:/Users/adf44/source/r/wt-ncgp-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true", R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win", FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/stan-cache")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
ns <- asNamespace("frmtmb")
put <- function(nm, f) { environment(f) <- ns; unlockBinding(nm, ns); assign(nm, f, envir = ns); lockBinding(nm, ns) }
of <- ns$gp_field; oc <- ns$gp_corr
switch(m,
  control = NULL,
  nug11 = { unlockBinding("gp_nugget", ns); assign("gp_nugget", 1e-11, envir = ns) },
  nug13 = { unlockBinding("gp_nugget", ns); assign("gp_nugget", 1e-13, envir = ns) },
  sd1e4 = put("gp_field", function(z, theta, blk) (1 + 1e-4) * of(z, theta, blk)),
  sd1e2 = put("gp_field", function(z, theta, blk) (1 + 1e-2) * of(z, theta, blk)),
  rho1e4 = put("gp_corr", function(theta, blk) { theta[-1] <- theta[-1] + 1e-4; oc(theta, blk) }),
  flip = put("gp_field", function(z, theta, blk) { v <- of(z, theta, blk); v[1] <- 1.001 * v[1]; v }))
if (m == "zdens") {
  core <- get("covstruct_core", envir = ns)
  core[["gp"]]$nll <- function(b, theta, blk) sum(RTMB::dnorm(b, 0, 1.0001, log = TRUE))
  unlockBinding("covstruct_core", ns); assign("covstruct_core", core, envir = ns); lockBinding("covstruct_core", ns)
}
f <- "C:/Users/adf44/source/r/frmtmb-wt-ncgp/tests/testthat/test-brms-likelihood.R"
r <- test_file(f, package = "frmtmb", env = test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE, desc = "check C: row 10a, the exact gp() against brms's own program")
for (t in r) {
  res <- t$results
  cat(sprintf("R10 %-8s pass %d fail %d err %d\n", m, sum(vapply(res, inherits, NA, "expectation_success")), sum(vapply(res, inherits, NA, "expectation_failure")), sum(vapply(res, inherits, NA, "expectation_error"))))
  for (x in res) if (!inherits(x, "expectation_success")) cat("  ", substr(gsub("\n", " ", conditionMessage(x)), 1, 150), "\n")
}
cat("DONE\n")
