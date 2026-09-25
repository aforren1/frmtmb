# The two registry and refusal pins of lane fams, run against a chosen
# build, so each can be seen to fail on the base build and pass on the
# lane build.
#
#   FRMTMB_LIB=base Rscript dev/fams-pins.R
#   FRMTMB_LIB=/opt/rlib/lane-fams Rscript dev/fams-pins.R
lib <- Sys.getenv("FRMTMB_LIB", "base")
.libPaths(c(if (lib != "base") lib, "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")

# 1. frm_compat() on simulate() for the families that have a simulator
for (f in c("hurdle_poisson", "compois", "tweedie", "cox")) {
  cat(sprintf("  simulate x %-15s registry %-8s sim_can %s\n", f,
              frm_compat(f, "simulate")$status,
              frmtmb:::sim_can(frmtmb:::family_registry[[f]]())))
}

# 2. residuals(type = "osa") on a hurdle fit
set.seed(1)
d <- data.frame(x = rnorm(80))
d$y <- ifelse(runif(80) < 0.3, 0, rpois(80, 2) + 1)
fit <- frm(bf(y ~ x) + hurdle_poisson(), data = d)
cat("  residuals_osa x hurdle_poisson registry",
    frm_compat("hurdle_poisson", "residuals_osa")$status, "\n")
cat("  residuals(type = \"osa\"):",
    tryCatch({
      residuals(fit, type = "osa")
      "ran"
    }, error = function(e) conditionMessage(e)), "\n")
