# Lane seflags, punch round 1: run one test file against the lane's
# build as it was before the punch round. Only R/se-check.R changed in
# the round's code, so the lane library is loaded and the pre-punch
# se-check.R (dev/seflags-out/prepunch-se-check.R, the file as install
# 9 had it, which the reviewer compared with the installed namespace)
# is evaluated into the namespace over it.
# SEFLAGS_OVERLAY names another se-check.R to evaluate instead: the
# punch round 1 file is dev/seflags-out/punch1-se-check.R, today's file
# with the round-2 changes to se_tier3(), jc_nonest() and
# jc_drop_boundary() taken out. Punch round 2b: dev/seflags-out/
# install22-overlay.R is R/se-check.R, R/predict.R, R/brms-shapes.R and
# R/confint.R as install 22 had them, one after the other. An active
# binding of the namespace is left as it is.
#   Rscript dev/seflags-prepunch.R <test file> [package]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-seflags-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
pkg <- if (length(a) >= 2) a[2] else "frmtmb"
suppressMessages({library(testthat); library(pkg, character.only = TRUE)})
ns <- asNamespace("frmtmb")
src <- Sys.getenv("SEFLAGS_OVERLAY",
  "C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-out/prepunch-se-check.R")
tmp <- new.env()
sys.source(src, envir = tmp)
n <- 0L
for (nm in ls(tmp, all.names = TRUE)) {
  if (exists(nm, envir = ns, inherits = FALSE) && !bindingIsActive(nm, ns)) {
    unlockBinding(nm, ns)
    f <- get(nm, envir = tmp)
    if (is.function(f)) environment(f) <- ns
    assign(nm, f, envir = ns)
    n <- n + 1L
  }
}
cat("pre-punch overlay:", n, "bindings from", src, "\n")
res <- test_file(a[1], package = pkg, env = test_env(pkg), reporter = "silent")
for (t in res) for (x in t$results) {
  k <- if (inherits(x, "expectation_failure")) "FAIL" else
    if (inherits(x, "expectation_error")) "ERROR" else NA
  if (!is.na(k)) cat("DETAIL", k, "[", t$test, "]",
                     substr(gsub("\n", " | ", conditionMessage(x)), 1, 300), "\n")
}
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            basename(a[1]), sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning)))
