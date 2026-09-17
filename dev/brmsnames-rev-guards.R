## Reviewer, claim 8: revert one fix at a time in memory and run the
## pinning test file.
##   Rscript dev/brmsnames-rev-guards.R none|psfit|refuse|residse
##   psfit   posterior_summary() on a fit falls through to the default
##   refuse  fit_refuse_draws_args() refuses nothing
##   residse residual__ Est.Error loses its betad derivative (reads a
##           constant), so it prints 0
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
q <- function(e) suppressWarnings(suppressMessages(e))
mut <- commandArgs(trailingOnly = TRUE)[1L]
q(library(testthat)); q(library(frmtmb))
ns <- asNamespace("frmtmb")
if (mut == "psfit") {
  f <- function(x, ...) get("posterior_summary.default", ns)(x, ...)
  registerS3method("posterior_summary", "frmtmb_fit", f, envir = ns)
  if (requireNamespace("brms", quietly = TRUE))
    registerS3method("posterior_summary", "frmtmb_fit", f, envir = asNamespace("brms"))
}
if (mut == "refuse") {
  assignInNamespace("fit_refuse_draws_args", function(...) invisible(NULL), "frmtmb")
}
if (mut == "residse") {
  orig <- get("varcorr_values", ns)
  g <- function(fit, vals, comp, layout = varcorr_layout(fit)) {
    out <- orig(fit, vals, comp, layout)
    if (!is.null(out$residual__)) {
      out$residual__$sd <- layout$residual$linkinv(
        fit$estimates$betad[layout$residual$betad])
    }
    out
  }
  environment(g) <- ns
  environment(g) <- list2env(list(orig = orig), parent = ns)
  assignInNamespace("varcorr_values", g, "frmtmb")
}
set.seed(1)
dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
dd$y <- 1 + dd$x + rnorm(20)[dd$g] + rnorm(200)
fit <- q(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
vc <- VarCorr(fit)
V <- vcov(fit, full = TRUE)
bn <- grep("sigma", rownames(V), value = TRUE)
cat("mutation", mut, ": residual Est.Error", vc$residual__$sd[, "Est.Error"],
    " delta method", exp(fit$estimates$betad[1]) * sqrt(V[bn, bn]), "\n")
res <- q(test_file("tests/testthat/test-brms-names.R", reporter = "silent",
                   package = "frmtmb"))
df <- as.data.frame(res)
cat(sprintf("test-brms-names.R: blocks %d PASS %d FAIL %d ERROR %d\n", nrow(df),
            sum(df$passed), sum(df$failed), sum(df$error)))
for (i in which(df$failed > 0 | df$error)) cat("  failing:", df$test[i], "\n")
