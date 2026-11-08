# Every family of test-nl-family.R's nl_fams_single, its three fits, on
# the file's own data builders (sourced from the test file, so the probe
# reads what the test reads): the distributional parameters each fit
# stops at, the warnings it raises, and the file's ratio. A dpar that
# runs off to |value| > 8 on the log scale is a maximum at the boundary,
# where the stopping point, and so the warning, is the platform's
# rounding.
#
#   Rscript dev/ci-nlfamprobe.R <lib1;lib2;...> <test file> <out tsv>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[[1L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), "\n")
src <- parse(args[[2L]])
# the helpers and the family vector, not the test_that() blocks
keep <- vapply(src, function(e) {
  is.call(e) && identical(e[[1L]], as.name("<-"))
}, NA)
env <- new.env(parent = asNamespace("frmtmb"))
for (e in src[keep]) eval(e, env)
fit1 <- function(form, f, d) {
  w <- character()
  fit <- withCallingHandlers(
    suppressMessages(frm(form, family = env$nl_fam_of(f), data = d)),
    warning = function(cnd) {
      w <<- c(w, substr(conditionMessage(cnd), 1, 40))
      invokeRestart("muffleWarning")
    })
  p <- fit$opt$par
  dp <- p[names(p) == "betad"]
  list(ll = as.numeric(logLik(fit)), w = paste(unique(w), collapse = " | "),
       dmax = if (length(dp)) max(abs(dp)) else NA_real_,
       conv = fit$opt$convergence)
}
rows <- list()
for (f in env$nl_fams_single) {
  z <- env$nl_fam_data()
  d <- data.frame(x = z$x, ntr = 10L)
  d$y <- env$nl_fam_response(f, z)
  L <- env$nl_fam_lhs(f)
  r <- list(
    lin = fit1(bf(stats::as.formula(paste(L, "~ x"))), f, d),
    int = fit1(bf(stats::as.formula(paste(L, "~ 1"))), f, d),
    nl = fit1(bf(stats::as.formula(paste(L, "~ b0 + b1 * x")),
                 b0 + b1 ~ 1, nl = TRUE), f, d))
  ratio <- abs(r$nl$ll - r$lin$ll) / (r$lin$ll - r$int$ll)
  for (m in names(r)) {
    rows[[length(rows) + 1L]] <- data.frame(
      family = f, model = m, ll = sprintf("%.10f", r[[m]]$ll),
      dmax = sprintf("%.3g", r[[m]]$dmax), conv = r[[m]]$conv,
      ratio = sprintf("%.3g", ratio), warnings = r[[m]]$w)
  }
}
out <- do.call(rbind, rows)
utils::write.table(out, args[[3L]], sep = "\t", quote = FALSE,
                   row.names = FALSE)
cat("families:", length(unique(out$family)), "rows:", nrow(out), "\n")
