# Reviewer: a null direction spread over two crossed factors, so that no
# parameter's loading reaches se_load_tol = 0.1. y ~ a + b with
# a ~ 0 + f and b ~ 0 + g (k levels each, crossed): only a_i + b_j is
# identified; a shift a + t, b - t leaves the likelihood unchanged.
#   Rscript dev/nanse-rev-spread2.R base|lane|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r5",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
for (k in c(20L, 60L)) {
  set.seed(3)
  d <- expand.grid(f = factor(seq_len(k)), g = factor(seq_len(k)))
  d <- d[sample(nrow(d), 6 * k), ]
  d$y <- rnorm(k)[d$f] + rnorm(k)[d$g] + rnorm(nrow(d), 0, 0.5)
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(y ~ a + b, a ~ 0 + f, b ~ 0 + g, nl = TRUE), data = d),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  se <- suppressWarnings(fixef(fit)[, "Est.Error"])
  lost <- ns$sdr_of(fit)$se_lost
  cat(sprintf(paste0("k = %d: code %d, %d coefficients, finite SEs %d, ",
                     "median SE %.3g, lost %d, warnings %d\n"),
              k, fit$opt$convergence, length(se), sum(is.finite(se)),
              median(se[is.finite(se)]), length(lost), length(w)))
  for (x in w) cat("   warn:", substr(x, 1, 150), "\n")
  # what a user reads: a "significant" a_1 whose value is arbitrary
  s <- fixef(fit)
  cat("   ", rownames(s)[1], "estimate", signif(s[1, "Estimate"], 3),
      "Est.Error", signif(s[1, "Est.Error"], 3), "\n")
}
