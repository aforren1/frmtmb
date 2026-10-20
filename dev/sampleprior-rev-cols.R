# Reviewer of lane sampleprior: the draws matrix's columns on rellib-r9
# and on the lane (same model, seed), and downstream readers.
lane <- !identical(Sys.getenv("REV_ARM"), "r9")
.libPaths(c(if (lane) "C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb.sample"), "\n")
set.seed(20261118)
G <- 8
n <- 120
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), each = n / G)))
d$y <- 1 + rnorm(G)[d$g] + (1 + rnorm(G, 0, .5)[d$g]) * d$x + rnorm(n)
fit <- frm(y ~ x + (1 + x | g), data = d)
for (lap in c(FALSE, TRUE)) {
  ds <- suppressMessages(suppressWarnings(frm_sample(
    fit, chains = 1, iter = 300, refresh = 0, seed = 2, laplace = lap)))
  cn <- colnames(ds$draws)
  cat("laplace", lap, ": ", paste(cn[!grepl("^r_", cn)], collapse = " "), "\n")
  cat("  positions of theta_:", match(grep("^theta_", cn, value = TRUE), cn),
      " ncol", length(cn), "\n")
  r <- tryCatch({
    pp <- posterior_predict(ds, ndraws = 20)
    ll <- log_lik(ds)
    s <- summary(ds)
    vc <- VarCorr(ds)
    sprintf("posterior_predict %s, log_lik %s, summary ok, VarCorr sd %s",
            paste(dim(pp), collapse = "x"), paste(dim(ll), collapse = "x"),
            paste(sprintf("%.4f", vc$g$sd[, 1]), collapse = " "))
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  cat("  ", r, "\n")
}
