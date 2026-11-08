# Review of lane arstep: a constant sigma (bf(..., sigma = 0.3)) passes
# gp_marginal_gate() with has_sigma FALSE. Does frm() let warnings escape,
# and is the marginal start used, with and without ar()?
# Usage: Rscript review-arstep-sigconst.R lane|base
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
arm <- commandArgs(TRUE)[1]
rv_libs(arm)
ns <- asNamespace("frmtmb")
set.seed(12)
n <- 50
d <- data.frame(x = round(runif(n, 0, 10), 3), t = seq_len(n))
d$y <- sin(d$x) + 2 * (d$x > 5) +
  as.numeric(arima.sim(list(ar = 0.5), n, sd = 0.2))
for (f in list(bf(y ~ gp(x), sigma = 0.3),
               bf(y ~ gp(x) + ar(time = t), sigma = 0.3),
               bf(y ~ gp(x) + ar(time = t), sigma = 0.2))) {
  r <- rv_quiet(frm(f, data = d, family = gaussian()))
  fit <- r$value
  lbl <- paste(deparse(f$formula), collapse = "")
  if (inherits(fit, "error")) {
    cat("SC", arm, lbl, "error", conditionMessage(fit), "\n"); next
  }
  cat(sprintf("SC %s %-34s sigma %s ll %.6f conv %d marginal-used %s warnings %d\n",
              arm, lbl, format(attr(f, "sigma") %||% NA), as.numeric(logLik(fit)),
              fit$opt$convergence, !is.null(fit$obj$env$.frmtmb_gp_marginal),
              length(r$warnings)))
  for (w in unique(r$warnings)) cat("SCW", arm, lbl, substr(w, 1, 120), "\n")
}
cat("DONE\n")
