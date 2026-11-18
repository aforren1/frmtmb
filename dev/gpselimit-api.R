# Lane gpselimit: can the API report what a white field's data determine
# (sdgp^2 times the product of the length scales), and do its fitted values
# and predictions keep their standard errors? Sweep row 3 (white) and row 6
# (resolved).
# Usage: Rscript gpselimit-api.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-lib.R")
a <- commandArgs(TRUE)
gs_libs(a[1])
try1 <- function(lab, expr) {
  v <- tryCatch(gs_quiet(expr), error = function(e) e)
  if (inherits(v, "error")) {
    cat(lab, "ERROR:", conditionMessage(v), "\n")
  } else {
    cat(lab, "\n"); print(v)
    if (length(attr(v, "conds"))) cat("  conds:", attr(v, "conds"), "\n")
  }
}
for (row in c(3, 6)) {
  s <- gs_sweep_row(row)
  fit <- gs_fit(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian())$fit
  cat("==== row", row, "\n")
  print(names(fit$opt$par))
  print(frmtmb::variables(fit))
  try1("hyp sdgp", hypothesis(fit, "sdgp_gpx = 0", class = ""))
  try1("hyp product", hypothesis(fit,
                                 "2 * log(sdgp_gpx) + log(lscale_gpx) = 0",
                                 class = ""))
  try1("fitted", head(fitted(fit), 3))
  try1("predict", head(predict(fit, newdata = data.frame(x = c(2.5, 12)),
                               ndraws = 400, seed = 1), 2))
  try1("linpred se", head(frm_linpred(fit, se.fit = TRUE)$se.fit, 3))
}
