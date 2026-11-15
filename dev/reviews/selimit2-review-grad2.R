# Review of lane selimit2, item 4824 (continued): a start far along the joint
# end of poisson `y ~ g` (group a all zero) with group b's mean off, and a
# small iteration cap, so the stop has a gradient across the end.
#   Rscript selimit2-review-grad2.R <base|lane>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "\n")
set.seed(1)
d <- data.frame(g = factor(rep(c("a", "b"), each = 20)),
                y = c(rep(0, 20), stats::rpois(20, 5)))
mb <- mean(d$y[d$g == "b"])
one <- function(st, k, rt = NULL) {
  ws <- character(0)
  oc <- list(iter.max = k, eval.max = 3 * k)
  if (!is.null(rt)) oc <- list(rel.tol = rt)
  fit <- withCallingHandlers(
    frm(y ~ g, data = d, family = poisson(), start = list(beta = st),
        control = frmtmb_control(optCtrl = oc)),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  b <- fit$opt$par
  g <- drop(fit$obj$gr(b))
  # gradient across the end: along (1, 1)/sqrt(2) in raw units
  cat(sprintf(paste0("start (%g, %g) iter %s | code %d | par %.4f %.4f | ",
                     "grad %.3g %.3g | b mean %.5f (MLE %.5f, rel %.3g) | ",
                     "nwarn %d\n"),
              st[1], st[2], if (is.null(rt)) k else paste("rel.tol", rt),
              fit$opt$convergence, b[1], b[2], g[1], g[2], exp(sum(b)), mb,
              exp(sum(b)) / mb - 1, length(ws)))
  for (w in ws) cat("     W:", substr(w, 1, 200), "\n")
}
for (st in list(c(-25, 24.5), c(-25, 26), c(-25, 25.2), c(-30, 31.5))) {
  for (k in c(1, 2, 3)) one(st, k)
}
