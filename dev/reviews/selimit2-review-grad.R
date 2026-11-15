# Review of lane selimit2, item 4824: can se_boundary_stop() with a gradient
# above grad_tol accept a fit that is not at its limit? A joint end marks
# every coefficient of the combination "end", and check_convergence() then
# drops their whole gradient, not only its component along the end.
#   Rscript selimit2-review-grad.R <base|lane>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "\n")
run <- function(label, fo, d, fam, ctl, truth) {
  ws <- character(0)
  ms <- character(0)
  fit <- withCallingHandlers(
    frm(fo, data = d, family = fam, control = ctl),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      ms <<- c(ms, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  g <- drop(fit$obj$gr(fit$opt$par))
  est <- fit$opt$par
  cat(sprintf("%s | code %d | obj %.10f | par %s | grad %s\n", label,
              fit$opt$convergence, fit$opt$objective,
              paste(format(est, digits = 5), collapse = " "),
              paste(format(g, digits = 3), collapse = " ")))
  cat("   truth check:", truth(fit), "\n")
  cat("   warnings:", length(ws), "\n")
  for (w in ws) cat("     W:", substr(w, 1, 160), "\n")
  for (m in ms) cat("     M:", substr(m, 1, 160), "\n")
  se <- tryCatch(sqrt(diag(vcov(fit))), warning = function(w) NA,
                 error = function(e) NA)
  cat("   se:", paste(format(se, digits = 4), collapse = " "), "\n")
  invisible(fit)
}
set.seed(1)
d <- data.frame(g = factor(rep(c("a", "b"), each = 20)),
                y = c(rep(0, 20), stats::rpois(20, 5)))
mb <- mean(d$y[d$g == "b"])
truth <- function(fit) {
  b <- fit$opt$par
  sprintf("group b mean %.6f vs MLE %.6f (rel %.3g)", exp(sum(b)), mb,
          exp(sum(b)) / mb - 1)
}
for (k in c(2:12, 15, 20)) {
  run(sprintf("y ~ g iter.max %d", k), y ~ g, d, poisson(),
      frmtmb_control(optCtrl = list(iter.max = k, eval.max = 2 * k)), truth)
}
run("y ~ g default", y ~ g, d, poisson(), frmtmb_control(), truth)
# rel.tol loose: the optimizer's own stop, no cap
for (rt in c(1e-3, 1e-4, 1e-5)) {
  run(sprintf("y ~ g rel.tol %g", rt), y ~ g, d, poisson(),
      frmtmb_control(optCtrl = list(rel.tol = rt)), truth)
}
