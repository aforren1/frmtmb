# Short check of lane squant fix round 1, item 5: the standard-error
# routes ?quantile_elf documents ("Standard errors", "Group-level terms",
# the example) run and say what the page says. Split out of
# squant-short-compat.R, whose gp() fit did not finish.
# Run (worktree root):
#   Rscript dev/reviews/squant-short-se.R > dev/reviews/squant-short-se.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
run <- function(lab, expr) {
  ws <- character()
  ms <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, substr(conditionMessage(w), 1, 150))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    ms <<- c(ms, substr(conditionMessage(m), 1, 150))
    invokeRestart("muffleMessage")
  }), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%s: ERROR %s\n", lab, substr(conditionMessage(r), 1, 300)))
    return(invisible(NULL))
  }
  cat(sprintf("%s: conv %s | logLik %.4f | warnings %d | messages %d\n",
              lab, r$opt$convergence, as.numeric(logLik(r)), length(ws),
              length(ms)))
  if (length(ws)) cat("   W:", paste(unique(ws), collapse = "\n   W: "), "\n")
  if (length(ms)) cat("   M:", paste(unique(ms), collapse = "\n   M: "), "\n")
  invisible(r)
}
# the page's example
set.seed(1)
nn <- 300
dd <- data.frame(x = runif(nn), row = seq_len(nn))
dd$y <- 1 + 2 * dd$x + (0.5 + dd$x) * rnorm(nn)
fit <- run("example fit", frm(bf(y ~ x, quantile = 0.9),
                              family = quantile_elf(), data = dd))
V <- vcov_cluster(fit, ~ row, "CR0", full = TRUE)
cat("default summary:\n")
print(summary(fit)$fixed)
cat("summary(fit, vcov = V):\n")
print(summary(fit, vcov = V)$fixed)
cat("confint(fit, vcov = V):\n")
print(confint(fit, vcov = V))
vc <- vcov(fit, cluster = ~ row)
Vd <- vcov_cluster(fit, ~ row)
k <- seq_len(nrow(vc))
cat("vcov(fit, cluster = ~ row) SEs:", sprintf("%.6f", sqrt(diag(vc))), "\n")
cat("vcov_cluster(fit, ~ row) default type, same rows:",
    sprintf("%.6f", sqrt(diag(Vd))[k]), "\n")
cat("vcov_cluster(fit, ~ row, 'CR0', full = TRUE), same rows:",
    sprintf("%.6f", sqrt(diag(V))[k]), "\n")
cat("Hessian SEs:", sprintf("%.6f", sqrt(diag(vcov(fit)))), "\n")

set.seed(3)
dg <- data.frame(x = runif(120), g = factor(rep(1:30, each = 4)))
dg$y <- 1 + dg$x + rnorm(30, 0, 0.5)[dg$g] + rnorm(120)
fq <- run("quadrature (1 | g), smoothing 1",
          frm(bf(y ~ x + (1 | g), quantile = 0.5),
              family = quantile_elf(smoothing = 1), data = dg,
              quadrature = TRUE))
if (!is.null(fq)) {
  r <- tryCatch(vcov_cluster(fq, ~ g, "CR0", full = TRUE),
                error = function(e) e)
  cat("vcov_cluster on the quadrature fit:",
      if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
        "returned a matrix", "\n")
  cat("quadrature fit Hessian SEs:",
      sprintf("%.4f", sqrt(diag(vcov(fq)))), "\n")
}
fl <- run("Laplace (1 | g), smoothing 1",
          frm(bf(y ~ x + (1 | g), quantile = 0.5),
              family = quantile_elf(smoothing = 1), data = dg))
if (!is.null(fl)) {
  r <- tryCatch(vcov_cluster(fl, ~ g, "CR0", full = TRUE),
                error = function(e) e)
  cat("vcov_cluster(fit, ~ g, 'CR0', full = TRUE) on the Laplace fit:",
      if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
        sprintf("SEs %s", paste(sprintf("%.4f", sqrt(diag(r))[1:2]),
                                collapse = " ")), "\n")
  cat("Laplace fit Hessian SEs:", sprintf("%.4f", sqrt(diag(vcov(fl)))),
      "\n")
}
cat("DONE\n")
