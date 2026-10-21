# Reviewer: asym_laplace with an observation-level effect. Does the
# Laplace objective run away without weights too, and does a one-block
# quadrature fit (exact marginal) stay bounded?
#   Rscript dev/fuzzfix-rev-asym.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261020)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                g = factor(rep(1:20, each = 4)))
d$y <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5) + rnorm(n, 0, 0.7)
d$w9 <- 0.9
run <- function(lbl, expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-30s ERROR %s\n", lbl, substr(conditionMessage(r), 1, 140)))
  } else {
    cat(sprintf("%-30s logLik %.6g betad %s warnings %d %s\n", lbl,
                as.numeric(logLik(r)),
                paste(sprintf("%.4g", r$estimates$betad), collapse = " "),
                length(ws), if (length(ws)) substr(ws[1], 1, 70) else ""))
  }
}
AL <- asym_laplace()
run("obs, no weights, Laplace", frm(bf(y ~ x + (1 | obs)), family = AL, data = d))
run("obs, w = 0.9, Laplace", frm(bf(y | weights(w9) ~ x + (1 | obs)),
                                 family = AL, data = d))
run("obs, no weights, quadrature", frm(bf(y ~ x + (1 | obs)), family = AL,
                                       data = d, quadrature = TRUE))
run("group of 4, no weights, Laplace", frm(bf(y ~ x + (1 | g)), family = AL,
                                           data = d))
run("no effect", frm(bf(y ~ x), family = AL, data = d))
