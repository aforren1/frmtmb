# Review of lane quadseam: what lane boundfam's coefficient-limit verdict
# (se_edge_limit(), gated off under quadrature) leaves on the one-block
# quadrature route. A poisson with a factor level whose every count is
# zero (its coefficient runs off) and (1 | g), Laplace and quadrature.
#   Rscript quadseam-review-limit.R <libpath>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(5)
d <- data.frame(g = factor(rep(1:20, each = 9)),
                f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
d$y <- rpois(180, exp(0.5 + 0.3 * d$x + rnorm(20, 0, 0.5)[d$g]))
d$y[d$f == "c"] <- 0L
for (quad in c(FALSE, TRUE)) {
  ms <- character(); w <- character()
  f <- withCallingHandlers(
    frm(y ~ f + x + (1 | g), family = poisson(), data = d,
        quadrature = quad),
    warning = function(x) {
      w <<- c(w, paste(class(x)[1], substr(conditionMessage(x), 1, 160)))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      ms <<- c(ms, paste(class(x)[1], substr(conditionMessage(x), 1, 160)))
      invokeRestart("muffleMessage")
    })
  v <- tryCatch(sqrt(diag(vcov(f))), error = function(e) NA)
  cat(sprintf("QUAD %s fc %.4g se(fc) %s conv %d\n", quad,
              f$estimates$beta[["fc"]], format(v[["fc"]], digits = 4),
              f$opt$convergence))
  for (x in c(w, ms)) cat("   ", gsub("\\s+", " ", x), "\n")
}
cat("END\n")
