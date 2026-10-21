# Release reviewer, 0.73.0: the asym_laplace observation-level refusal
# against cdfs's new censored families: zero_inflated_asym_laplace, and
# asym_laplace under cens() with no row censored.
#   Rscript dev/rel073-rev-alap.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")), "\n")
run <- function(label, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) cat(label, ": ERROR:", substr(gsub("[[:space:]]+", " ", r), 1, 150), "\n")
  else cat(label, sprintf(": logLik %.4f conv %d sigma %.3g", as.numeric(logLik(r)),
                          r$opt$convergence, exp(r$estimates$betad[1])),
           "| warnings", length(w), "\n")
}
set.seed(7)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)), cn = "none")
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$yz <- ifelse(runif(n) < 0.2, 0, d$y)
run("asym_laplace (1 | obs)", frm(y ~ x + (1 | obs), family = asym_laplace(), data = d))
run("asym_laplace cens(none) (1 | obs)", frm(y | cens(cn) ~ x + (1 | obs),
    family = asym_laplace(), data = d))
run("zero_inflated_asym_laplace (1 | obs)", frm(yz ~ x + (1 | obs),
    family = zero_inflated_asym_laplace(), data = d))
run("zero_inflated_asym_laplace, no obs term", frm(yz ~ x,
    family = zero_inflated_asym_laplace(), data = d))
