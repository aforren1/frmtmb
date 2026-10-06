# Reviewer: is the objective itself deterministic across R processes?
# Fits the n = 140, seed 11 data set of dev/cifix-scan2.R, then
# evaluates the objective and gradient at two fixed points, printing 17
# digits. Usage: Rscript dev/cifixrev-determ2.R <lib or "base">
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
n <- 140L
s <- 11L
set.seed(s)
d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
  stats::rnorm(n, 0, 0.3)
fit <- suppressWarnings(
  frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)), data = d))
obj <- fit$obj
p0 <- c(-0.08, 0.78, -1.21, 1.84, -1.25, -11.4, -12.4, -0.47, 0.32)
p1 <- c(-0.08, 0.78, -1.21, 1.84, -1.25, -3, -3, -0.47, 0.32)
f17 <- function(v) paste(sprintf("%.17g", v), collapse = " ")
cat("FIT", f17(fit$opt$par), "\n")
init <- obj$env$par
cat("INIT", digest::digest(init), digest::digest(obj$env$data), "
")
for (k in 1:2) {
  obj$env$last.par <- init; obj$env$last.par.best <- init
  cat("FN0", f17(obj$fn(p0)), "\n")
  cat("GR0", f17(obj$gr(p0)), "\n")
  cat("FN1", f17(obj$fn(p1)), "\n")
  cat("GR1", f17(obj$gr(p1)), "\n")
}
# the inner Hessian at a fixed point, to see whether the inner solve's
# input differs
cat("RAND", digest::digest(obj$env$last.par.best), "\n")
