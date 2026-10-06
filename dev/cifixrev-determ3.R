# Reviewer: build the same objective several times in one process and
# evaluate each fresh object at one fixed point from its own initial
# inner state. Different values from identical inputs mean the tape or
# the inner solve depends on something other than the inputs.
# Usage: Rscript dev/cifixrev-determ3.R <lib or "base">
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
n <- 140L
set.seed(11)
d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
  stats::rnorm(n, 0, 0.3)
p0 <- c(-0.08, 0.78, -1.21, 1.84, -1.25, -11.4, -12.4, -0.47, 0.32)
f17 <- function(v) paste(sprintf("%.17g", v), collapse = " ")
ctl <- frmtmb_control(optCtrl = list(eval.max = 1, iter.max = 1),
                      check_se = "ignore")
for (k in 1:6) {
  fit <- suppressWarnings(
    frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)), data = d,
        control = ctl))
  obj <- fit$obj
  init <- obj$env$par
  obj$env$last.par <- init
  obj$env$last.par.best <- init
  cat("OBJ", k, "init", digest::digest(init), "\n")
  cat("FN0", k, f17(obj$fn(p0)), "\n")
  obj$env$last.par <- init
  cat("FN0b", k, f17(obj$fn(p0)), "\n")
  # the inner solution reached, and the inner Hessian at a fixed point
  cat("UHAT", k, digest::digest(obj$env$last.par), "\n")
  he <- obj$env$spHess(init, random = TRUE)
  cat("SPH", k, digest::digest(as.matrix(he)), "\n")
  cat("F0", k, f17(obj$env$f(init, order = 0)), "\n")
  cat("G0", k, digest::digest(obj$env$f(init, order = 1)), "\n")
}
