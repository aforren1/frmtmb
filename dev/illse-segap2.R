# Lane illse, item b: what the SE check's ridge rule reads on a silent
# fit of dev/illse-segap.R (base build): every se_ridge_curvature() and
# se_curvature_real() call with its result.
#   Rscript dev/illse-segap2.R <seed> <formula 1|2|4> [base|dev]
args <- commandArgs(TRUE)
seed <- as.integer(args[1]); which_fo <- as.integer(args[2])
arm <- if (length(args) > 2) args[3] else "base"
R <- "C:/Users/adf44/source/r"
.libPaths(c(file.path(R, "rellib-r8"), "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "dev") {
  suppressMessages(pkgload::load_all(file.path(R, "frmtmb-wt-illse"), quiet = TRUE))
} else suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
se_mo_gap <- function(seed) {
  set.seed(seed)
  inc <- sample(c(0L, 1L, 3L), 100, TRUE)
  ls <- c(30, 60, 70, 75)[inc + 1L] + rnorm(100, sd = 7)
  data.frame(inc, ls, age = rnorm(100, mean = 40, sd = 10))
}
src <- readLines(file.path(R, "frmtmb-wt-illse/dev/relrev070-interplay.R"))
i <- grep("^mk <- function", src); j <- grep("^designs <- list", src)
eval(parse(text = src[i:(j - 1L)]))
fo <- switch(which_fo, ls ~ mo(inc), ls ~ mo(inc) + age, ls ~ mo(inc) * age, y ~ mo(x2))
d <- if (which_fo == 4L) mk(seed) else se_mo_gap(seed)
suppressMessages({
trace("se_ridge_curvature", where = ns, print = FALSE, exit = quote(
  cat(sprintf("ridge: k %d ev %.3e c0 %.3e thr_flat %.3e\n", k, e$values[k],
              returnValue(), ns$se_flat_tol * max(abs(e$values))))))
trace("se_curvature_real", where = ns, print = FALSE, exit = quote(
  cat(sprintf("curvature_real: lambda %.3e -> %s\n", lambda, returnValue()))))
trace("se_tier3", where = ns, print = FALSE, tracer = quote(cat("tier3 entered\n")))
})
f <- suppressWarnings(frm(fo, data = d))
cat("lost:", names(ns$sdr_of(f)$se_lost), "\n")
