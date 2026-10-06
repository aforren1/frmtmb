# Is the non-PD scaled Hessian on saturated-simplex seeds real or the
# finite-difference step? Compares optimHess (what sdreport uses) with
# the exact AD Hessian obj$he() on given seeds.
#   Rscript dev/nanse-mo-hess.R [lib] [seeds]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else c(12, 18, 7)
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
sc <- function(H) {
  D <- sqrt(abs(diag(H)))
  ok <- D > 0
  S <- (H / outer(D, D))[ok, ok]
  signif(eigen(S, TRUE, only.values = TRUE)$values, 3)
}
for (s in seeds) {
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = mk(s)))
  p <- f$opt$par
  nm <- frmtmb:::outer_par_names(f)
  H1 <- optimHess(p, f$obj$fn, f$obj$gr)
  H2 <- f$obj$he(p)
  cat("\n== seed", s, "\n")
  cat("grad:", paste(sprintf("%s=%.2g", nm, f$obj$gr(p)), collapse = " "),
      "\n")
  cat("diag optimHess:", paste(signif(diag(H1), 3), collapse = " "), "\n")
  cat("diag he       :", paste(signif(diag(H2), 3), collapse = " "), "\n")
  cat("scaled eig optimHess:", paste(sc(H1), collapse = " "), "\n")
  cat("scaled eig he       :", paste(sc(H2), collapse = " "), "\n")
  D <- sqrt(abs(diag(H2)))
  S <- H2 / outer(D, D)
  dimnames(S) <- list(nm, nm)
  print(round(S, 3))
}
