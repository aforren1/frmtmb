# Defect 8 cause probe: `ls ~ mo(income) * age` on brms_monotonic's
# data code. Per seed: the SEs, the outer Hessian's eigenvalues, the
# parameters loading on the smallest eigenvector, the simplex values.
#
#   Rscript dev/nanse-mo-probe.R [lib] [seeds]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:40
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
for (s in seeds) {
  d <- mk(s)
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = d))
  sdr <- frmtmb:::sdr_of(f)
  se <- sqrt(diag(sdr$cov.fixed))
  p <- f$opt$par
  H <- optimHess(p, f$obj$fn, f$obj$gr)
  Hr <- tryCatch(f$obj$he(p), error = function(e) NULL)
  ev <- eigen(H, symmetric = TRUE)
  nm <- frmtmb:::outer_par_names(f)
  v <- ev$vectors[, length(p)]
  cat(sprintf("\n== seed %d code %d allNaN %s pdHess %s units %s\n", s,
              f$opt$convergence, all(!is.finite(se)), sdr$pdHess,
              paste(unique(signif(f$par_units %||% 1, 3)), collapse = ",")))
  cat("par:", paste(sprintf("%s=%.4g", nm, p), collapse = " "), "\n")
  cat("eig(optimHess):", paste(signif(ev$values, 3), collapse = " "), "\n")
  if (!is.null(Hr)) {
    cat("eig(obj$he):", paste(signif(eigen(Hr, TRUE, only.values = TRUE)$values,
                                      3), collapse = " "), "\n")
  }
  cat("null dir:", paste(sprintf("%s=%.3f", nm, v)[abs(v) > 0.05],
                         collapse = " "), "\n")
  cat("sdr cov finite:", sum(is.finite(sdr$cov.fixed)), "of",
      length(sdr$cov.fixed), "; rcond(H)", signif(rcond(H), 3), "\n")
  zt <- f$estimates
  print(lapply(zt[grepl("simo|zeta", names(zt))], signif, 4))
}
