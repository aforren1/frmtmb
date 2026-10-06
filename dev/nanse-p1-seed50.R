# Punch round 1: tier 3 on mo() seed 50 (and 12), the seeds whose lost
# sets changed: spectrum, eigenvectors of the removed directions, the
# line probe's gain, and the verdicts.
#   Rscript dev/nanse-p1-seed50.R [seeds]
args <- commandArgs(TRUE)
seeds <- if (length(args)) eval(parse(text = args[1])) else c(50, 12)
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
for (s in seeds) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = d))
  nm <- ns$outer_par_names(f)
  p <- f$opt$par
  Hx <- f$obj$he(p)
  H <- (Hx + t(Hx)) / 2
  D <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(D, D), TRUE)
  cat("\n== seed", s, "\n eig:", paste(signif(e$values, 3), collapse = " "),
      "\n diag:", paste(signif(diag(H), 3), collapse = " "), "\n")
  for (k in which(e$values < 1e-3 * max(e$values))) {
    v <- e$vectors[, k]
    gain <- {
      tol <- 1e-3
      dd <- v / D * sqrt(4 * tol / abs(e$values[k]))
      f$obj$fn(p) - min(f$obj$fn(p + dd), f$obj$fn(p - dd))
    }
    cat(sprintf(" dir %.3g (probe gain %.3g):", e$values[k], gain),
        paste(sprintf("%s=%.3f", nm, v)[abs(v) > 0.01], collapse = " "),
        "\n")
  }
  cat(" lost:", paste(names(ns$sdr_of(f)$se_lost), ns$sdr_of(f)$se_lost,
                      sep = "=", collapse = " "), "\n")
}
