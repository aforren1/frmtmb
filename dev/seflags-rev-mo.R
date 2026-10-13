# Reviewer of lane seflags: moparity's new mo() structures through the
# new estimability rule. mo(x) * f and (1 + mo(x) | g), every category
# observed and the third of four unobserved. Lost set, warnings,
# fitted(newdata) SEs and their warnings, conditional_effects() bands.
#   Rscript dev/seflags-rev-mo.R <merge|moparity>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(switch(a[1],
  merge = "C:/Users/adf44/source/r/seflags-rev-lib/mo",
  moparity = "C:/Users/adf44/source/r/wt-moparity-lib"),
  "C:/Users/adf44/source/r/rellib-r7",
  "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", a[1], find.package("frmtmb"), "\n")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e)
    structure(conditionMessage(e), class = "rev_err")),
    warning = function(x) { w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning") },
    message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
for (seed in 1:4) for (miss in c(FALSE, TRUE)) {
  set.seed(seed)
  n <- 160
  cats <- if (miss) c(0L, 1L, 3L) else 0:3
  x <- sample(cats, n, TRUE)
  f <- factor(sample(c("a", "b"), n, TRUE))
  g <- factor(rep(1:8, length.out = n))
  eff <- c(0, 2, 3, 3.5)
  y <- 1 + eff[x + 1] * ifelse(f == "a", 1, 1.5) +
    stats::rnorm(8, 0, 0.5)[g] + rnorm(8, 0, 0.3)[g] * eff[x + 1] +
    stats::rnorm(n, 0, 0.7)
  d <- data.frame(y, x, f, g)
  for (fo in list(y ~ mo(x) * f, y ~ mo(x) + (1 + mo(x) | g))) {
    r <- cap(frm(fo, data = d))
    fit <- r$v
    lab <- sprintf("seed %d miss %-5s %-26s", seed, miss, deparse(fo))
    if (inherits(fit, "rev_err")) { cat(lab, "ERROR", fit, "\n"); next }
    lost <- suppressWarnings(ns$sdr_of(fit)$se_lost)
    nd <- d[c(1, 2, 3, 4), ]
    p <- cap(fitted(fit, newdata = nd, re_formula = NA))
    ce <- cap(conditional_effects(fit, effects = "x"))
    cen <- if (inherits(ce$v, "rev_err")) NA else
      sum(!is.finite(ce$v[[1]]$se__))
    cat(sprintf("%s code %d | lost %s | fit warn %d | fitted se %s warn %d | ce NaN se %s warn %d\n",
                lab, fit$opt$convergence,
                if (length(lost)) paste(names(lost), lost, sep = ":",
                                        collapse = ",") else "-",
                length(r$w),
                if (inherits(p$v, "rev_err")) p$v else
                  paste(signif(p$v[, 2], 3), collapse = " "),
                length(p$w), cen, length(ce$w)))
  }
}
