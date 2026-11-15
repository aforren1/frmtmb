# Review of lane slopescale, attack item 7: frmtmb.sample on an
# autoscaled random-slope fit. frm_sample(laplace = TRUE) on the fit of
# y ~ 1 + (1 + x | g) with x * 1e-4 (the default engages) and on the same
# data with x (it does not). The slope's sd draws, times 1e-4, and the
# other outer draws must describe the same posterior if the reported fit
# is on the original scale; with a flat prior on the log sd the sampler
# is translation invariant there, so the summaries should agree closely.
# Usage: Rscript slopescale-review-sample.R
.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", find.package("frmtmb"), "\nfrmtmb.sample",
    find.package("frmtmb.sample"), "\n")
S <- 1e-4
set.seed(2)
g <- factor(rep(1:20, each = 15))
x0 <- rnorm(300)
y <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * x0 + rnorm(300)
d1 <- data.frame(g = g, x = x0, y = y)
ds <- d1
ds$x <- x0 * S
f1 <- frm(y ~ 1 + (1 + x | g), data = d1)
fs <- frm(y ~ 1 + (1 + x | g), data = ds)
cat("autoscaled: spread 1", isTRUE(f1$autoscaled), "; spread 1e-4",
    isTRUE(fs$autoscaled), "\n")
cat("logLik diff", as.numeric(logLik(f1)) - as.numeric(logLik(fs)), "\n")
smp <- function(f) {
  suppressMessages(suppressWarnings(
    frm_sample(f, chains = 2, iter = 600, refresh = 0, seed = 11,
               laplace = TRUE)))
}
s1 <- smp(f1)
ss <- smp(fs)
m1 <- s1$draws
ms <- ss$draws
cat("columns:", paste(colnames(m1), collapse = " "), "\n")
sc <- ifelse(grepl("__x$", colnames(ms)), S, 1)
ms2 <- sweep(ms, 2, sc, "*")
qs <- function(m) apply(m, 2, stats::quantile, c(0.05, 0.5, 0.95))
q1 <- qs(m1)
q2 <- qs(ms2)
print(signif(q1, 4))
print(signif(q2, 4))
keep <- setdiff(colnames(m1), "lp__")
cat("max relative difference of the 5/50/95 % quantiles (lp__ out):",
    signif(max(abs(q2[, keep] / q1[, keep] - 1)), 3), "\n")
cat("identical draws after the scale:", isTRUE(all.equal(ms2[, keep],
                                                        m1[, keep],
                                                        tolerance = 1e-6)),
    "\n")
