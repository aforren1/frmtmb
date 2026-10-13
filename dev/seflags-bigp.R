# Lane seflags: the cost of the SE check on dev/setier-rev2-sepmiss.R's
# 2050-level sparse bernoulli fit (seed 1): the fit with and without the
# check, and how many directions se_ridge_curvature() is asked of.
#   Rscript dev/seflags-bigp.R <lane|base>
a <- commandArgs(trailingOnly = TRUE)
libs <- switch(a[1],
  lane = c("C:/Users/adf44/source/r/wt-seflags-lib",
           "C:/Users/adf44/source/r/rellib-r7"),
  base = "C:/Users/adf44/source/r/rellib-r7")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
set.seed(1)
L <- 2050
d <- data.frame(f = factor(rep(seq_len(L), each = 4)), x = rnorm(4 * L))
d$y <- rbinom(4 * L, 1, plogis(0.5 + 0.5 * d$x))
n_ridge <- 0L
if (exists("se_ridge_curvature", ns)) {
  trace("se_ridge_curvature", quote(n_ridge <<- n_ridge + 1L),
        print = FALSE, where = ns)
}
for (ch in c("ignore", "warning", "ignore", "warning")) {
  t <- system.time(f <- suppressWarnings(suppressMessages(
    frm(y ~ f + x, family = bernoulli(), data = d,
        control = frmtmb_control(sparse_x = TRUE, check_se = ch)))))
  cat(sprintf("check_se %-7s fit %.1f s | evals %s | lost %d\n", ch,
              t[["elapsed"]], format(f$opt$evals),
              length(suppressWarnings(ns$sdr_of(f))$se_lost)))
}
cat("se_ridge_curvature calls:", n_ridge, "\n")
