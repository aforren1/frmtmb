# Lane selimit2: zi ~ 0 + f and zi ~ f with level a all zero, on the
# test's data (no x column drawn first).
#   Rscript dev/selimit2-dbg3.R <lib or base>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-selimit2/dev/selimit2-probe.R")
eval(parse(text = src[grep("^cap <- ", src):(grep("^if [(]want", src)[1] - 1L)]))
for (seed in c(21, 22, 23)) {
  set.seed(seed)
  n <- 300
  d <- data.frame(f = factor(rep(c("a", "b"), each = n / 2)))
  d$y <- ifelse(d$f == "a", 0L, stats::rpois(n, 4))
  report(paste("seed", seed, "zi ~ 0 + f"),
         cap(frm(bf(y ~ 1, zi ~ 0 + f), family = zero_inflated_poisson(),
                 data = d)))
  f <- report(paste("seed", seed, "zi ~ f"),
              cap(frm(bf(y ~ 1, zi ~ f), family = zero_inflated_poisson(),
                      data = d)))
  str(ns$se_limit_info(f)[c("kind", "dir", "sep", "sep1", "joint", "gain")])
  report(paste("seed", seed, "hu ~ f"),
         cap(frm(bf(y ~ 1, hu ~ f), family = hurdle_poisson(), data = d)))
}
