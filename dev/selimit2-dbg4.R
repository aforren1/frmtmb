# Lane selimit2: boundfam's separating zi and hu slopes on the lane build.
#   Rscript dev/selimit2-dbg4.R <lib or base>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-selimit2/dev/selimit2-probe.R")
eval(parse(text = src[grep("^cap <- ", src):(grep("^if [(]want", src)[1] - 1L)]))
for (s in 1:2) {
  set.seed(s)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- ifelse(d$x < 0, 0L, stats::rpois(200, exp(2)))
  f <- report(paste("zisep", s), cap(frm(bf(y ~ 1, zi ~ x),
                                         family = zero_inflated_poisson(),
                                         data = d)))
  str(ns$se_limit_info(f)[c("kind", "dir", "sep", "sep1", "joint", "gain")])
}
