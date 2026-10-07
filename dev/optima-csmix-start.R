# Lane optima, item 3: the objective and gradient of a cs() ordinal
# mixture at its starting values, on base and lane, for one seed.
#   Rscript dev/optima-csmix-start.R base|lane seed
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seed <- as.integer(args[2])
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
d <- mk(seed, mix = TRUE)
u <- suppressWarnings(frm(bf(y ~ x + cs(z)),
                          family = mixture(cumulative(), sratio()),
                          data = d, dry_run = "objective"))
p <- u$obj$par
cat("fn", format(u$obj$fn(p), digits = 12), "\n")
g <- u$obj$gr(p)
print(stats::setNames(as.numeric(g), names(p)))
set.seed(1)
for (i in 1:5) {
  q <- p + stats::rnorm(length(p), 0, 0.3)
  cat("jitter", i, "fn", format(u$obj$fn(q), digits = 10), "grad finite",
      all(is.finite(u$obj$gr(q))), "\n")
}
