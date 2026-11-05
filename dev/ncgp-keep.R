# m1: how many points the marginal should be polished from, and what each one costs.
# keep = 1, 2, 3 on the review's stationary designs and curve set, and at n 500 for the
# time. Reads the marginal's own value, so no tape is built past the warm start.
#   Rscript dev/ncgp-keep.R
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
frame_of <- function(d) suppressWarnings(frm(y ~ gp(x), d))$frame
read <- function(tag, d) {
  fr <- frame_of(d)
  s <- ns$gp_marginal_setup(fr)
  v <- numeric(3); t <- numeric(3)
  for (k in 1:3) {
    t0 <- proc.time()[[3]]
    v[k] <- -ns$gp_marginal_opt(s, k)$value
    t[k] <- proc.time()[[3]] - t0
  }
  cat(sprintf("KEEP %-24s ll %.5f | keep 1 short %.2e %.2fs | keep 2 short %.2e %.2fs | keep 3 %.2fs\n",
              tag, v[3], v[3] - v[1], t[1], v[3] - v[2], t[2], t[3]))
}
for (seed in 1:3) {
  set.seed(seed); x <- sort(runif(80, 0, 10))
  read(sprintf("log1p(x) seed %d", seed), data.frame(x = x, y = log1p(x) + 1e-3 * rnorm(80)))
}
set.seed(131); u <- sort(runif(20, 0, 6)); x <- rep(u, each = 30); e <- rnorm(600)
read("rep 20x30 s=1e-2 seed 131", data.frame(x = x, y = sin(2 * x) + 1e-2 * e))
fns <- list(mix = function(x) 2 * sin(3 * x) + x, mix2 = function(x) sin(4 * x) + 0.5 * x,
            mix3 = function(x) 3 * cos(2 * x) - 2 * x, mix4 = function(x) sin(5 * x) + x^2 / 10,
            two = function(x) sin(x) + 0.3 * sin(10 * x), sin3 = function(x) sin(3 * x),
            sin10 = function(x) sin(10 * x))
for (seed in 1:2) for (fn in names(fns)) for (sg in c(0.3, 1e-2)) {
  set.seed(31000 + seed * 7 + 120L)
  x <- sort(round(runif(120, 0, 10), 4))
  read(sprintf("%s s=%g seed %d", fn, sg, seed), data.frame(x = x, y = fns[[fn]](x) + sg * rnorm(120)))
}
set.seed(20261103)
x <- sort(runif(500, 0, 10))
read("sin(x) n 500", data.frame(x = x, y = sin(x) + 0.3 * rnorm(500)))
