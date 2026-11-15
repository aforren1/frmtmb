# Re-check of lane famend after fix round 2: the 4 tapes of
# famend-short-subnorm.R that differ (weight on h only, z = +-0.001):
# which part differs, and by how much.
#   Rscript dev/reviews/famend-short-subnorm2.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
src <- readLines("dev/reviews/famend-short-out/fe-base.R")
i0 <- grep("^gev_parts <- function", src)
i1 <- i0 - 1L + grep("^}", src[i0:length(src)])[1]
old_parts <- eval(parse(text = src[i0:i1]))
environment(old_parts) <- ns
run <- function(parts, q, p, w) {
  tp <- RTMB::MakeTape(function(p) {
    g <- parts(q, list(mu = p[1], sigma = exp(p[2]), xi = p[3]))
    w[1] * g$l + w[2] * g$h
  }, p)
  list(v = tp(p), g = tp$jacobian(p), h = tp$jacfun()$jacobian(p))
}
for (cs in list(c(-0.9, 0.001), c(0.3, 0.001), c(0.6, 0.001), c(10, -0.001),
                c(0.3, 0.0011), c(0.3, 0.005), c(0.3, 0.01), c(0.3, 0.02))) {
  xi <- cs[1]
  x <- cs[2] / xi
  q <- 0.3 + exp(-0.2) * x
  p <- c(0.3, -0.2, xi)
  for (w in list(c(0, 1), c(1, 0))) {
    o <- run(old_parts, q, p, w)
    n <- run(ns$gev_parts, q, p, w)
    rd <- function(u, v) max(c(0, abs(u - v)[u != v] / pmax(abs(u[u != v]),
                                                           1e-300)))
    cat(sprintf("xi %5g z %7g w %s: v same %s, g same %s (max rel %.3g), H same %s (max rel %.3g, entries differing %d of 9)\n",
                xi, cs[2], paste(w, collapse = ","), identical(o$v, n$v),
                identical(o$g, n$g), rd(o$g, n$g), identical(o$h, n$h),
                rd(o$h, n$h), sum(o$h != n$h)))
  }
}
