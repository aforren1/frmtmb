# Re-check of lane famend after fix round 2 (xi on the identity, not log1p): the three slots on
# inside rows with |z| below 0.016 (log1p_ratio()'s blend), xi through
# log1p as in dev/famend-gevbit.R, lane gev_parts() against 0.82.0's:
# value, gradient and Hessian bit for bit, and the largest Hessian
# difference where they are not.
#   Rscript dev/reviews/famend-short-blend.R <lib>
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
fam <- gen_extreme_value()
with_parts <- function(parts, expr) {
  orig <- get("gev_parts", ns)
  unlockBinding("gev_parts", ns)
  assign("gev_parts", parts, ns)
  on.exit({ assign("gev_parts", orig, ns); lockBinding("gev_parts", ns) })
  expr
}
ev <- function(slot, y, xi) {
  p <- c(0.3, -0.2, xi)
  tp <- RTMB::MakeTape(function(p) {
    fam[[slot]](y, list(mu = p[1], sigma = exp(p[2]), xi = p[3] + 0,
                        .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  list(v = tp(p), g = tp$jacobian(p), h = tp$jacfun()$jacobian(p))
}
n <- 0; vs <- 0; gs <- 0; hs <- 0; hmax <- 0
for (xi in c(-0.9, -0.5, -0.1, 0.001, 0.01, 0.06, 0.3, 0.6, 1, 3, 10)) {
  for (z in c(-0.0159, -0.012, -0.008, -0.004, -0.001, -1e-5, 1e-5, 0.001,
              0.004, 0.008, 0.012, 0.0159)) {
    y <- 0.3 + exp(-0.2) * (z / xi)
    for (slot in c("lpdf", "lcdf", "lccdf")) {
      o <- with_parts(old_parts, ev(slot, y, xi))
      nw <- ev(slot, y, xi)
      n <- n + 1
      vs <- vs + identical(o$v, nw$v)
      gs <- gs + identical(o$g, nw$g)
      hs <- hs + identical(o$h, nw$h)
      d <- o$h != nw$h
      if (any(d)) hmax <- max(hmax, abs(o$h - nw$h)[d] / abs(o$h)[d])
    }
  }
}
cat(sprintf("BLEND |z| < 0.016: %d rows x slots; value identical %d, gradient identical %d, Hessian identical %d; largest Hessian relative difference %.3g\n",
            n, vs, gs, hs, hmax))
