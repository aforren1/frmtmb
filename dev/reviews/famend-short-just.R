# Re-check of lane famend after fix round 2: the "JUST INSIDE the clamp"
# block of famend-short-gevclamp.R (138 rows x slots, 14 differ). Classify
# each row as the tape forms it (x = (y - mu) / sigma, xi = expm1(log1p(xi)),
# out = gev_parts()'s own out on doubles) and count identical rows per class.
#   Rscript dev/reviews/famend-short-just.R <lib>
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
mu0 <- 0.3
ls0 <- -0.2
ev <- function(slot, y, xi) {
  p <- c(mu0, ls0, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam[[slot]](y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                        .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  list(v = tp(p), g = tp$jacobian(p), h = tp$jacfun()$jacobian(p))
}
cl <- -1 + 1e-15
tab <- matrix(0, 2, 2, dimnames = list(c("tape inside", "tape outside"),
                                       c("identical", "differ")))
for (xi in c(0.06, 0.3, 1, 3, 10, -0.5)) {
  for (k in 1:8) {
    zt <- cl + k * .Machine$double.eps / 2
    x <- zt / xi
    y <- mu0 + exp(ls0) * x
    if (!(xi * ((y - mu0) / exp(ls0)) > cl)) next
    # z as the tape itself forms it
    p0 <- c(mu0, ls0, log1p(xi))
    z_tape <- RTMB::MakeTape(function(p) {
      expm1(p[3]) * ((y - p[1]) / exp(p[2]))
    }, p0)(p0)
    cls <- if (z_tape <= cl) "tape outside" else "tape inside"
    for (slot in c("lpdf", "lcdf", "lccdf")) {
      same <- identical(ev(slot, y, xi), with_parts(old_parts, ev(slot, y, xi)))
      tab[cls, if (same) "identical" else "differ"] <-
        tab[cls, if (same) "identical" else "differ"] + 1
    }
  }
}
print(tab)
