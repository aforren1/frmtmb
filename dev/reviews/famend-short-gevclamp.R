# Short check of lane famend, F1 attack: who reads gev_ratio_top?
#   Rscript dev/reviews/famend-short-gevclamp.R <lib>
# 1. Inside rows at xi near 0 (+-1e-8, +-1e-12, 0, +-1e-4) and rows k
#    ulps above the clamp: lane gev_parts() against 0.82.0's on the tape,
#    value, gradient and Hessian bit for bit.
# 2. Outside rows: on doubles, replay gev_parts()'s clamp arithmetic and
#    count rows whose zc is NOT the clamp (z + out rounded above it), so
#    that the ratio is read past gev_ratio_top's ad_min() with its live
#    slope; tape-evaluate those rows whose h is in (-700, -600) and count
#    non-finite gradients.
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "| BLAS",
    extSoftVersion()[["BLAS"]], "\n")
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
eval_slot <- function(slot, y, xi) {
  p <- c(mu0, ls0, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam[[slot]](y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                        .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  list(v = tp(p), g = tp$jacobian(p), h = tp$jacfun()$jacobian(p))
}
cl <- -1 + 1e-15
top <- ns$gev_ratio_top
rat <- function(z) log1p(z) / z
cat(sprintf("gev_ratio_top %.17g; ratio at clamp %.17g; (clamp - top) / top %.3g\n",
            top, rat(cl), (rat(cl) - top) / top))
for (k in 1:4) {
  zk <- cl + k * .Machine$double.eps / 2
  cat(sprintf("  zc = clamp + %d ulp: ratio %.17g, below top by %.4g\n",
              k, rat(zk), top - rat(zk)))
}
# --- 1. inside rows
same <- 0
dif <- 0
n <- 0
fin_bad <- 0
for (xi in c(-1e-4, -1e-8, -1e-12, 0, 1e-12, 1e-8, 1e-4)) {
  lo <- if (xi > 0) -1 / xi else -60
  hi <- if (xi < 0) -1 / xi else 60
  lo <- max(lo, -1e9)
  hi <- min(hi, 1e9)
  xs <- c(seq(lo, hi, length.out = 201)[-c(1, 201)], -5, -1, -1e-3, 0, 1e-3,
          1, 5, 20, -20)
  if (xi > 0) xs <- c(xs, lo + c(1e-6, 1e-3, 1) / xi)
  if (xi < 0) xs <- c(xs, hi - c(1e-6, 1e-3, 1) / abs(xi))
  xs <- xs[xi * xs > -1]
  for (slot in c("lpdf", "lcdf", "lccdf")) {
    for (x in xs) {
      y <- mu0 + exp(ls0) * x
      o <- with_parts(old_parts, eval_slot(slot, y, xi))
      nw <- eval_slot(slot, y, xi)
      n <- n + 1
      if (identical(o, nw)) same <- same + 1 else dif <- dif + 1
      fin_bad <- fin_bad + !all(is.finite(unlist(nw)))
    }
  }
}
cat(sprintf("INSIDE near xi 0: %d rows x slots, bit-identical %d, differ %d, non-finite (lane) %d\n",
            n, same, dif, fin_bad))
# rows k ulps above the clamp, z = xi * x exact by construction: choose xi
# and set y so that the double z lands on cl + k ulp
same <- 0
dif <- 0
n <- 0
zgot <- 0
for (xi in c(0.06, 0.3, 1, 3, 10, -0.5)) {
  for (k in 1:8) {
    zt <- cl + k * .Machine$double.eps / 2
    x <- zt / xi
    y <- mu0 + exp(ls0) * x
    xx <- (y - mu0) / exp(ls0)
    zz <- xi * xx
    zgot <- zgot + (zz > cl)
    if (!(zz > cl)) next
    for (slot in c("lpdf", "lcdf", "lccdf")) {
      o <- with_parts(old_parts, eval_slot(slot, y, xi))
      nw <- eval_slot(slot, y, xi)
      n <- n + 1
      if (identical(o, nw)) same <- same + 1 else dif <- dif + 1
    }
  }
}
cat(sprintf("JUST INSIDE the clamp (z above it): %d rows x slots, bit-identical %d, differ %d\n",
            n, same, dif))
# --- 2. outside rows whose zc is not the clamp after rounding
ad_pos <- function(u) 0.5 * (u + abs(u))
replay <- function(x, xi) {
  z <- xi * x
  out <- ad_pos(cl - z)
  zc <- z + out
  zc <- zc + ad_pos(cl - zc)
  list(z = z, zc = zc, h = x * min(rat(zc), top))
}
set.seed(1)
xis <- c(0.06, 0.1, 0.3, 0.5, 1, 2, 3, 5, 10, 10^stats::runif(40, -1.3, 1))
tot <- 0
above <- 0
inband <- 0
tested <- 0
bad <- 0
badlist <- character()
for (xi in xis) {
  # x where h = x * 34.5 spans (-705, -600) plus the whole outside range
  xs <- c(seq(-1 / xi - 1e-6, -1 / xi - 60, length.out = 4000),
          seq(-17.3, -20.5, length.out = 4000))
  xs <- xs[xi * xs < cl]
  for (x in xs) {
    y <- mu0 + exp(ls0) * x
    xx <- (y - mu0) / exp(ls0)
    r <- replay(xx, expm1(log1p(xi)))
    tot <- tot + 1
    if (r$zc > cl) {
      above <- above + 1
      if (r$h > -700 && r$h < -600) {
        inband <- inband + 1
        if (tested < 400) {
          tested <- tested + 1
          nw <- eval_slot("lpdf", y, xi)
          if (any(!is.finite(nw$g))) {
            bad <- bad + 1
            badlist <- c(badlist, sprintf("xi %.6g x %.6g zc-cl %.3g h %.2f",
                                          xi, xx, r$zc - cl, r$h))
          }
        }
      }
    }
  }
}
cat(sprintf(paste0("OUTSIDE: %d rows over %d xi; zc above the clamp after ",
                   "rounding %d; of those with h in (-700, -600) %d; ",
                   "tape-tested %d, non-finite gradient %d\n"),
            tot, length(xis), above, inband, tested, bad))
if (length(badlist)) cat(head(badlist, 10), sep = "\n")
# direct construction: z just below the clamp by a fraction of an ulp of
# d = cl - z, so z + fl(cl - z) rounds above cl, at x where h is near -690
tested2 <- 0
above2 <- 0
bad2 <- 0
for (xi in c(0.5, 1, 2, 3, 5, 10, 20, 34)) {
  for (x in seq(-19.4, -20.3, length.out = 3000)) {
    xi <- expm1(log1p(xi))
    y <- mu0 + exp(ls0) * x
    r <- replay((y - mu0) / exp(ls0), xi)
    if (r$zc > cl && r$h > -700) {
      above2 <- above2 + 1
      if (tested2 < 300) {
        tested2 <- tested2 + 1
        nw <- eval_slot("lpdf", y, xi)
        bad2 <- bad2 + any(!is.finite(nw$g))
      }
    }
  }
}
cat(sprintf("BAND x in [-20.3, -19.4], xi 0.5 to 34: zc above clamp with h > -700: %d; tape-tested %d, non-finite gradient %d\n",
            above2, tested2, bad2))
