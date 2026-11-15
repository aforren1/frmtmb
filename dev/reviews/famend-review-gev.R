# Review of lane famend, items 3067 and 3075: the GEV cap against Rmpfr
# and evd, the gradient band over a wider xi and x range, and the
# non-regular message on refits.
#   Rscript dev/reviews/famend-review-gev.R <base|lane>
a <- commandArgs(trailingOnly = TRUE)
lanelib <- "C:/Users/adf44/source/r/wt-famend-lib"
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] == "lane") lanelib, base, ulib))
suppressPackageStartupMessages({
  library(frmtmb)
  library(Rmpfr)
})
cat("lib", dirname(find.package("frmtmb")), "\n")
fam <- gen_extreme_value()
tape <- function(fn, y, xi) {
  p <- c(0, 0, log1p(xi))
  RTMB::MakeTape(function(p) {
    fn(y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
               .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
}
# reference log density at a double z, sigma 1, mu 0, in 200-bit floats
ref_ld <- function(z, xi) {
  zm <- mpfr(z, 200)
  if (xi == 0) return(asNumeric(-zm - exp(-zm)))
  t <- 1 + mpfr(xi, 200) * zm
  asNumeric(-(1 + 1 / mpfr(xi, 200)) * log(t) - t^(-1 / mpfr(xi, 200)))
}
ref_cdf <- function(z, xi) {
  zm <- mpfr(z, 200)
  if (xi == 0) return(asNumeric(exp(-exp(-zm))))
  t <- 1 + mpfr(xi, 200) * zm
  asNumeric(exp(-t^(-1 / mpfr(xi, 200))))
}
cat("\n== tape log density and CDF against 200-bit Rmpfr; h = log(t) / xi\n")
hs <- c(-800, -705, -690, -600, -510, -500.5, -499, -300, -50, -5, 0, 5, 30,
        300)
for (xi in c(-0.9, 0, 0.06, 0.3, 1)) {
  for (h in hs) {
    z <- if (xi == 0) h else expm1(xi * h) / xi
    if (!is.finite(z)) next
    tp <- tape(fam$lpdf, z, xi)
    v <- tp(c(0, 0, log1p(xi)))
    g <- tp$jacobian(c(0, 0, log1p(xi)))
    r <- ref_ld(z, xi)
    tc <- tape(fam$lcdf, z, xi)
    vc <- tc(c(0, 0, log1p(xi)))
    rc <- ref_cdf(z, xi)
    ev <- if (xi == 0) evd::dgev(z, 0, 1, 0, log = TRUE) else
      evd::dgev(z, 0, 1, xi, log = TRUE)
    rel <- if (is.finite(r)) abs(v - r) / max(abs(r), 1) else
      if (identical(v, r)) 0 else NA
    cat(sprintf(paste0("xi %5.2f h %8.1f z %11.4g: tape %.6g ref %.6g ",
                       "evd %.6g relerr %.3g | grad finite %s | cdf tape %.6g ",
                       "ref %.6g\n"), xi, h, z, v, r, ev, rel,
                all(is.finite(g)), vc, rc))
  }
}

cat("\n== outside rows, sigma 1: non-finite gradient / Hessian counts\n")
for (xi in c(0.001, 0.01, 0.06, 0.3, 1, 3, 10)) {
  lo <- -1 / xi
  xs <- c(lo - seq(0.01, 60, by = 0.01), lo - 10^seq(log10(60), 6, by = 0.05))
  bg <- 0L
  bh <- 0L
  nh <- 0L
  first <- NA
  for (j in seq_along(xs)) {
    tp <- tape(fam$lpdf, xs[j], xi)
    p <- c(0, 0, log1p(xi))
    gb <- any(!is.finite(tp$jacobian(p)))
    bg <- bg + gb
    if (gb && is.na(first)) first <- xs[j] - lo
    if (j %% 20 == 0) {
      nh <- nh + 1L
      bh <- bh + any(!is.finite(tp$jacfun()$jacobian(p)))
    }
  }
  cat(sprintf(paste0("xi %6.3f: x - lower end from -0.01 to -1e6, %d points: ",
                     "gradient non-finite %d (first at %s), Hessian ",
                     "non-finite %d of %d\n"),
              xi, length(xs), bg, format(first), bh, nh))
}
cat("\ndone\n")
