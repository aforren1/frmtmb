# Short check of lane famend, F1 attack: outside GEV rows whose clamped
# zc rounds ABOVE the clamp, so the ratio skips gev_ratio_top's ad_min().
# The test's setting (mu 0, sigma 1, xi through log1p), x from -19.4 to
# -20.3 by 0.001, lane gev_parts() and 0.82.0's, at the given xi.
#   Rscript dev/reviews/famend-short-gevband.R <lib> <xi> [<xi> ...]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "| BLAS",
    extSoftVersion()[["BLAS"]], "\n")
src <- readLines("dev/rel083-log/kt/famend/families-extreme-ec8178d3.R")
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
grad <- function(y, xi) {
  p <- c(0, 0, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam$lpdf(y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                     .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  list(v = tp(p), g = tp$jacobian(p))
}
# zc as the tape forms it, replayed on doubles
cl <- -1 + 1e-15
ad_pos <- function(u) 0.5 * (u + abs(u))
zc_of <- function(x, xi) {
  z <- xi * x
  zc <- z + ad_pos(cl - z)
  zc + ad_pos(cl - zc)
}
xs <- seq(-19.4, -20.3, by = -0.001)
for (xi in as.numeric(a[-1])) {
  xit <- expm1(log1p(xi))
  n_new <- n_old <- n_above <- n_above_bad <- n_at_bad <- 0
  ex <- character()
  vdiff <- 0
  for (x in xs) {
    nw <- grad(x, xi)
    o <- with_parts(old_parts, grad(x, xi))
    bn <- any(!is.finite(nw$g))
    n_new <- n_new + bn
    n_old <- n_old + any(!is.finite(o$g))
    vdiff <- vdiff + !identical(nw$v, o$v)
    zc <- zc_of(x, xit)
    if (zc > cl) {
      n_above <- n_above + 1
      n_above_bad <- n_above_bad + bn
      if (bn && length(ex) < 4) {
        ex <- c(ex, sprintf("   x %.3f zc - clamp %.3g gradient %s", x,
                            zc - cl, paste(format(nw$g, digits = 4),
                                           collapse = " ")))
      }
    } else {
      n_at_bad <- n_at_bad + bn
    }
  }
  cat(sprintf(paste0("xi %-6g rows %d | non-finite gradient: lane %d, ",
                     "0.82.0 %d | rows with zc above the clamp %d, of them ",
                     "non-finite (lane) %d; at the clamp, non-finite %d | ",
                     "value differs from 0.82.0 %d\n"),
              xi, length(xs), n_new, n_old, n_above, n_above_bad, n_at_bad,
              vdiff))
  if (length(ex)) cat(ex, sep = "\n")
}
