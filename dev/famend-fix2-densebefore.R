# Lane famend, fix round 2: the dense xi grid of dev/famend-gevbit.R
# (0.4 to 0.9 by 0.01, same rows) with the gev_parts() of 0.82.0 and of
# the first fix round (0.82.0's clamp, ratio held by gev_ratio_top),
# swapped into this build: the count of non-finite gradients, before.
#   Rscript dev/famend-fix2-densebefore.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
src <- readLines("dev/reviews/famend-short-out/fe-base.R")
i0 <- grep("^gev_parts <- function", src)
i1 <- i0 - 1L + grep("^}", src[i0:length(src)])[1]
base_parts <- eval(parse(text = src[i0:i1]))
environment(base_parts) <- ns
round1_parts <- function(q, dpars) {
  x <- (q - dpars[["mu"]]) / dpars[["sigma"]]
  z <- dpars[["xi"]] * x
  out <- ad_pos(-1 + 1e-15 - z)
  zc <- ad_max(z + out, -1 + 1e-15)
  list(l = log1p(zc), h = x * ad_min(log1p_ratio(zc), gev_ratio_top),
       out = out)
}
environment(round1_parts) <- ns
fam <- gen_extreme_value()
with_parts <- function(parts, expr) {
  orig <- get("gev_parts", ns)
  unlockBinding("gev_parts", ns)
  assign("gev_parts", parts, ns)
  on.exit({ assign("gev_parts", orig, ns); lockBinding("gev_parts", ns) })
  expr
}
grad <- function(y, xi) {
  p <- c(0.3, -0.2, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam$lpdf(y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                     .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  tp$jacobian(p)
}
for (nm in c("0.82.0", "round 1")) {
  parts <- if (nm == "0.82.0") base_parts else round1_parts
  tot <- rows <- nxi <- 0
  for (xi in seq(0.4, 0.9, by = 0.01)) {
    xs <- c(seq(-1 / xi - 0.01, -1 / xi - 60, by = -0.05),
            seq(-19.4, -20.8, by = -0.005))
    b <- 0
    for (x in xs) {
      b <- b + any(!is.finite(with_parts(parts,
                                         grad(0.3 + exp(-0.2) * x, xi))))
    }
    tot <- tot + b
    rows <- rows + length(xs)
    nxi <- nxi + (b > 0)
  }
  cat(sprintf("DENSE %s: %d rows, non-finite gradient %d, at %d of 51 xi\n",
              nm, rows, tot, nxi))
}
