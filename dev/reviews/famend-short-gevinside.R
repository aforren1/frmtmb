# Short check of lane famend, F1 attack: rows k ulps above the clamp,
# built so that the tape's own z (xi through log1p and expm1, x through
# mu and sigma) is above -1 + 1e-15. Lane gev_parts() against 0.82.0's,
# value, gradient and Hessian bit for bit; and the value change on
# outside rows away from the band.
#   Rscript dev/reviews/famend-short-gevinside.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
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
  p <- c(0, 0, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam[[slot]](y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                        .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  list(v = tp(p), g = tp$jacobian(p), h = tp$jacfun()$jacobian(p))
}
cl <- -1 + 1e-15
n <- 0
same <- 0
skipped <- 0
for (xi in c(0.06, 0.3, 0.5, 1, 3, 10, -0.5, -0.9)) {
  xit <- expm1(log1p(xi))
  for (k in 1:12) {
    zt <- cl + k * .Machine$double.eps / 2
    x <- zt / xit
    # nudge x until the tape's product lands strictly above the clamp
    for (j in 1:4) if (!(xit * x > cl)) x <- x * (1 - sign(xit) * 1e-16)
    if (!(xit * x > cl)) {
      skipped <- skipped + 1
      next
    }
    for (slot in c("lpdf", "lcdf", "lccdf")) {
      n <- n + 1
      same <- same + identical(ev(slot, x, xi),
                               with_parts(old_parts, ev(slot, x, xi)))
    }
  }
}
cat(sprintf("JUST INSIDE (tape z above the clamp): %d rows x slots, bit-identical %d, differ %d; rows skipped %d\n",
            n, same, n - same, skipped))
# outside rows away from the band: relative change of the value
for (xi in c(0.3, 1, 3)) {
  rel <- c()
  for (x in c(-1 / xi - 0.5, -1 / xi - 2, -5 / xi, -15, -30)) {
    vn <- ev("lpdf", x, xi)$v
    vo <- with_parts(old_parts, ev("lpdf", x, xi))$v
    rel <- c(rel, abs(vn - vo) / abs(vo))
  }
  cat(sprintf("OUTSIDE xi %g, lpdf relative change at 5 rows: %s\n", xi,
              paste(format(rel, digits = 3), collapse = " ")))
}
