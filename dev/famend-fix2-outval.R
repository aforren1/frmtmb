# Lane famend, fix round 2: how the tape log density of an outside GEV
# row moves against 0.82.0's gev_parts(). Rows as in dev/famend-gevbit.R
# (mu 0.3, sigma exp(-0.2), x from the lower end to 60 below it by
# 0.01, every 5th row). Rows split by whether 0.82.0's clamped z (replayed
# on doubles) was the clamp itself: the relative change on each group,
# and the smallest |log density| among the rows that moved more than
# 1e-10 relative.
#   Rscript dev/famend-fix2-outval.R <lib>
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
val <- function(y, xi) {
  p <- c(0.3, -0.2, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam$lpdf(y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                     .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  tp(p)
}
cl <- -1 + 1e-15
ad_pos <- function(u) 0.5 * (u + abs(u))
old_zc <- function(z) {
  zc <- z + ad_pos(cl - z)
  zc + ad_pos(cl - zc)
}
for (xi in c(0.06, 0.3, 0.44, 0.5, 0.6, 0.7, 0.84, 1, 3, 10)) {
  xs <- seq(-1 / xi - 0.01, -1 / xi - 60, by = -0.01)
  xs <- xs[seq(5, length(xs), by = 5)]
  xit <- expm1(log1p(xi))
  at <- above <- 0
  rel_at <- rel_above <- 0
  min_moved <- Inf
  for (x in xs) {
    y <- 0.3 + exp(-0.2) * x
    vn <- val(y, xi)
    vo <- with_parts(old_parts, val(y, xi))
    r <- abs(vn - vo) / abs(vo)
    if (old_zc(xit * ((y - 0.3) / exp(-0.2))) == cl) {
      at <- at + 1
      rel_at <- max(rel_at, r)
    } else {
      above <- above + 1
      rel_above <- max(rel_above, r)
    }
    if (r > 1e-10) min_moved <- min(min_moved, abs(vo), abs(vn))
  }
  cat(sprintf(paste0("xi %5g: %d rows; 0.82.0 on the clamp %d, max relative ",
                     "change %.3g; above it %d, max relative change %.3g; ",
                     "smallest |lpdf| of a row moved > 1e-10: %.3g\n"),
              xi, length(xs), at, rel_at, above, rel_above, min_moved))
}
