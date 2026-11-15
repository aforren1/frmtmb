# Lane famend, fix round F1: the GEV's density, CDF and log survivor on
# the tape, with this build's gev_parts() against 0.82.0's (read from
# dev/famend-out/families-extreme-head.R, which the caller writes with
# git show HEAD:R/families-extreme.R), in one process: values and
# gradients bit for bit on rows inside the support, and the count of
# non-finite gradients and Hessians on outside rows. Fix round 2 adds
# xi 0.44 to 0.84, where 0.82.0's clamped z rounded above the clamp, a
# dense xi grid 0.4 to 0.9, and the outside rows' value change.
#   Rscript dev/famend-gevbit.R <lib>
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
new_parts <- ns$gev_parts
fam <- gen_extreme_value()
# the family's three slots with gev_parts() swapped in place
with_parts <- function(parts, expr) {
  orig <- get("gev_parts", ns)
  unlockBinding("gev_parts", ns)
  assign("gev_parts", parts, ns)
  on.exit({ assign("gev_parts", orig, ns); lockBinding("gev_parts", ns) })
  expr
}
eval_slot <- function(slot, y, xi) {
  p <- c(0.3, -0.2, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam[[slot]](y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                        .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  list(v = tp(p), g = tp$jacobian(p), h = tp$jacfun()$jacobian(p))
}
xis <- c(-0.9, -0.5, -0.1, 0, 0.001, 0.01, 0.06, 0.3, 0.44, 0.5, 0.6, 0.7,
         0.84, 1, 3, 10)
same <- 0
diff_in <- 0
n_in <- 0
for (xi in xis) {
  # rows inside the support, t from 1e-14 up, and far out on the open side
  lo <- if (xi > 0) -1 / xi else -60
  hi <- if (xi < 0) -1 / xi else 60
  ys <- 0.3 + exp(-0.2) * c(seq(lo, hi, length.out = 401)[-c(1, 401)],
                             if (xi > 0) lo + c(1e-13, 1e-10, 1e-6) / xi,
                             if (xi < 0) hi - c(1e-13, 1e-10, 1e-6) / abs(xi))
  for (slot in c("lpdf", "lcdf", "lccdf")) {
    for (y in ys) {
      o <- with_parts(old_parts, eval_slot(slot, y, xi))
      nw <- eval_slot(slot, y, xi)
      n_in <- n_in + 1
      if (identical(o, nw)) same <- same + 1 else diff_in <- diff_in + 1
    }
  }
}
cat(sprintf("INSIDE rows x slots: %d, bit-identical value, gradient and Hessian %d, differ %d\n",
            n_in, same, diff_in))
# outside rows: the 6000-point grid of item 3067 and a wider one; on
# every 5th row also 0.82.0's gradient and the relative change of the
# value (fix round 2)
for (xi in c(0.001, 0.01, 0.06, 0.3, 0.44, 0.5, 0.6, 0.7, 0.84, 1, 3, 10)) {
  xs <- seq(-1 / xi - 0.01, -1 / xi - 60, by = -0.01)
  if (xi < 0.05) xs <- seq(-1 / xi - 0.01, -1 / xi - 60, length.out = 6000)
  bad_n <- bad_o <- badh_n <- nanh_n <- 0
  rel <- 0
  for (j in seq_along(xs)) {
    y <- 0.3 + exp(-0.2) * xs[j]
    nw <- eval_slot("lpdf", y, xi)
    bad_n <- bad_n + any(!is.finite(nw$g))
    badh_n <- badh_n + any(!is.finite(nw$h))
    nanh_n <- nanh_n + any(is.nan(nw$h))
    if (j %% 5 == 0) {
      o <- with_parts(old_parts, eval_slot("lpdf", y, xi))
      bad_o <- bad_o + any(!is.finite(o$g))
      rel <- max(rel, abs(nw$v - o$v) / abs(o$v))
    }
  }
  cat(sprintf(paste0("OUTSIDE xi %5g: %d rows; non-finite gradient %d, ",
                     "Hessian %d (NaN %d) (this build); 0.82.0's ",
                     "gev_parts, every 5th row: %d of %d, lpdf relative ",
                     "change max %.3g\n"),
              xi, length(xs), bad_n, badh_n, nanh_n, bad_o,
              length(xs) %/% 5, rel))
}
# a dense xi grid, 0.4 to 0.9 by 0.01, rows 0.05 apart from the lower
# end to 60 below it and 0.005 apart over the band x in [-20.8, -19.4]
tot_g <- tot_nan <- tot_rows <- 0
for (xi in seq(0.4, 0.9, by = 0.01)) {
  xs <- c(seq(-1 / xi - 0.01, -1 / xi - 60, by = -0.05),
          seq(-19.4, -20.8, by = -0.005))
  bg <- bh <- 0
  for (x in xs) {
    nw <- eval_slot("lpdf", 0.3 + exp(-0.2) * x, xi)
    bg <- bg + any(!is.finite(nw$g))
    bh <- bh + any(is.nan(nw$h))
  }
  tot_g <- tot_g + bg
  tot_nan <- tot_nan + bh
  tot_rows <- tot_rows + length(xs)
  if (bg || bh) {
    cat(sprintf("DENSE xi %.2f: %d rows, non-finite gradient %d, NaN Hessian %d\n",
                xi, length(xs), bg, bh))
  }
}
cat(sprintf("DENSE xi 0.40 to 0.90 by 0.01: %d rows, non-finite gradient %d, NaN Hessian %d\n",
            tot_rows, tot_g, tot_nan))
