# Lane famend, fix round 2: the inside rows near xi = 0 of
# dev/reviews/famend-short-gevclamp.R (part 1) whose value, gradient or
# Hessian differ from 0.82.0's gev_parts(): print each with the entries
# that differ.
#   Rscript dev/famend-fix2-xi0diff.R <lib>
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
fmt <- function(v) paste(format(v, digits = 6), collapse = " ")
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
      if (!identical(o, nw)) {
        xx <- (y - mu0) / exp(ls0)
        cat(sprintf("xi %g %s x %.10g (tape-side z %.17g)\n", xi, slot, x,
                    expm1(log1p(xi)) * xx))
        cat("  v old", fmt(o$v), "| new", fmt(nw$v), "\n")
        cat("  g old", fmt(o$g), "\n  g new", fmt(nw$g), "\n")
        cat("  h old", fmt(o$h), "\n  h new", fmt(nw$h), "\n")
      }
    }
  }
}
