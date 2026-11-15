# Re-check of lane famend after fix round 2: on inside rows the new zc
# (r - (ad_pos(ex) - ad_pos(ex - out))) must pass an adjoint through
# unchanged, also a subnormal one, where fl(0.5 A) is not A / 2. Tapes
# w1 * l + w2 * h of gev_parts() with odd subnormal weights, lane against
# 0.82.0's gev_parts() (where zc = z on an inside row): gradient and
# Hessian bit for bit.
#   Rscript dev/reviews/famend-short-subnorm.R <lib>
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
new_parts <- ns$gev_parts
tiny <- 2^-1074
ws <- list(c(tiny, 0), c(3 * tiny, 0), c(0, tiny), c(0, 7 * tiny),
           c(5 * tiny, 9 * tiny), c(1001 * tiny, 0), c(0, 2^-1022 + tiny),
           c(1, 0), c(0, 1), c(1e-300, 3e-310))
run <- function(parts, q, p, w) {
  tp <- RTMB::MakeTape(function(p) {
    g <- parts(q, list(mu = p[1], sigma = exp(p[2]), xi = p[3]))
    w[1] * g$l + w[2] * g$h
  }, p)
  list(v = tp(p), g = tp$jacobian(p), h = tp$jacfun()$jacobian(p))
}
n <- 0
same <- 0
sub_seen <- 0
diffs <- character()
cl <- -1 + 1e-15
for (xi in c(-0.9, -0.5, -1e-4, 0, 1e-8, 1e-4, 0.3, 0.6, 1, 10)) {
  # z from just inside the clamp to well inside, and z = 0
  zs <- c(cl + (1:6) * 2^-53, -0.999, -0.5, -1e-3, 0, 1e-3, 0.5, 3, 40)
  for (z in zs) {
    x <- if (xi == 0) z * 50 else z / xi
    q <- 0.3 + exp(-0.2) * x
    # the row as the tape forms it must be inside
    if (!(xi * ((q - 0.3) / exp(-0.2)) > cl)) next
    p <- c(0.3, -0.2, xi)
    for (w in ws) {
      o <- run(old_parts, q, p, w)
      nw <- run(new_parts, q, p, w)
      n <- n + 1
      sub_seen <- sub_seen + any(abs(c(o$g, o$h)) > 0 &
                                 abs(c(o$g, o$h)) < 2^-1022)
      if (identical(o, nw)) {
        same <- same + 1
      } else if (length(diffs) < 10) {
        diffs <- c(diffs, sprintf("xi %g z %.17g w %s: g %s vs %s", xi, xi * x,
                                  paste(format(w), collapse = ","),
                                  paste(format(o$g), collapse = " "),
                                  paste(format(nw$g), collapse = " ")))
      }
    }
  }
}
cat(sprintf(paste0("SUBNORMAL ADJOINT: %d inside tapes (rows x weights), ",
                   "bit-identical value, gradient and Hessian %d, differ %d; ",
                   "tapes with a subnormal gradient or Hessian entry %d\n"),
            n, same, n - same, sub_seen))
if (length(diffs)) cat(diffs, sep = "\n")
