# Item 3.5: is the contaminant density's gradient accurate for a row
# just above the non-decision time? Compares the tape gradient of the
# anchored form against a central difference of the same function, and
# against the data-anchored form, at rows 1e-2 to 1e-7 s above ndt.
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib2",
            "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all("extensions/frmtmb.eam", quiet = TRUE,
                                   export_all = TRUE))
cr <- c(0.1, 5)
lg <- -log(2 * diff(cr))
old <- function(lf, l, l1m) {
  d <- lg - lf
  m <- 0.5 * (d + abs(d))
  lf + m + log(ddm_floor(exp(l1m - m) + exp(l + d - m), 1e-300))
}
new <- function(lf, l, l1m) {
  lg + log(ddm_floor(exp(l1m + lf - lg) + exp(l), 1e-300))
}
for (gap in 10^-(2:7)) {
  y <- 0.2 + gap
  for (form in c("old", "new")) {
    fn <- get(form)
    f <- function(p) {
      lf <- ddm_lpdf_both(y - p[1], 0.8, exp(p[2]), 0.5, 1)
      eta <- -3
      fn(lf, -log1p(exp(-eta)), -log1p(exp(eta)))
    }
    tp <- RTMB::MakeTape(f, c(0.2, log(1.4)))
    p <- c(0.2, log(1.4))
    g <- tp$jacobian(p)
    h <- c(1e-3 * gap, 1e-6)
    fd <- c((tp(p + c(h[1], 0)) - tp(p - c(h[1], 0))) / (2 * h[1]),
            (tp(p + c(0, h[2])) - tp(p - c(0, h[2]))) / (2 * h[2]))
    cat(sprintf("gap %.0e %s: value %.6g  grad (ndt, log a) = %s  fd = %s\n",
                gap, form, tp(p), paste(format(g, digits = 4), collapse = " "),
                paste(format(fd, digits = 4), collapse = " ")))
  }
}

# A row BELOW the non-decision time, under the floored density a stated
# max_ndt selects: the value of each form as the boundary separation
# moves by one part in 1e9. The contaminant alone explains the row, so
# the true value does not depend on the separation at all.
delta <- 1e-10
cat("\nrow 0.05 s below ndt, value minus its own first entry:\n")
for (form in c("old", "new")) {
  fn <- get(form)
  v <- vapply(log(1.4) + (0:5) * 1e-9, function(la) {
    lf <- ddm_lpdf_both(ddm_floor(0.15 - 0.2 - delta, delta), 0.8, exp(la),
                        0.5, 1)
    fn(lf, -log1p(exp(3)), -log1p(exp(-3)))
  }, 0)
  cat(form, ":", format(v - v[1], digits = 3), "\n")
}
