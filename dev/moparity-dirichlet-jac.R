# The Dirichlet on a mo() simplex in each chart (dirichlet_logdens()):
# its log-Jacobian against finite differences, and its normalization by
# integration over the chart for D = 3.
# Usage: Rscript dev/moparity-dirichlet-jac.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args) && nzchar(args[[1]])) args[[1]] else
  "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ms <- frmtmb:::mo_simplex
dd <- frmtmb:::dirichlet_dist
dl <- frmtmb:::dirichlet_logdens
# |det d(w_2..w_D)/dz| by central differences
fd_logjac <- function(z, chart) {
  f <- function(z) ms(z, chart)[-1L]
  h <- 1e-6
  J <- vapply(seq_along(z), function(j) {
    e <- numeric(length(z))
    e[j] <- h
    (f(z + e) - f(z - e)) / (2 * h)
  }, numeric(length(z)))
  log(abs(det(J)))
}
set.seed(1)
worst <- 0
for (D in 3:5) {
  a <- stats::runif(D, 0.5, 4)
  for (chart in list(NULL, "softmax")) {
    for (i in 1:20) {
      z <- stats::rnorm(D - 1)
      w <- ms(z, chart)
      ld_w <- lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * log(w))
      sheets <- if (is.null(chart)) D * log(2) else 0
      want <- ld_w + fd_logjac(z, chart) - sheets
      got <- dl(z, dd(a, D, chart))
      worst <- max(worst, abs(got - want) / max(1, abs(want)))
    }
  }
}
cat(sprintf("JAC worst relative gap (D 3..5, both charts, 20 points each): %.3g\n",
            worst))
# normalization over the chart, D = 3: the density integrates to one
a <- c(2, 1.3, 0.7)
for (chart in list(NULL, "softmax")) {
  dist <- dd(a, 3L, chart)
  # polar grid in z: radius r = tan(t) maps (0, pi/2) onto (0, Inf)
  nt <- 1200
  nphi <- 1200
  tt <- (seq_len(nt) - 0.5) / nt * (pi / 2)
  ph <- (seq_len(nphi) - 0.5) / nphi * 2 * pi
  tot <- 0
  for (t in tt) {
    r <- tan(t)
    jac <- r / cos(t)^2
    v <- vapply(ph, function(p) exp(dl(c(r * cos(p), r * sin(p)), dist)),
                0)
    tot <- tot + sum(v) * jac
  }
  tot <- tot * (pi / 2 / nt) * (2 * pi / nphi)
  cat(sprintf("NORM chart %-14s integral %.6f\n",
              if (is.null(chart)) "stereographic" else chart, tot))
}
