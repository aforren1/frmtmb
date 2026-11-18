# Review of lane slopetrig, check 5: the fits the ratio rule newly
# engages that reach the maximum and warn (36 in the lane's grids). For
# each: the default's warnings in full, VarCorr()'s warnings, and whether
# the correlation in x's own units is determined: its value and a
# delta-method SE from the fit on the centered, standardized column (Z),
# whose correlation parameter is interior.
# Usage: Rscript slopetrig-review-warn36.R [max cases]
a <- commandArgs(TRUE)
nmax <- if (length(a)) as.integer(a[1]) else 36L
.libPaths(c("C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
rd <- function(p) utils::read.delim(file.path("dev/slopetrig-log", p))
B <- rbind(rd("calib-gau-base.tsv"), rd("calib-pois-base.tsv"),
           rd("calib2-gau-base.tsv"), rd("calib2-pois-base.tsv"))
L <- rbind(rd("calib-gau-lane.tsv"), rd("calib-pois-lane.tsv"),
           rd("calib2-gau-lane.tsv"), rd("calib2-pois-lane.tsv"))
best <- pmax(B$ll_ref, B$ll_F, B$ll_T, B$ll_D, L$ll_D, na.rm = TRUE)
sel <- which(L$as_D & B$spread >= 0.05 & best - L$ll_D <= 1e-3 & L$w_D > 0)
cat("cases", length(sel), "\n")
mk <- function(cell, seed, fam) {
  set.seed(seed)
  n <- 300
  g <- factor(rep(1:20, each = 15))
  trt <- rep(c(TRUE, FALSE), length.out = n)
  z <- rnorm(n)
  b0 <- rnorm(20, 0, 0.5)[g]
  b1 <- rnorm(20, 0, 0.5)[g]
  c0 <- cell$c0
  r <- cell$r
  x <- switch(cell$design,
    unc = c0 + r * c0 * z, year = c0 + r * c0 * z,
    binlev = c0 + trt + 0.01 * z)
  eta <- if (cell$design == "binlev") 0.5 + b0 + (0.3 + b1) * (x - c0) else
    0.5 + b0 + b1 * z + 0.3 * z
  y <- if (fam == "gau") eta + rnorm(n) else rpois(n, exp(eta))
  data.frame(g = g, x = x, y = y)
}
grab <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
rho <- function(t) t / sqrt(1 + t^2)
raw_cor <- function(th, m, sx) {
  s0 <- exp(th[1]); s1 <- exp(th[2]); r <- rho(th[3])
  S <- matrix(c(s0^2, r * s0 * s1, r * s0 * s1, s1^2), 2)
  A <- matrix(c(1, 0, -m / sx, 1 / sx), 2)
  Sr <- A %*% S %*% t(A)
  Sr[1, 2] / sqrt(Sr[1, 1] * Sr[2, 2])
}
wtab <- character(0)
for (i in utils::head(sel, nmax)) {
  cell <- B[i, ]
  d <- mk(list(design = cell$design, c0 = cell$c0, r = cell$r), cell$seed,
          cell$fam)
  fo <- if (cell$form == "x") y ~ x + (1 + x | g) else y ~ 1 + (1 + x | g)
  fam <- if (cell$fam == "gau") gaussian() else poisson()
  m <- mean(d$x); sx <- sd(d$x)
  dz <- d; dz$x <- (d$x - m) / sx
  fd <- grab(frm(fo, data = d, family = fam))
  fz <- grab(frm(fo, data = dz, family = fam))
  vc <- grab(VarCorr(fd$v))
  vz <- grab(VarCorr(fz$v))
  # theta of Z and its covariance
  pz <- fz$v$opt$par
  ti <- which(names(pz) == "theta")
  Vz <- tryCatch(ns$sdr_of(fz$v)$cov.fixed[ti, ti], error = function(e) NULL)
  thz <- pz[ti]
  cr <- raw_cor(thz, m, sx)
  se <- NA
  if (!is.null(Vz) && all(is.finite(Vz))) {
    gr <- vapply(1:3, function(k) {
      h <- 1e-6 * max(1, abs(thz[k])); e <- replace(numeric(3), k, h)
      (raw_cor(thz + e, m, sx) - raw_cor(thz - e, m, sx)) / (2 * h)
    }, 0)
    se <- sqrt(drop(t(gr) %*% Vz %*% gr))
  }
  cd <- fd$v$opt$par[which(names(fd$v$opt$par) == "theta")]
  cat(sprintf(paste0("\n%s c0=%g r=%g %s form=%s seed=%d ratio=%.4g ll D-Z %.3g\n",
                     "  raw cor D %.10f ; mapped from Z %.10f ; delta SE from Z %.3g ;",
                     " 1+cor %.3g ; t D %.4g ; Z warnings %d, VarCorr(Z) warnings %d\n"),
              cell$design, cell$c0, cell$r, cell$fam, cell$form, cell$seed,
              cell$ratio, as.numeric(logLik(fd$v)) - as.numeric(logLik(fz$v)),
              rho(cd[3]), cr, se, 1 - abs(cr), cd[3], length(fz$w),
              length(vz$w)))
  for (w in fd$w) cat("  FIT WARN:", substr(gsub("\n", " ", w), 1, 400), "\n")
  for (w in vc$w) cat("  VARCORR WARN:", substr(gsub("\n", " ", w), 1, 400), "\n")
  vd <- vc$v$g
  cat("  VarCorr D cor:", format(vd$cor[1, , 2], digits = 8), "\n")
  wtab <- c(wtab, substr(fd$w, 1, 70))
}
cat("\nwarning heads over the cases:\n")
print(sort(table(wtab), decreasing = TRUE))
cat("DONE\n")
