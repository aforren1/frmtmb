# Lane optrap, punch round 1 (m10): raw polynomials y ~ x + ... + x^deg
# (the review's probe design, x on [1, 2], residual sd near 0.36),
# degrees 3 to 7, seeds 1 to 20, against lm(). Per degree: fits more
# than 1e-4 below lm(), the largest gap, nonzero codes, warned fits.
#   Rscript dev/optrap-p1-polystudy.R <lib> <label> [round0]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
if (identical(a[3], "round0")) source("dev/optrap-round0-unit.R")
# OPTRAP_POLY_ARM=nospread: every column spread 1 (the response's unit
# alone); =unitonly1: predictor_units() as is but the coefficients'
# units left out of par_units (the sd start and inner tolerance only)
arm <- Sys.getenv("OPTRAP_POLY_ARM")
if (arm == "nospread") {
  utils::assignInNamespace("predictor_col_spread", function(lp) {
    rep(1, NCOL(lp[["X"]]))
  }, "frmtmb")
}
if (arm == "nocoef") {
  utils::assignInNamespace("predictor_par_units", function(...) NULL,
                           "frmtmb")
}
cat("arm", arm, "
")
sdr <- as.numeric(Sys.getenv("OPTRAP_POLY_SD", "0.2"))
# OPTRAP_POLY_YK: the response multiplied by this (its log-likelihood
# and lm()'s shift alike, so the gaps compare across multipliers)
yk <- as.numeric(Sys.getenv("OPTRAP_POLY_YK", "1"))
cat("lib", dirname(find.package("frmtmb")), "label", a[2], "noise sd", sdr, "y x", yk,
    "\n")
X <- NULL
for (deg in 3:7) for (s in 1:20) {
  set.seed(1000 * deg + s)
  d <- data.frame(x = runif(200, 1, 2))
  d$y <- sin(3 * d$x) + rnorm(200, 0, sdr)
  d$y <- d$y * yk
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + ")))
  w <- 0L
  f <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
    w <<- w + 1L
    invokeRestart("muffleWarning")
  })
  gap <- as.numeric(logLik(f)) - as.numeric(logLik(lm(fo, data = d)))
  X <- rbind(X, data.frame(deg = deg, seed = s, gap = gap,
                           code = f$opt$convergence, warn = w,
                           unit = (f$par_units %||% 1)[1]))
}
for (dg in 3:7) {
  x <- X[X$deg == dg, ]
  sh <- x$gap < -1e-4
  cat(sprintf("deg %d | short %2d of 20 (min %9.3g) | silent short %2d | code!=0 %2d | warned %2d | units %s\n",
              dg, sum(sh), min(x$gap), sum(sh & x$code == 0 & x$warn == 0),
              sum(x$code != 0), sum(x$warn > 0),
              paste(sort(unique(signif(x$unit, 3))), collapse = ",")))
}
