# Reviewer, punch round 2, item 2: the value across both ties at the
# largest shapes, against Rmpfr at the worst points of a 401-point scan.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(RTMB); library(Rmpfr) })
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-rev3-ibeta.R")
eval(parse(text = src[grep("^prec <- 256", src):(grep("^# derivatives of g", src) - 1)]))
lib <- frmtmb:::log_ibeta_half
for (ab in list(c(3e4, 7e4), c(1e6, 3e6), c(1e7, 1e7), c(3e6, 1e7))) {
  a <- ab[1]; b <- ab[2]; s <- a + b; t1 <- a / s; gap <- b / (s * (s + 1))
  u <- seq(-1, 2, length.out = 401)
  xs <- t1 + u * gap; ok <- xs < 0.5; xs <- xs[ok]; u <- u[ok]
  e <- abs(lib(xs, a, b) - pbeta(xs, a, b, log.p = TRUE))
  w <- order(-e)[1:3]
  for (i in w) {
    r <- as.numeric(mp_I(xs[i], a, b))
    cat(sprintf("a %g b %g x = t1 %+.3f gap: ours vs mpfr %.1e, pbeta vs mpfr %.1e\n", a, b, u[i],
                abs(lib(xs[i], a, b) - r) / max(1, abs(r)),
                abs(pbeta(xs[i], a, b, log.p = TRUE) - r) / max(1, abs(r))))
  }
}
