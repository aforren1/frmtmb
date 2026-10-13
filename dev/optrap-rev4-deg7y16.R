# Reviewer of lane optrap, consolidation check: raw polynomials of degree
# 6 and 7 with (1 | g), y x 1 and 16 (dev/optrap-rev4-gate.R's data),
# on any build, against lmer: short, silent, code != 0, warned.
#   Rscript dev/optrap-rev4-deg7y16.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
lab <- a[2]
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
for (deg in 6:7) for (yk in c(1, 16)) {
  sh <- si <- cd <- wn <- 0L
  mn <- 0
  for (s in 1:20) {
    set.seed(1000 * deg + s)
    d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
    d$y <- (sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)) * yk
    for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
    fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                           "+ (1 | g)"))
    ref <- as.numeric(logLik(suppressMessages(suppressWarnings(
      lmer(fo, data = d, REML = FALSE)))))
    w <- 0L
    f <- withCallingHandlers(suppressMessages(frm(fo, data = d)),
      warning = function(x) {w <<- w + 1L; invokeRestart("muffleWarning")})
    gap <- as.numeric(logLik(f)) - ref
    if (gap < -1e-4) {
      sh <- sh + 1L
      if (f$opt$convergence == 0 && w == 0) si <- si + 1L
    }
    mn <- min(mn, gap)
    cd <- cd + (f$opt$convergence != 0)
    wn <- wn + (w > 0)
  }
  cat(sprintf("%s deg %d (1 | g) y x %g | short %d (silent %d, min %.3g) | code != 0 %d | warned %d\n",
              lab, deg, yk, sh, si, mn, cd, wn))
}
