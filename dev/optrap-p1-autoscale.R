# Lane optrap, punch round 1: test-autoscale.R's "autoscale covers dpar
# formulas" fits, the two coefficients the test compares, and the units.
#   Rscript dev/optrap-p1-autoscale.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(606)
n <- 300
g <- factor(rep(1:30, each = 10))
xbig <- stats::rnorm(n) * 1e5
b0 <- stats::rnorm(30, 0, 0.4)
dd <- data.frame(xbig = xbig, g = g)
dd$y <- stats::rnorm(n, 1 + 5e-6 * xbig + b0[g], exp(-0.5 + 2e-6 * xbig))
form <- bf(y ~ xbig + (1 | g), sigma ~ xbig) + gaussian()
fa <- frm(form, data = dd, control = frmtmb_control(autoscale = TRUE))
sb <- stats::sd(dd$xbig)
dd2 <- dd
dd2$xz <- (dd$xbig - mean(dd$xbig)) / sb
fr <- frm(bf(y ~ xz + (1 | g), sigma ~ xz) + gaussian(), data = dd2)
for (f in list(fa, fr)) {
  cat("units", deparse(frmtmb:::predictor_units(f$frame)), "par_units",
      deparse(signif(f$par_units, 3)), "obj", format(f$opt$objective,
                                                    digits = 12),
      "code", f$opt$convergence, "\n")
}
print(c(fa_mu = fixef_by_dpar(fa)$mu[["xbig"]],
        fr_mu = fixef_by_dpar(fr)$mu[["xz"]] / sb,
        fa_sig = fixef_by_dpar(fa)$sigma[["xbig"]],
        fr_sig = fixef_by_dpar(fr)$sigma[["xz"]] / sb), digits = 10)
