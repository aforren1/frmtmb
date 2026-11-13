# Lane boundfam: test-latres.R's lognormal ar(cov = TRUE) fit, whose
# latent takes the residual, on one build: verdicts and the latres reason.
#   Rscript dev/boundfam-latres.R <lib or "base">
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(3L)
d <- expand.grid(t = 1:8, g = factor(1:25))
d <- d[order(d$g, d$t), ]
d$x <- stats::rnorm(200L)
e <- as.vector(replicate(25, stats::arima.sim(list(ar = 0.6), 8, sd = 0.5)))
d$y <- stats::rlnorm(200L, e, 2)
w <- character(); m <- character()
f <- withCallingHandlers(frm(y ~ x + ar(t, g, cov = TRUE), family = lognormal(), data = d),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") },
  message = function(x) { m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage") })
for (x in w) cat("W:", substr(x, 1, 400), "\n")
for (x in m) cat("M:", substr(x, 1, 400), "\n")
print(f$opt$par)
print(suppressWarnings(sqrt(diag(vcov(f, full = TRUE)))))
print(ns$sdr_of(f)$se_lost)
bk <- ns$frame_latres(f$frame)[[1L]]
cat("shift", ns$latres_shift_tk(f, bk), "\n")
