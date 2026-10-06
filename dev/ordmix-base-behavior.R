# What the base build (rellib-r5, 0.67.0) and the lane build do with
# the models this lane adds, and one silent defect of the base it
# found. Data seed 20261005. Usage: Rscript dev/ordmix-base-behavior.R
# [lane|base]; outputs dev/ordmix-log-base-behavior-<arm>.txt.
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "base"
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", find.package("frmtmb"), "\n")
set.seed(20261005)
n <- 400
x <- rnorm(n)
z <- rnorm(n)
g <- factor(sample(c("a", "b"), n, TRUE))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
yh <- ifelse(runif(n) < 0.2, 0L, y)
d <- data.frame(y, yh, x, z, g)
tryf <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-40s ERROR: %s\n", label, conditionMessage(r)))
  } else {
    cat(sprintf("%-40s logLik %.10f\n", label, as.numeric(logLik(r))))
  }
  invisible(r)
}
tryf("mixture(cumulative, cumulative)",
     frm(bf(y ~ x), family = mixture(cumulative(), cumulative()), data = d))
tryf("mixture(..., order = 'mu')",
     frm(bf(y ~ x), family = mixture(cumulative(), cumulative(),
                                     order = "mu"), data = d))
tryf("hurdle_cumulative + thres(gr = g)",
     frm(bf(yh | thres(gr = g) ~ x), family = hurdle_cumulative(), data = d))
tryf("hurdle_cumulative + cs(x)",
     frm(bf(yh ~ cs(x)), family = hurdle_cumulative("probit"), data = d))
# The base build put the cs() offsets of ANY predictor into the one
# `.cs` slot the ordinal densities read, so cs() in the formula of disc
# moved the thresholds as if it were written in mu's: the fit below is
# the fit of y ~ x + cs(z), to the last bit, under another name.
a <- tryf("sratio: disc ~ cs(z)",
          frm(bf(y ~ x, disc ~ cs(z)), family = sratio(), data = d))
b <- tryf("sratio: y ~ x + cs(z)",
          frm(bf(y ~ x + cs(z)), family = sratio(), data = d))
if (!inherits(a, "error") && !inherits(b, "error")) {
  cat("identical logLik:", identical(as.numeric(logLik(a)),
                                     as.numeric(logLik(b))), "\n")
  print(fixef(a))
}
