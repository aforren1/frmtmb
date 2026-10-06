# A cumulative fit whose response skips an interior category, on the
# base build: the increment between the two thresholds that bound the
# empty category runs to 0, the thresholds collapse in floating point,
# and log(F(a) - F(b)) meets log(0). Seed 20261005.
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "base"
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(20261005)
n <- 300
x <- rnorm(n)
lat <- 0.8 * x + rlogis(n)
y <- 1L + (lat > -1) + (lat > 0.5) + (lat > 1.5)
y[y == 3L] <- 4L
d <- data.frame(y = y, x = x)
print(table(d$y))
for (lk in c("logit", "probit")) {
  r <- tryCatch(frm(bf(y | thres(3) ~ x), family = cumulative(lk), data = d),
                error = function(e) e, warning = function(w) w)
  cat(lk, ": ", if (inherits(r, "condition")) conditionMessage(r) else
    paste("logLik", format(as.numeric(logLik(r)), digits = 12),
          "tau_raw", paste(format(r$estimates$tau_raw, digits = 4),
                           collapse = " ")), "\n")
}
