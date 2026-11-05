# Which warnings does frm() give as the noise shrinks? Lane build, so the
# nugget warning shows beside the SE check. Seed 20261102.
#   Rscript dev/gpgap-lowwarn.R
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib","C:/Users/adf44/source/r/rellib-r16","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
for (s in c(1e-2,5e-3,3e-3,1e-3,1e-4)) {
set.seed(20261102)
n <- 60
x <- round(stats::runif(n, 0, 6), 3)
d <- data.frame(x = x, y = sin(2 * x) + s * stats::rnorm(n))
w <- character()
fit <- withCallingHandlers(frm(bf(y ~ gp(x)), data = d), warning = function(c) {w <<- c(w, class(c)[1], substr(conditionMessage(c),1,400)); invokeRestart("muffleWarning")})
cat("true", s, "sigma", sigma(fit), "\n"); print(w)
}
