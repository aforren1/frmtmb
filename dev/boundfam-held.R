# Lane boundfam: fitted() and predict() standard errors on 3269 (a
# zero_inflated_binomial whose zi intercept ran to its limit) and on the
# all-zero poisson, on one build. Prints values to 12 digits for a diff.
#   Rscript dev/boundfam-held.R <lib or "base">
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
set.seed(7)
d <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
d$y <- stats::rbinom(200, d$N, stats::plogis(0.3 + 0.5 * d$x))
fz <- cap(frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
              data = d))$v
fb <- frm(bf(y | trials(N) ~ x), family = binomial(), data = d)
r <- cap(fitted(fz))
cat("fitted warnings:", length(r$w), "NaN rows",
    sum(is.nan(r$v[, "Est.Error"])), "\n")
cat("fitted Est.Error[1:6]:", sprintf("%.12g", r$v[1:6, "Est.Error"]), "\n")
eb <- fitted(fb)[, "Est.Error"]
cat("max |zib - binomial| fitted Est.Error:",
    format(max(abs(r$v[, "Est.Error"] - eb)), digits = 4), "\n")
set.seed(11)
p <- cap(predict(fz, ndraws = 400))
cat("predict warnings:", length(p$w), "NaN rows",
    sum(is.nan(p$v[, "Est.Error"])), "\n")
cat("predict Est.Error[1:6]:", sprintf("%.12g", p$v[1:6, "Est.Error"]), "\n")
writeLines(sprintf("%.15g", c(r$v[, "Est.Error"], p$v[, "Est.Error"])),
           file.path(a[2], paste0("held-", if (a[1] == "base") "base" else "lane", ".txt")))
# the all-zero poisson: x is undetermined, a fitted value reads it
set.seed(2)
de <- data.frame(x = stats::rnorm(50), y = 0L)
f0 <- cap(frm(y ~ x, family = poisson(), data = de))$v
r0 <- cap(fitted(f0, newdata = data.frame(x = c(-1, 0, 1))))
cat("pois0 fitted Est.Error:", r0$v[, "Est.Error"], " warnings:",
    length(r0$w), substr(r0$w[1], 1, 80), "\n")
