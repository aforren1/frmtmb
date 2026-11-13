# Lane boundfam, consolidation item 1: fitted() Est.Error of the two fits
# of ppdraw's test-predict-draw-edge.R (zi ~ x at zi's end; hu ~ x where
# no y is 0), and of a separating bernoulli and the all-zero poisson,
# written to a file for a diff between builds.
#   Rscript dev/boundfam-cons-fitted.R <lib or "base"> <out file>
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
q <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
set.seed(1)
d <- data.frame(x = stats::rnorm(200), N = 10L)
p <- stats::plogis(-0.5 + 0.8 * d$x)
d$y <- ifelse(stats::runif(200) < 0, 0L, stats::rbinom(200, d$N, p))
f1 <- q(frm(bf(y | trials(N) ~ x, zi ~ x), family = zero_inflated_binomial(),
            data = d))$v
set.seed(12)
dh <- data.frame(x = stats::rnorm(300))
dh$y <- stats::rpois(300, exp(2.5 + 0.2 * dh$x))
f4 <- q(frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(), data = dh))$v
set.seed(514)
ds <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
ds$yb <- as.integer(ds$z > 0)
fs <- q(frm(yb ~ z + x, family = bernoulli(), data = ds))$v
set.seed(2)
de <- data.frame(x = stats::rnorm(50), y = 0L)
f0 <- q(frm(y ~ x, family = poisson(), data = de))$v
out <- character()
for (nm in c("zi", "hu", "sepbern", "pois0")) {
  f <- switch(nm, zi = f1, hu = f4, sepbern = fs, pois0 = f0)
  r <- q(fitted(f))
  e <- r$v[, "Est.Error"]
  cat(sprintf("%-8s NaN %d of %d, warnings %d, lost %s\n", nm,
              sum(is.nan(e)), length(e), length(r$w),
              paste(names(frmtmb:::sdr_of(f)$se_lost),
                    frmtmb:::sdr_of(f)$se_lost, collapse = " ")))
  out <- c(out, paste(nm, sprintf("%.15g", e)))
}
writeLines(out, a[2])
