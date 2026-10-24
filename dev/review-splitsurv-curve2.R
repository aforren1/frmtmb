# frm_curve() cases the lane did not test: a polynomial (no smooth) and
# a cox() fit's log hazard ratio; print() with dots. Lane lib.
.libPaths(c("C:/Users/adf44/source/r/wt-splitsurv-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.survival)})
cat("frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")), "\n")
set.seed(11)
d <- data.frame(x = runif(200))
d$y <- 1 + 2 * d$x - 3 * d$x^2 + rnorm(200, 0, 0.3)
fp <- frm(y ~ poly(x, 2, raw = TRUE), data = d)
g <- data.frame(x = seq(0.05, 0.95, length.out = 7))
cp <- frm_curve(fp, newdata = g, nsim = 2000, seed = 1)
lp <- frm_linpred(fp, newdata = g, se.fit = TRUE)
cat("poly curve vs frm_linpred: max |est diff| =",
    max(abs(cp$.estimate - lp$fit)), " max se ratio-1 =",
    max(abs(cp$.se / lp$se.fit - 1)), "\n")
dp <- frm_curve_deriv(fp, var = "x", newdata = g, simultaneous = FALSE)
cat("poly derivative vs 2 - 6x: max |diff| / se =",
    max(abs(dp$.estimate - (fixef(fp)[2] + 2 * fixef(fp)[3] * g$x)) / dp$.se),
    "\n")
fm <- frm_curve_feature(fp, var = "x", type = "maximum",
                        newdata = data.frame(x = seq(0.05, 0.95, length.out = 41)))
b <- fixef(fp)
cat("poly maximum:", fm$x, " closed form:", -b[2] / (2 * b[3]), "\n")
r <- tryCatch(print(cp, foo = 1), error = function(e) conditionMessage(e))
cat("print(curve, foo = 1):", if (is.character(r)) r else "printed", "\n")
r <- tryCatch({capture.output(print(cp, digits = 3)); "printed"},
              error = function(e) conditionMessage(e))
cat("print(curve, digits = 3):", r, "\n")
# cox(): the log hazard ratio as a curve in a smooth covariate
set.seed(2)
n <- 400
dd <- data.frame(x = runif(n))
dd$time <- rexp(n, exp(sin(2 * pi * dd$x)))
dd$cens <- rbinom(n, 1, 0.2)
fc <- frm(bf(time | cens(cens) ~ s(x, k = 6)), family = cox(), data = dd)
gc <- data.frame(x = seq(0.1, 0.9, length.out = 9))
cc <- tryCatch(frm_curve(fc, newdata = gc, nsim = 2000, seed = 2),
               error = function(e) conditionMessage(e))
if (is.character(cc)) cat("cox frm_curve ERROR:", cc, "\n") else {
  cat("cox log-HR curve (centered) vs truth sin(2 pi x):\n")
  print(round(cbind(x = gc$x, est = cc$.estimate - mean(cc$.estimate),
                    truth = sin(2 * pi * gc$x) - mean(sin(2 * pi * gc$x)),
                    se = cc$.se), 3))
}
