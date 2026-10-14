# Lane illse: frm_curve_deriv() on the degree-8 raw polynomial of
# frmtmb.spline's test-illcond.R, against the derivative's exact SE
# from lm()'s QR (ML sigma).
.libPaths(c("C:/Users/adf44/source/r/wt-illse-lib", "C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spline)})
set.seed(1)
d <- data.frame(x = stats::runif(200, 1, 2))
d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
fo <- y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8)
fit <- frm(fo, data = d)
g <- data.frame(x = seq(1.1, 1.9, length.out = 9))
d1 <- frm_curve_deriv(fit, var = "x", order = 1, newdata = g, simultaneous = FALSE)
lmf <- lm(fo, data = d); n <- 200; k <- 9
D <- outer(g$x, 0:8, function(x, j) ifelse(j == 0, 0, j * x^pmax(j - 1, 0)))
Ri <- backsolve(qr.R(lmf$qr), diag(k)); s <- sqrt(sum(residuals(lmf)^2) / n)
ref <- sqrt(rowSums((D %*% Ri)^2)) * s
print(d1$.se / ref)
