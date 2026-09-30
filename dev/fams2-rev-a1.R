# Reviewer: xbeta on interior-only data from an xbeta truth (kappa = 1,
# phi = 200, mu = 0.5, seed 4101): where does the gradient go NaN?
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(4101)
n <- 2000
z <- rbeta(n, 0.5 * 200, 0.5 * 200)
y <- 3 * z - 1
d <- data.frame(y = y)
lp <- function(p) {
  mu <- plogis(p[1]); phi <- exp(p[2]); k <- exp(p[3]); dd <- 1 + 2 * k
  sum(dbeta((y + k) / dd, mu * phi, (1 - mu) * phi, log = TRUE) - log(dd))
}
o <- optim(c(0, log(5), log(0.1)), function(p) -lp(p), method = "BFGS",
           control = list(maxit = 5000, reltol = 1e-14))
cat("reference optimum: mu", plogis(o$par[1]), "phi", exp(o$par[2]),
    "kappa", exp(o$par[3]), "logLik", -o$value, "\n")
fb <- frm(y ~ 1, family = Beta(), data = d)
cat("Beta logLik", as.numeric(logLik(fb)), "\n")
cat("profile of the reference in log kappa (phi, mu re-optimized):\n")
for (lk in c(-6, -3, -1, 0, 0.5, 1)) {
  oo <- optim(c(0, 3), function(q) -lp(c(q, lk)), method = "BFGS")
  cat(sprintf("  log kappa %5.1f  logLik %.4f  phi %.1f\n", lk, -oo$value, exp(oo$par[2])))
}
f <- tryCatch(frm(y ~ 1, family = xbeta(), data = d, verbose = TRUE), error = function(e) e)
print(f)
f2 <- tryCatch(frm(y ~ 1, family = xbeta(), data = d,
                   control = frmtmb_control(optimizer = "optim")), error = function(e) e)
print(if (inherits(f2, "error")) f2 else confint(f2))
f3 <- tryCatch(frm(bf(y ~ 1, kappa = 1), family = xbeta(), data = d), error = function(e) e)
print(if (inherits(f3, "error")) f3 else logLik(f3))
