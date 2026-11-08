# Probe: at an unseen gp() position, which Jacobian is right: marginaleffects'
# forward difference over get_coef(), or frm_linpred(se.fit)'s analytic one?
# A central difference with a fresh cache per step is the referee. Seed 8150.
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
ns <- asNamespace("frmtmb")
set.seed(8150)
m <- 90
dg <- data.frame(x = runif(m, 0, 6), z = runif(m, 0, 6))
dg$y <- sin(dg$x) + 0.3 * dg$z + rnorm(m, sd = 0.3)
fit <- quiet(frm(y ~ z + gp(x), data = dg))
ndg <- data.frame(x = c(dg$x[5], 3.33, 7, 12), z = 3)
cf <- get_coef(fit)
V <- get_vcov(fit)[names(cf), names(cf)]
f <- function(v) {
  g <- set_coef(fit, v)
  g$cache <- new.env(parent = emptyenv())
  as.numeric(frm_linpred(g, newdata = ndg))
}
for (h in c(1e-3, 1e-4, 1e-5)) {
  J <- vapply(seq_along(cf), function(i) {
    e <- numeric(length(cf)); e[i] <- h * max(1, abs(cf[i]))
    (f(cf + e) - f(cf - e)) / (2 * e[i])
  }, numeric(nrow(ndg)))
  cat(sprintf("central h %g: se %s\n", h,
              paste(sprintf("%.6f", sqrt(rowSums((J %*% V) * J))),
                    collapse = " ")))
}
l <- frm_linpred(fit, newdata = ndg, se.fit = TRUE)
lb <- frm_lp_basis(fit, newdata = ndg)
cat("frm_linpred se:", sprintf("%.6f", l$se.fit), "\n")
cat("basis extra_var:", lb$extra_var, "\n")
cat("basis A V A' + extra:", sprintf("%.6f",
    sqrt(rowSums((as.matrix(lb$A) %*% lb$V) * as.matrix(lb$A)) +
           lb$extra_var)), "\n")
cat("basis A V A':", sprintf("%.6f",
    sqrt(rowSums((as.matrix(lb$A) %*% lb$V) * as.matrix(lb$A)))), "\n")
p <- quiet(predictions(fit, newdata = ndg))
cat("predictions se:", sprintf("%.6f", p$std.error), "\n")
