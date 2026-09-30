# Lane wt-defects: the two identities dev/defects-findings.md quotes,
# on the constructions of test-brms-parity-defects.R.
#   Rscript dev/defects-probe-identities.R
source("dev/defects-pre.R")
# 1. fitted(scale = "linear") on sratio + cs(): P(Y = 1) = F(tau_1 -
#    eta_1), P(Y = 2) = (1 - P(Y = 1)) F(tau_2 - eta_2)
set.seed(11)
n <- 150
dd <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
y <- ifelse(runif(n) < plogis(-0.5 - 0.8 * dd$x1 - 0.3 * dd$x2), 1L, NA)
p2 <- plogis(0.7 - 0.8 * dd$x1 + 0.4 * dd$x2)
y[is.na(y)] <- ifelse(runif(sum(is.na(y))) < p2[is.na(y)], 2L, 3L)
dd$y <- y
fit <- frm(bf(y ~ x1 + cs(x2)), family = sratio(), data = dd)
fl <- fitted(fit, scale = "linear")
fr <- fitted(fit)
tau <- fixef(fit)[c("Intercept[1]", "Intercept[2]"), "Estimate"]
q1 <- plogis(tau[[1]] - fl[, "Estimate", 1])
q2 <- (1 - q1) * plogis(tau[[2]] - fl[, "Estimate", 2])
cat(sprintf("sratio P1 residual %.3g, P2 residual %.3g\n",
            max(abs(q1 - fr[, "Estimate", 1])),
            max(abs(q2 - fr[, "Estimate", 2]))))
# 2. summary()$gp's lscale times the largest distance is
#    confint_varcorr()'s range in data units
set.seed(12)
n <- 40
dg <- data.frame(x = seq(0, 10, length.out = n))
dg$y <- sin(dg$x) + rnorm(n, 0, 0.3)
fg <- frm(y ~ gp(x), data = dg)
ls <- summary(fg)$gp["lscale(gpx)", "Estimate"]
cv <- confint_varcorr(fg)
rg <- cv$estimate[cv$term == "range(gp)"]
cat(sprintf("lscale %.7f, x 10 = %.7f, range %.7f, difference %.3g\n",
            ls, 10 * ls, rg, 10 * ls - rg))
