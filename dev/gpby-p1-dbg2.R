# Punch round 1: the order-2 standard errors on y ~ t + (1 + t | g) at an
# unseen level, final build. The term is linear in t, so the true value
# of every order-2 .se is 0 and what is measured is rounding.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), "\n")
set.seed(23)
ng <- 25
d <- data.frame(g = factor(rep(seq_len(ng), each = 12)),
                t = rep(seq(0, 1, length.out = 12), ng))
u0 <- stats::rnorm(ng, 0, 0.5)
u1 <- stats::rnorm(ng, 0, 0.3)
d$y <- -1 + 2 * d$t + u0[d$g] + u1[d$g] * d$t +
  stats::rnorm(nrow(d), 0, 0.2)
fit <- frm(bf(y ~ t + (1 + t | g)), family = gaussian(), data = d)
nd <- data.frame(t = seq(0, 1, length.out = 21),
                 g = factor("n1", levels = c(levels(d$g), "n1")))
d2 <- frm_curve_deriv(fit, var = "t", order = 2, newdata = nd,
                      allow_new_levels = TRUE, simultaneous = FALSE)
p2 <- frm_curve_deriv(fit, var = "t", order = 2, newdata = nd,
                      re_formula = NA, simultaneous = FALSE)
d1 <- frm_curve_deriv(fit, var = "t", order = 1, newdata = nd,
                      allow_new_levels = TRUE, simultaneous = FALSE)
E2 <- frm_extra_cov_deriv(fit, nd, var = "t", order = 2,
                          allow_new_levels = TRUE)
cat(sprintf("order 2 .se, new level: max %.3e\n", max(d2$.se)))
cat(sprintf("order 2 .se, population: max %.3e\n", max(p2$.se)))
cat(sprintf("order 2 extra cov: max |E2| %.3e\n", max(abs(E2))))
cat(sprintf("order 1 .se, new level: min %.3e max %.3e\n", min(d1$.se),
            max(d1$.se)))
cat("eps order 2:", attr(d2, "eps"), "\n")
