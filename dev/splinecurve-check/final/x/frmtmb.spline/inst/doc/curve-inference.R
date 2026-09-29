## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.2)
has_tinyplot <- requireNamespace("tinyplot", quietly = TRUE)
library(frmtmb)
library(frmtmb.spline)

## ----data---------------------------------------------------------------------
set.seed(4)
n_sub <- 20
n_rep <- 12
n_t <- 30
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
d <- data.frame(
  subject = factor(sub),
  trial = rep(seq_len(n_sub * n_rep), each = n_t),
  t = rep(seq(0, 1, length.out = n_t), times = n_sub * n_rep))
h_sub <- rnorm(n_sub, 1, 0.12)
s_sub <- rnorm(n_sub, 0.5, 0.04)
d$v <- peak(d$t, h_sub[sub], s_sub[sub]) + rnorm(nrow(d), 0, 0.06)
str(d, max.level = 1)

## ----truth--------------------------------------------------------------------
pop_curve <- function(tt) {
  rowMeans(vapply(seq_len(n_sub),
                  function(i) peak(tt, h_sub[i], s_sub[i]),
                  numeric(length(tt))))
}
c(peak_of_one_subject = max(peak(seq(0, 1, 0.001), 1, 0.5)),
  peak_of_the_average = max(pop_curve(seq(0, 1, 0.001))))

## ----fit----------------------------------------------------------------------
fit <- frm(bf(v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)),
           family = gaussian(), data = d)
confint_varcorr(fit)

## ----curve--------------------------------------------------------------------
pop <- factor("population", levels = c(levels(d$subject), "population"))
grid <- data.frame(t = seq(0, 1, length.out = 80), subject = pop)
cv <- frm_curve(fit, newdata = grid, re_formula = NA,
                allow_new_levels = TRUE, nsim = 20000, seed = 1)
cv

## ----crit---------------------------------------------------------------------
c(pointwise = cv$.crit[1], simultaneous = cv$.crit_sim[1],
  ratio = cv$.crit_sim[1] / cv$.crit[1])

## ----mcse---------------------------------------------------------------------
attr(cv, "check")$crit_mcse

## ----check--------------------------------------------------------------------
attr(cv, "check")$cov_rel_error

## ----ncall--------------------------------------------------------------------
attr(cv, "check")$n_predict

## ----fig-curve, eval = has_tinyplot, fig.alt = "The fitted population speed profile against time, a bell shape peaking near 1.0 at t = 0.5. A wide pale band is the simultaneous interval and a narrower darker band the pointwise interval. A dashed line, the truth, lies inside both across the whole range."----
tt <- cv$t
est <- cv$.estimate
tinyplot::tinyplot(x = tt, y = est, type = "l", lwd = 2,
                   col = "steelblue4", theme = "clean2",
                   ylim = range(cv$.lower_sim, cv$.upper_sim),
                   xlab = "time", ylab = "speed",
                   main = "Population speed profile")
tinyplot::tinyplot_add(x = c(tt, rev(tt)),
                       y = c(cv$.lower_sim, rev(cv$.upper_sim)),
                       type = "polygon", col = "#4682b422", border = NA)
tinyplot::tinyplot_add(x = c(tt, rev(tt)),
                       y = c(cv$.lower_ci, rev(cv$.upper_ci)),
                       type = "polygon", col = "#4682b444", border = NA)
tinyplot::tinyplot_add(x = tt, y = est, type = "l", lwd = 2,
                       col = "steelblue4")
tinyplot::tinyplot_add(x = grid$t, y = pop_curve(grid$t), type = "l",
                       lty = 2, col = "grey30")

## ----cover--------------------------------------------------------------------
truth <- pop_curve(grid$t)
c(inside_simultaneous = mean(truth >= cv$.lower_sim & truth <= cv$.upper_sim),
  inside_pointwise = mean(truth >= cv$.lower_ci & truth <= cv$.upper_ci))

## ----deriv--------------------------------------------------------------------
g2 <- data.frame(t = seq(0.05, 0.95, length.out = 40), subject = pop)
d1 <- frm_curve_deriv(fit, var = "t", order = 1, newdata = g2,
                      re_formula = NA, allow_new_levels = TRUE,
                      nsim = 20000, seed = 2)
rising <- d1$.lower_sim > 0
falling <- d1$.upper_sim < 0
c(rising_from = min(g2$t[rising]), rising_to = max(g2$t[rising]),
  falling_from = min(g2$t[falling]), falling_to = max(g2$t[falling]))

## ----eps----------------------------------------------------------------------
c(order_1 = attr(d1, "eps"),
  order_2 = attr(frm_curve_deriv(fit, var = "t", order = 2, newdata = g2,
                                 re_formula = NA, allow_new_levels = TRUE,
                                 simultaneous = FALSE),
                 "eps"))

## ----fig-deriv, eval = has_tinyplot, fig.alt = "The first derivative of the fitted speed profile against time. It rises from zero, peaks near 3 at t = 0.35, crosses zero at t = 0.5 and falls to about minus 3 at t = 0.65 before returning toward zero. A shaded band around it excludes zero on both sides of the crossing."----
t2 <- d1$t
tinyplot::tinyplot(x = t2, y = d1$.estimate, type = "l", lwd = 2,
                   col = "firebrick", theme = "clean2",
                   ylim = range(d1$.lower_sim, d1$.upper_sim),
                   xlab = "time", ylab = "d speed / d t",
                   main = "Slope of the population profile")
tinyplot::tinyplot_add(x = c(t2, rev(t2)),
                       y = c(d1$.lower_sim, rev(d1$.upper_sim)),
                       type = "polygon", col = "#b2222222", border = NA)
tinyplot::tinyplot_add(x = t2, y = d1$.estimate, type = "l", lwd = 2,
                       col = "firebrick")
abline(h = 0, lty = 3)

## ----peak---------------------------------------------------------------------
pk <- frm_curve_feature(fit, var = "t", type = "maximum", newdata = g2,
                        re_formula = NA, allow_new_levels = TRUE)
pk

## ----onset--------------------------------------------------------------------
frm_curve_feature(fit, var = "t", type = "crossing", at = 0.2,
                  newdata = g2, re_formula = NA, allow_new_levels = TRUE)

## ----subject------------------------------------------------------------------
gs <- data.frame(t = seq(0.05, 0.95, length.out = 40),
                 subject = factor(3, levels = levels(d$subject)))
pk3 <- frm_curve_feature(fit, var = "t", type = "maximum", newdata = gs)
pk3[, c(".estimate", ".se", ".value", ".value_se")]

