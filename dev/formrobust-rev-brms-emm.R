# Reviewer: brms 2.23.0's emmeans() and conditional_effects() on
# offset models, fitted (Stan seed 7, 1 chain, 2000 iter), against the
# lane's claim that emmeans leaves the offset out and epred = TRUE
# puts it in at mean(time). Data seed 21.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(brms); library(emmeans)})
set.seed(21)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 1, 3),
                f = factor(sample(c("a", "b"), n, TRUE)),
                z = runif(n, 0, 1))
d$y <- rpois(n, exp(0.3 + 0.4 * d$x + 0.2 * (d$f == "b")) * d$time)
cat("allvars:", deparse(brms:::brmsterms(bf(y ~ x + f +
                                             offset(log(time))))$allvars),
    "\n")
fb <- brm(y ~ x + f + offset(log(time)), data = d, family = poisson(),
          chains = 1, iter = 2000, seed = 7, refresh = 0)
dr <- as_draws_df(fb)
em <- summary(emmeans(fb, ~ f))
em_ep <- summary(emmeans(fb, ~ f, epred = TRUE))
b0 <- dr$b_Intercept; bx <- dr$b_x; bf_ <- dr$b_fb
mx <- mean(d$x)
lp_a <- b0 + bx * mx
lp_b <- b0 + bx * mx + bf_
cat("emmean a, b:", em$emmean, "\n")
cat("median of draws b0 + bx*mean(x) (+ bfb):", median(lp_a), median(lp_b),
    "\n")
cat("difference a (0 means offset left out):", em$emmean[1] - median(lp_a),
    " log(mean(time)) =", log(mean(d$time)), "\n")
# the grid holds time at its mean; epred at that grid
ep_a <- exp(lp_a + log(mean(d$time)))
ep_b <- exp(lp_b + log(mean(d$time)))
cat("epred emmean a, b:", em_ep$emmean, "\n")
cat("median exp(lp + log(mean(time))):", median(ep_a), median(ep_b), "\n")
ce <- conditional_effects(fb, effects = "x")[[1]]
i <- which.min(abs(ce$x - 0))
cat("CE at x =", ce$x[i], " estimate", ce$estimate__[i], " time in grid",
    unique(ce$time), " mean(time)", mean(d$time), "\n")
cat("CE / median exp(b0 + bx x + log(mean time)) (f at ref a):",
    ce$estimate__[i] / median(exp(b0 + bx * ce$x[i] + log(mean(d$time)))),
    "\n")

# a nonlinear model with an offset in a nonlinear parameter's formula
d$y2 <- rnorm(n, 1 + 0.5 * d$x + d$z, 0.5)
fnl <- brm(bf(y2 ~ a + b * x, a ~ 1 + offset(z), b ~ 1, nl = TRUE),
           data = d, family = gaussian(), chains = 1, iter = 2000, seed = 7,
           refresh = 0,
           prior = c(prior(normal(0, 5), nlpar = "a"),
                     prior(normal(0, 5), nlpar = "b")))
drn <- as_draws_df(fnl)
em2 <- summary(emmeans(fnl, ~ 1))
cat("nl emmean:", em2$emmean, " median(a + b*mean(x)):",
    median(drn$b_a_Intercept + drn$b_b_Intercept * mx),
    " + mean(z) would add", mean(d$z), "\n")
# a sigma offset
d$y3 <- rnorm(n, 0.5 * d$x, exp(0.2 + log(d$time)))
fs <- brm(bf(y3 ~ x, sigma ~ 1 + offset(log(time))), data = d,
          chains = 1, iter = 2000, seed = 7, refresh = 0)
drs <- as_draws_df(fs)
em3 <- summary(emmeans(fs, ~ 1, dpar = "sigma"))
cat("sigma emmean:", em3$emmean, " median(b_sigma_Intercept):",
    median(drs$b_sigma_Intercept), " log(mean(time))", log(mean(d$time)),
    "\n")
saveRDS(list(d = d, pois = fixef(fb), nl = fixef(fnl), sig = fixef(fs)),
        "C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-log/brms-emm.rds")
