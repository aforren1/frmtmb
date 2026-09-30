# Reviewer: does brms 2.23.0's posterior_linpred(offset = FALSE) drop
# the offset on newdata? Data seed 21 as formrobust-rev-brms-emm.R,
# Stan seed 7, 1 chain, 2000 iter.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(brms); library(emmeans)})
set.seed(21)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 1, 3),
                f = factor(sample(c("a", "b"), n, TRUE)),
                z = runif(n, 0, 1))
d$y <- rpois(n, exp(0.3 + 0.4 * d$x + 0.2 * (d$f == "b")) * d$time)
fb <- brm(y ~ x + f + offset(log(time)), data = d, family = poisson(),
          chains = 1, iter = 2000, seed = 7, refresh = 0)
nd <- data.frame(x = 0, f = factor("a", levels = c("a", "b")), time = 2)
l1 <- posterior_linpred(fb, newdata = nd, offset = TRUE)
l0 <- posterior_linpred(fb, newdata = nd, offset = FALSE)
dr <- as_draws_df(fb)
cat("offset=TRUE  - b0:", range(l1[, 1] - dr$b_Intercept), " log 2 =",
    log(2), "\n")
cat("offset=FALSE - b0:", range(l0[, 1] - dr$b_Intercept), "\n")
rg <- ref_grid(fb)
print(rg@grid)
cat("emmeans' own offset column in the grid:", ".offset." %in%
      names(rg@grid), "\n")
cat("emm (link):\n"); print(summary(emmeans(fb, ~ f)))
cat("emm at time = 1:\n"); print(summary(emmeans(fb, ~ f,
                                                 at = list(time = 1))))
# a model whose offset is a plain variable, offset(time)
fb2 <- update(fb, formula. = y ~ x + f + offset(time), newdata = d,
              refresh = 0, seed = 7)
l0 <- posterior_linpred(fb2, newdata = nd, offset = FALSE)
dr2 <- as_draws_df(fb2)
cat("offset(time), offset=FALSE - b0:", range(l0[, 1] - dr2$b_Intercept),
    "\n")
print(summary(emmeans(fb2, ~ f)))
cat("median b0, b0+bfb, bx:", median(dr2$b_Intercept),
    median(dr2$b_Intercept + dr2$b_fb), median(dr2$b_x), " mean(x)",
    mean(d$x), " mean(time)", mean(d$time), "\n")
