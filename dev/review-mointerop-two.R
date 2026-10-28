## Reviewer: two-value mo() against brms 2.23.0 (no compile: standata,
## default_prior, stancode) and frmtmb consumers.
##   Rscript dev/review-mointerop-two.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
suppressMessages(library(brms))
set.seed(5)
n <- 200
d <- data.frame(x = sample(1:2, n, TRUE), z = rnorm(n),
                g = factor(sample(1:8, n, TRUE)))
d$y <- 1 + 0.8 * d$x + 0.5 * d$z + rnorm(8)[d$g] * 0.5 + rnorm(n)
d$x1 <- 4L
f2 <- frm(y ~ mo(x) + z, data = d)
show("frm default_prior two", as.data.frame(frmtmb::default_prior(
  y ~ mo(x) + z, data = d))[, c("prior", "class", "coef")])
show("brms default_prior two", as.data.frame(brms::default_prior(
  y ~ mo(x) + z, data = d))[, c("prior", "class", "coef")])
show("frm default_prior one", as.data.frame(frmtmb::default_prior(
  y ~ mo(x1) + z, data = d))[, c("prior", "class", "coef")])
show("brms default_prior one", as.data.frame(brms::default_prior(
  y ~ mo(x1) + z, data = d))[, c("prior", "class", "coef")])
show("frm default_prior two, interaction", as.data.frame(
  frmtmb::default_prior(y ~ mo(x) * z, data = d))[, c("prior", "class",
                                                       "coef")])
show("brms default_prior two, interaction", as.data.frame(
  brms::default_prior(y ~ mo(x) * z, data = d))[, c("prior", "class",
                                                     "coef")])
show("frm variables", variables(f2))
show("frm prior_summary", prior_summary(f2))
show("frm summary", summary(f2))
show("frm ce x", conditional_effects(f2, effects = "x")$x[, c("x", "z",
                                                             "estimate__")])
show("frm ce z held x", unique(conditional_effects(f2, effects = "z")$z$x))
## data where the median of x is 2 (brms holds at round(median()))
d2 <- d
d2$x <- ifelse(runif(n) < 0.8, 2L, 1L)
f2b <- frm(y ~ mo(x) + z, data = d2)
show("frm ce z held x, median 2", unique(conditional_effects(
  f2b, effects = "z")$z$x))
show("new data x = 3", fitted(f2, newdata = data.frame(x = 3, z = 0)))
show("new data x = 0", fitted(f2, newdata = data.frame(x = 0, z = 0)))
show("new data x = 1.5", fitted(f2, newdata = data.frame(x = 1.5, z = 0)))
show("emmeans ~ x", emmeans(f2, ~ x, at = list(x = 1:2)))
show("emmeans default (no at)", emmeans(f2, ~ x))
show("pairs emmeans", pairs(emmeans(f2, ~ x, at = list(x = 1:2))))
show("lm pairs", pairs(emmeans(lm(y ~ factor(x) + z, data = d), ~ x)))
show("predictions x 1, 2", predictions(f2, newdata = data.frame(
  x = 1:2, z = 0))[, c("estimate", "std.error")])
show("anova vs linear", anova(frm(y ~ z, data = d), f2))
show("drop1", drop1(f2))
show("update to factor", fixef(update(f2, . ~ . - mo(x) + factor(x))))
tf <- tempfile(fileext = ".rds")
saveRDS(f2, tf)
show("saved fit summary", summary(readRDS(tf)))
show("hypothesis bsp", hypothesis(f2, "mox = 0"))
## brms standata on 2 categories with a factor multiplier
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
show("brms standata Jmo mo(x) * f", brms::standata(y ~ mo(x) * f, data = d)[
  c("Jmo", "Imo")])
f2f <- show("frm mo(x) * f two values", frm(y ~ mo(x) * f, data = d))
show("frm fixef mo(x) * f", fixef(f2f))
show("lm I(x-1) * f", coef(lm(y ~ I(x - 1) * f, data = d)))
show("frm logLik", as.numeric(logLik(f2f)))
show("lm logLik", as.numeric(logLik(lm(y ~ I(x - 1) * f, data = d))))
show("summary mo(x) * f", summary(f2f))
show("frm variables mo(x)*f", variables(f2f))
## brms default_prior simo names for mo(x)*f
show("brms default_prior mo(x) * f", as.data.frame(brms::default_prior(
  y ~ mo(x) * f, data = d))[, c("prior", "class", "coef")])
show("frm default_prior mo(x) * f", as.data.frame(frmtmb::default_prior(
  y ~ mo(x) * f, data = d))[, c("prior", "class", "coef")])
