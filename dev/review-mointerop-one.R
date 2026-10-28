## Reviewer: a mo() term on one value (D = 0) through every consumer.
##   Rscript dev/review-mointerop-one.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
set.seed(5)
n <- 200
d <- data.frame(x = sample(1:3, n, TRUE), z = rnorm(n),
                g = factor(sample(1:8, n, TRUE)))
d$y <- 1 + 0.8 * d$x + 0.5 * d$z + rnorm(8)[d$g] * 0.5 + rnorm(n)
d$x1 <- 4L
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$z))

f1 <- show("fit y ~ mo(x1) + z", frm(y ~ mo(x1) + z, data = d))
ml <- lm(y ~ z, data = d)
show("logLik and df frm", c(ll = as.numeric(logLik(f1)),
                            df = attr(logLik(f1), "df")))
show("logLik and df lm(y ~ z)", c(ll = as.numeric(logLik(ml)),
                                  df = attr(logLik(ml), "df")))
show("lm with an aliased zero column", {
  d$zero <- 0
  m0 <- lm(y ~ zero + z, data = d)
  c(df = attr(logLik(m0), "df"), aic = AIC(m0))
})
show("AIC frm vs lm", c(frm = AIC(f1), lm = AIC(ml)))
show("BIC frm vs lm", c(frm = BIC(f1), lm = BIC(ml)))
show("summary", summary(f1))
show("fixef", fixef(f1))
show("vcov", vcov(f1))
show("confint", confint(f1))
show("predict newdata", predict(f1, newdata = data.frame(x1 = 4, z = c(0, 1))))
show("fitted newdata", fitted(f1, newdata = data.frame(x1 = 4, z = c(0, 1))))
show("lm predict se (ML sigma)", {
  s <- summary(ml)
  predict(ml, newdata = data.frame(z = c(0, 1)), se.fit = TRUE)$se.fit *
    sqrt(sum(resid(ml)^2) / n) / s$sigma
})
show("conditional_effects z", head(conditional_effects(f1, effects = "z")$z[
  , c("z", "x1", "estimate__", "se__")], 3))
show("conditional_effects x1", conditional_effects(f1, effects = "x1")$x1[
  , c("x1", "estimate__", "se__")])
show("emmeans ~ z", emmeans(f1, ~ z, at = list(z = c(0, 1))))
show("emmeans ~ x1", emmeans(f1, ~ x1))
show("emmeans lm ~ z", emmeans(ml, ~ z, at = list(z = c(0, 1))))
show("predictions", as.data.frame(predictions(
  f1, newdata = data.frame(x1 = 4, z = c(0, 1))))[, c("estimate", "std.error")])
show("avg_slopes z", avg_slopes(f1, variables = "z"))
show("avg_slopes lm z", avg_slopes(ml, variables = "z"))
show("predictions HC0", as.data.frame(predictions(
  f1, newdata = data.frame(x1 = 4, z = c(0, 1)), vcov = "HC0"))[
    , c("estimate", "std.error")])
show("predictions lm HC0", as.data.frame(predictions(
  ml, newdata = data.frame(z = c(0, 1)), vcov = "HC0"))[
    , c("estimate", "std.error")])
show("predictions ~ g", as.data.frame(predictions(
  f1, newdata = data.frame(x1 = 4, z = c(0, 1)), vcov = ~ g))[
    , c("estimate", "std.error")])
show("predictions vcov(full = TRUE)", as.data.frame(predictions(
  f1, newdata = data.frame(x1 = 4, z = c(0, 1)),
  vcov = vcov(f1, full = TRUE)))[, c("estimate", "std.error")])
show("vcov_cluster CR0 each row", vcov_cluster(f1, seq_len(n), "CR0"))
show("hypothesis mox1 = 0", hypothesis(f1, "mox1 = 0"))
show("hypothesis z = 0", hypothesis(f1, "z = 0"))
f0 <- frm(y ~ z, data = d)
show("anova(f0, f1)", anova(f0, f1))
show("drop1(f1)", drop1(f1))
show("update - mo(x1)", fixef(update(f1, . ~ . - mo(x1))))
show("update newdata 3 values", fixef(update(f1, newdata = transform(
  d, x1 = x))))
tf <- tempfile(fileext = ".rds")
saveRDS(f1, tf)
f1r <- readRDS(tf)
show("saved fit fitted", fitted(f1r, newdata = data.frame(x1 = 4, z = 0)))
show("saved fit summary", summary(f1r))
show("one-value interaction mo(x1) * z", {
  f <- frm(y ~ mo(x1) * z, data = d)
  list(fixef = fixef(f), se = sqrt(diag(vcov(f))))
})
show("one-value with random slope (mo(x1) | g)", {
  f <- frm(y ~ mo(x1) + z + (mo(x1) | g), data = d)
  list(fixef = fixef(f), VarCorr = VarCorr(f), ll = logLik(f))
})
show("one-value in sigma", {
  f <- frm(bf(y ~ z, sigma ~ mo(x1)), data = d)
  list(fixef = fixef(f), ll = logLik(f))
})
show("one-value poisson", {
  f <- frm(cnt ~ mo(x1) + z, data = d, family = poisson())
  list(fixef = fixef(f), ll = logLik(f),
       glm = logLik(glm(cnt ~ z, data = d, family = poisson())))
})
show("one-value with prior on b", {
  f <- frm(y ~ mo(x1) + z, data = d,
           prior = set_prior("normal(0, 1)", class = "b"))
  list(fixef = fixef(f), se = sqrt(diag(vcov(f))))
})
show("one-value ordered factor with one level", {
  dd <- d
  dd$xf <- factor(rep("a", n), levels = "a", ordered = TRUE)
  fixef(frm(y ~ mo(xf) + z, data = dd))
})
show("frm_sample refusal", {
  suppressMessages(library(frmtmb.sample))
  frm_sample(f1, chains = 1, iter = 10)
})
show("frm_sample from formula", frmtmb.sample::frm_sample(
  y ~ mo(x1) + z, data = d, chains = 1, iter = 10))
