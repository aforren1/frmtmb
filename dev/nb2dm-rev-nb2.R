# Reviewer probes: negbinomial2 through the post-fit paths, side by side
# with negbinomial(). Seed 2026 (data), n = 400, 20 groups.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
set.seed(2026)
n <- 400
dd <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(20, n / 20),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)))
dd$e <- exp(runif(n))
u <- rnorm(20, 0, 0.4)[dd$g]
dd$y <- rnbinom(n, mu = dd$e * exp(0.7 + 0.5 * dd$x + 0.3 * (dd$f == "b") + u),
                size = 1 / exp(-0.9 + 0.4 * dd$z))
NB2 <- brmsfamily("negbinomial2")

f1 <- frm(y ~ x + f, family = negbinomial(), data = dd)
f2 <- frm(y ~ x + f, family = NB2, data = dd)
probe("logLik nb / nb2", c(logLik(f1), logLik(f2)))
probe("summary nb2", summary(f2))
probe("print nb2", print(f2))
probe("confint wald nb2", confint(f2))
probe("confint wald nb", confint(f1))
probe("confint profile nb2", confint(f2, method = "profile"))
probe("confint profile nb", confint(f1, method = "profile"))
probe("profile nb2 sigma names", names(profile(f2, parm = "sigma")))
probe("sigma(nb2) and 1/exp(b_shape)",
      c(sigma(f2), 1 / exp(fixef_by_dpar(f1)$shape)))
probe("sigma(nb) (negbinomial: what does sigma() say)", sigma(f1))
probe("predict type disp nb2 head", head(predict(f2, type = "disp")))
probe("residuals pearson nb2 vs nb, max rel",
      max(abs(residuals(f2, type = "pearson")[, 1] /
                residuals(f1, type = "pearson")[, 1] - 1)))
probe("residuals quantile nb2", head(residuals(f2, type = "quantile")))
probe("residuals quantile nb", head(residuals(f1, type = "quantile")))
probe("dharma nb2", class(dharma_residuals(f2, nsim = 20)))
probe("residuals osa nb2", head(residuals(f2, type = "osa")))
probe("residuals osa nb", head(residuals(f1, type = "osa")))
probe("conditional_effects nb2", {
  ce <- conditional_effects(f2, effects = "x")
  head(ce[[1]])
})
probe("conditional_effects nb2 dpar sigma", {
  ce <- conditional_effects(f2, effects = "x", dpar = "sigma")
  head(ce[[1]])
})
probe("emmeans nb2 vs nb", {
  e2 <- summary(emmeans::emmeans(f2, ~ f, type = "response"))
  e1 <- summary(emmeans::emmeans(f1, ~ f, type = "response"))
  list(e2 = e2, max_rel = max(abs(e2$response / e1$response - 1)))
})
probe("emmeans nb2 dpar sigma", {
  summary(emmeans::emmeans(f2, ~ 1, dpar = "sigma"))
})
probe("marginaleffects avg_slopes nb2 vs nb", {
  s2 <- marginaleffects::avg_slopes(f2, variables = "x")
  s1 <- marginaleffects::avg_slopes(f1, variables = "x")
  c(nb2 = s2$estimate, nb = s1$estimate, se2 = s2$std.error,
    se1 = s1$std.error)
})
probe("simulate nb2 moments vs fitted", {
  s <- simulate(f2, nsim = 200, seed = 1)
  c(mean_sim = mean(as.matrix(s)), mean_fit = mean(fitted(f2)[, 1]))
})
probe("anova nb vs nb2 (same model)", anova(f1, f2))
probe("drop1 nb2", drop1(f2))
probe("update nb2 ~ . - f", logLik(update(f2, . ~ . - f)))

# rate(): fitted and predict on newdata
r1 <- frm(y | rate(e) ~ x + f, family = negbinomial(), data = dd)
r2 <- frm(y | rate(e) ~ x + f, family = NB2, data = dd)
nd <- data.frame(x = c(-1, 0, 1), f = c("a", "b", "c"), e = c(0.5, 1, 3))
probe("rate fitted newdata nb2 / nb",
      fitted(r2, newdata = nd)[, 1] / fitted(r1, newdata = nd)[, 1])
probe("rate predict newdata nb2", {
  set.seed(5)
  predict(r2, newdata = nd, ndraws = 4000)
})
probe("rate predict newdata nb", {
  set.seed(5)
  predict(r1, newdata = nd, ndraws = 4000)
})
probe("rate fitted newdata without e", fitted(r2, newdata = nd[, 1:2]))

# sigma ~ z and its consumers
s2 <- frm(bf(y ~ x, sigma ~ z), family = NB2, data = dd)
s1 <- frm(bf(y ~ x, shape ~ z), family = negbinomial(), data = dd)
probe("sigma ~ z: logLik", c(logLik(s1), logLik(s2)))
probe("sigma ~ z: fixef sigma vs -shape",
      rbind(fixef_by_dpar(s2)$sigma, -fixef_by_dpar(s1)$shape))
probe("sigma ~ z: fitted dpar sigma head",
      head(fitted(s2, dpar = "sigma")))
probe("sigma ~ z: simulate", dim(simulate(s2, nsim = 2, seed = 1)))
probe("shape ~ z under negbinomial2 (must refuse)",
      frm(bf(y ~ x, shape ~ z), family = NB2, data = dd))
probe("sigma ~ z under negbinomial (must refuse)",
      frm(bf(y ~ x, sigma ~ z), family = negbinomial(), data = dd))

# group term: VarCorr and summary
g2 <- frm(y ~ x + (1 | g), family = NB2, data = dd)
probe("VarCorr nb2", VarCorr(g2))
probe("ranef nb2 head", head(ranef(g2)$g[, , 1]))

# priors on sigma, the core (MAP) route
probe("default_prior nb2 (fit route)", default_prior(y ~ x, dd, family = NB2))
probe("frm with prior on sigma", {
  p <- frm(y ~ x, family = NB2, data = dd,
           prior = set_prior("gamma(2, 4)", class = "sigma"))
  c(sigma_prior = sigma(p), sigma_ml = sigma(frm(y ~ x, family = NB2,
                                                    data = dd)))
})
probe("frm with prior on shape under nb2 (must refuse)",
      frm(y ~ x, family = NB2, data = dd,
          prior = set_prior("gamma(2, 4)", class = "shape")))

# mixture and mvbf
probe("mixture(negbinomial2, poisson)", {
  m <- frm(y ~ x, family = mixture(NB2, poisson()), data = dd)
  list(ll = logLik(m), fx = fixef_by_dpar(m))
})
probe("mvbf(nb2, gaussian)", {
  dd$w <- 0.5 * dd$x + rnorm(n)
  m <- frm(mvbf(bf(y ~ x, family = NB2), bf(w ~ x, family = gaussian())),
           data = dd)
  logLik(m)
})
# Poisson data: sigma -> 0 is the boundary; what do the messages say?
set.seed(7)
dp <- data.frame(x = rnorm(300))
dp$y <- rpois(300, exp(1 + 0.3 * dp$x))
probe("poisson data, negbinomial", {
  m <- frm(y ~ x, family = negbinomial(), data = dp)
  list(conv = m$opt$convergence, shape = exp(fixef_by_dpar(m)$shape),
       se = sqrt(diag(vcov(m))))
})
probe("poisson data, negbinomial2", {
  m <- frm(y ~ x, family = NB2, data = dp)
  list(conv = m$opt$convergence, sigma = sigma(m), se = sqrt(diag(vcov(m))),
       ll = logLik(m), ll_pois = logLik(glm(y ~ x, poisson, dp)))
})
probe("poisson data, negbinomial2 summary", summary(frm(y ~ x, family = NB2,
                                                         data = dp)))
# a non-count response: the message names which family?
probe("non-integer response", frm(I(y + 0.5) ~ x, family = NB2, data = dd))
probe("negative response", frm(I(y - 100) ~ x, family = NB2, data = dd))
# other links on sigma
probe("link_sigma identity fit", {
  m <- frm(y ~ x, family = brmsfamily("negbinomial2", link_sigma = "identity"),
           data = dd)
  c(sigma(m), logLik(m), logLik(f2))
})
probe("link_sigma softplus with sigma ~ z", {
  m <- frm(bf(y ~ x, sigma ~ z),
           family = brmsfamily("negbinomial2", link_sigma = "softplus"),
           data = dd)
  logLik(m)
})
probe("frm_compat negbinomial2 rows", {
  tb <- frm_compat("negbinomial2")
  tb[, c("feature_a", "feature_b", "status")]
})
