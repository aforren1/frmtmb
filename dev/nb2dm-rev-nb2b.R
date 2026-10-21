# Reviewer probes, second batch: mixture and profile paths of
# negbinomial2 against negbinomial (same data as dev/nb2dm-rev-nb2.R).
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
probe("mixture(nb2, poisson) vs mixture(negbinomial, poisson)", {
  m2 <- frm(y ~ x, family = mixture(NB2, poisson()), data = dd)
  m1 <- frm(y ~ x, family = mixture(negbinomial(), poisson()), data = dd)
  c(nb2 = logLik(m2), nb = logLik(m1),
    sigma1 = exp(fixef_by_dpar(m2)$sigma1[[1]]),
    inv_shape1 = 1 / exp(fixef_by_dpar(m1)$shape1[[1]]))
})
f1 <- frm(y ~ x + f, family = negbinomial(), data = dd)
f2 <- frm(y ~ x + f, family = NB2, data = dd)
probe("confint profile sigma_(Intercept) vs -shape_(Intercept)", {
  c2 <- confint(f2, parm = "sigma_(Intercept)", method = "profile")
  c1 <- confint(f1, parm = "shape_(Intercept)", method = "profile")
  list(nb2 = c2, nb = c1)
})
probe("hypothesis sigma = 0.5", hypothesis(f2, "sigma = 0.5", class = NULL))
probe("emmeans response column names", {
  e2 <- summary(emmeans::emmeans(f2, ~ f, type = "response"))
  e1 <- summary(emmeans::emmeans(f1, ~ f, type = "response"))
  list(names2 = names(e2), names1 = names(e1),
       max_rel = max(abs(e2[[2]] / e1[[2]] - 1)))
})
probe("quadrature (1 | g) nb2 vs nb", {
  q2 <- frm(y ~ x + (1 | g), family = NB2, data = dd, quadrature = TRUE)
  q1 <- frm(y ~ x + (1 | g), family = negbinomial(), data = dd,
            quadrature = TRUE)
  c(nb2 = logLik(q2), nb = logLik(q1))
})
probe("VarCorr residual__ equals sigma()", {
  g2 <- frm(y ~ x + (1 | g), family = NB2, data = dd)
  c(varcorr = VarCorr(g2)$residual__$sd[1, 1], sigma = sigma(g2))
})
probe("VarCorr of negbinomial with (1 | g): residual__?", {
  g1 <- frm(y ~ x + (1 | g), family = negbinomial(), data = dd)
  names(VarCorr(g1))
})
probe("valid_y message on nb2 vs brms wording", {
  tryCatch(frm(I(y + 0.5) ~ x, family = NB2, data = dd),
           error = conditionMessage)
})
