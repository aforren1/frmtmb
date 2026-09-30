# Reviewer: try to falsify the dpar-equation claims on the lane's build.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
try_msg <- function(expr) {
  r <- tryCatch(withCallingHandlers(expr,
                  warning = function(w) {
                    cat("  WARN:", conditionMessage(w), "\n")
                    invokeRestart("muffleWarning")
                  }, message = function(m) invokeRestart("muffleMessage")),
                error = function(e) structure(conditionMessage(e),
                                              class = "ERR"))
  if (inherits(r, "ERR")) cat("  ERROR:", r, "\n") else cat("  ok\n")
  invisible(r)
}
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
k <- rbinom(n, 1, 0.4)
d$y <- ifelse(k == 1, rnorm(n, 3 + 0.5 * d$x, 1), rnorm(n, -1, 1))
d$y2 <- 1 + d$x + rnorm(n)

# an independent mixture likelihood in base R, no RTMB
mixll <- function(y, mu1, mu2, s1, s2, th1) {
  sum(log(th1 * dnorm(y, mu1, s1) + (1 - th1) * dnorm(y, mu2, s2)))
}
ll_at <- function(fit) {
  s1 <- fitted(fit, dpar = "sigma1")[, "Estimate"]
  s2 <- fitted(fit, dpar = "sigma2")[, "Estimate"]
  m1 <- fitted(fit, dpar = "mu1")[, "Estimate"]
  m2 <- fitted(fit, dpar = "mu2")[, "Estimate"]
  t1 <- fitted(fit, dpar = "theta1")[, "Estimate"]
  c(indep = mixll(d$y, m1, m2, s1, s2, t1), frm = as.numeric(logLik(fit)),
    max_s1_minus_s2 = max(abs(s1 - s2)))
}

fam <- mixture(gaussian(), gaussian())
cat("E1 sigma1 = sigma2\n")
f1 <- frm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d)
print(ll_at(f1)); print(logLik(f1))
f0 <- frm(bf(y ~ x), family = fam, data = d)
cat("  unequated df", attr(logLik(f0), "df"), "equated df",
    attr(logLik(f1), "df"), "\n")
print(fixef(f1)); print(variables(f1)); print(summary(f1)$spec_pars)

cat("E2 reverse direction sigma2 = sigma1\n")
f2 <- frm(bf(y ~ x, sigma2 = "sigma1"), family = fam, data = d)
print(ll_at(f2))
cat("  rel diff logLik E1 vs E2",
    (as.numeric(logLik(f1)) - as.numeric(logLik(f2))) /
      abs(as.numeric(logLik(f1))), "\n")
print(variables(f2))

cat("E3 different links: component 1 sigma softplus, component 2 log\n")
fam_l <- mixture(brmsfamily("gaussian", link_sigma = "softplus"), gaussian())
f3 <- try_msg(frm(bf(y ~ x, sigma1 = "sigma2"), family = fam_l, data = d))
if (!inherits(f3, "ERR")) {
  print(ll_at(f3)); print(variables(f3))
  cat("  rel diff logLik vs E1",
      (as.numeric(logLik(f1)) - as.numeric(logLik(f3))) /
        abs(as.numeric(logLik(f1))), "\n")
  print(f3$frame$linpreds[["sigma1"]]$link$name)
}
f3b <- try_msg(frm(bf(y ~ x, sigma2 = "sigma1"), family = fam_l, data = d))
if (!inherits(f3b, "ERR")) {
  print(ll_at(f3b))
  print(fixef(f3b))
  cat("  rel diff logLik vs E1",
      (as.numeric(logLik(f1)) - as.numeric(logLik(f3b))) /
        abs(as.numeric(logLik(f1))), "\n")
}

cat("E4 equated dpar that also has a formula\n")
try_msg(bf(y ~ x, sigma1 = "sigma2", sigma1 ~ z))
try_msg(bf(y ~ x, sigma1 ~ z, sigma1 = "sigma2"))
try_msg(bf(y ~ x, sigma1 = "sigma2") + lf(sigma1 ~ z))
try_msg(bf(y ~ x, sigma1 ~ z) + lf(sigma1 = "sigma2"))
try_msg(bf(bf(y ~ x, sigma1 = "sigma2"), sigma1 ~ z))
try_msg(lf(sigma1 ~ z, sigma1 = "sigma2"))
try_msg(frm(bf(y ~ x, sigma1 = "sigma2") + lf(sigma1 ~ z), family = fam,
            data = d))

cat("E5 equating in a multivariate model\n")
f5 <- try_msg(frm(bf(y ~ x, sigma1 = "sigma2") + fam + bf(y2 ~ x) +
                    gaussian() + set_rescor(FALSE), data = d))
if (!inherits(f5, "ERR")) {
  g2 <- frm(y2 ~ x, data = d)
  cat("  mv logLik", as.numeric(logLik(f5)), "sum of parts",
      as.numeric(logLik(f1)) + as.numeric(logLik(g2)), "rel",
      (as.numeric(logLik(f5)) - as.numeric(logLik(f1)) -
         as.numeric(logLik(g2))) / abs(as.numeric(logLik(f5))), "\n")
  cat("  df", attr(logLik(f5), "df"), "=", attr(logLik(f1), "df"), "+",
      attr(logLik(g2), "df"), "\n")
  print(variables(f5)); print(names(f5$frame$linpreds))
  print(get_prior(bf(y ~ x, sigma1 = "sigma2") + fam + bf(y2 ~ x) +
                    gaussian() + set_rescor(FALSE), data = d))
}
cat("E5b both responses mixtures, equation only in the second\n")
d$y3 <- d$y + rnorm(n, 0, 0.1)
f5b <- try_msg(frm(bf(y ~ x) + fam + bf(y3 ~ x, sigma1 = "sigma2") + fam +
                     set_rescor(FALSE), data = d))
if (!inherits(f5b, "ERR")) {
  print(names(f5b$frame$linpreds)); print(variables(f5b))
  s <- f5b$frame$linpreds
  cat("  y3_sigma1 idx", s[["y3 sigma1"]]$idx, "y3 sigma2 idx",
      s[["y3 sigma2"]]$idx, "y sigma1 idx", s[["y sigma1"]]$idx, "\n")
}

cat("E6 chains and fan-in\n")
try_msg(bf(y ~ x, sigma1 = "sigma2", sigma2 = "sigma3"))
try_msg(bf(y ~ x, sigma2 = "sigma3", sigma1 = "sigma2"))
try_msg(bf(y ~ x, sigma1 = "sigma2") + lf(sigma2 = "sigma3"))
try_msg(bf(y ~ x, sigma1 = "sigma2") + lf(sigma3 = "sigma1"))
fam3 <- mixture(gaussian(), gaussian(), gaussian())
f6 <- try_msg(frm(bf(y ~ x, sigma1 = "sigma3", sigma2 = "sigma3"),
                  family = fam3, data = d))
if (!inherits(f6, "ERR")) {
  print(variables(f6)); print(attr(logLik(f6), "df"))
}

cat("E7 target fixed to a constant\n")
try_msg(bf(y ~ x, sigma1 = "sigma2", sigma2 = 1))
try_msg(bf(y ~ x, sigma2 = 1, sigma1 = "sigma2"))
try_msg(bf(y ~ x, sigma2 = 1) + lf(sigma1 = "sigma2"))
try_msg(update(f1, bf(~ ., sigma2 = 1)))
cat("  update(f1, bf(~ ., sigma1 = 1)): replaces the equation\n")
u7 <- try_msg(update(f1, bf(~ ., sigma1 = 1)))
if (!inherits(u7, "ERR")) print(variables(u7))

cat("E8 predictions of the equated fit\n")
try_msg(print(head(predict(f1, newdata = d[1:5, ]))))
try_msg(print(dim(simulate(f1, nsim = 3, seed = 1))))
try_msg(print(head(frm_linpred(f1, dpar = "sigma1"))))
try_msg(print(head(frm_linpred(f1, dpar = "sigma2"))))
try_msg(print(confint(f1)))
try_msg(print(hypothesis(f1, "sigma1 = sigma2")))
try_msg(print(head(fitted(f1, newdata = d[1:3, ], dpar = "sigma1"))))
try_msg(print(coef(f1)))
try_msg(print(vcov(f1)))
try_msg(print(f1))
try_msg(print(summary(f1)))

cat("E9 priors\n")
print(get_prior(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d))
try_msg(frm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d,
            prior = prior(normal(0, 1), class = "sigma1")))
fp <- try_msg(frm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d,
                  prior = prior(normal(0.2, 0.05), class = "sigma2")))
if (!inherits(fp, "ERR")) {
  print(variables(fp)); print(prior_summary(fp))
}

cat("E10 par_template and frm_simulate\n")
try_msg(print(par_template(bf(y ~ x, sigma1 = "sigma2"), family = fam,
                           data = d)))
try_msg(print(head(frm_simulate(bf(y ~ x, sigma1 = "sigma2"), family = fam,
                                data = d, seed = 2))))

cat("E12 update of the equated fit\n")
u1 <- try_msg(update(f1, newdata = d[-(1:10), ]))
if (!inherits(u1, "ERR")) print(variables(u1))
u2 <- try_msg(update(f1, . ~ . + z))
if (!inherits(u2, "ERR")) print(variables(u2))
u3 <- try_msg(update(f0, bf(~ ., sigma1 = "sigma2")))
if (!inherits(u3, "ERR")) {
  cat("  update(f0, bf(~., sigma1 = sigma2)) logLik rel diff vs f1",
      (as.numeric(logLik(u3)) - as.numeric(logLik(f1))) /
        abs(as.numeric(logLik(f1))), "\n")
}

cat("E14 student mixture nu1 = nu2 and mixed families\n")
f14 <- try_msg(frm(bf(y ~ x, nu1 = "nu2"),
                   family = mixture(student(), student()), data = d))
if (!inherits(f14, "ERR")) print(variables(f14))
try_msg(frm(bf(y ~ x, nu1 = "nu2"),
            family = mixture(gaussian(), student()), data = d))
f14c <- try_msg(frm(bf(y ~ x, sigma1 = "sigma2"),
                    family = mixture(gaussian(), student()), data = d))
if (!inherits(f14c, "ERR")) print(variables(f14c))

cat("E15 non-mixture, and categorical mu classes\n")
try_msg(bf(y ~ x, sigma = "nu"))
d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
try_msg(bf(cat ~ x, mub = "muc"))
try_msg(frm(bf(cat ~ x, mub = "muc"), family = categorical(), data = d))

cat("E16 latent: hmm transition cells\n")
if (requireNamespace("frmtmb.latent", quietly = TRUE)) {
  cat("  frmtmb.latent from", find.package("frmtmb.latent"), "\n")
}

cat("E17 bf(sigma1 = 'sigma2') where the family has neither\n")
try_msg(frm(bf(y ~ x, sigma1 = "sigma2"), data = d))
cat("DONE\n")
