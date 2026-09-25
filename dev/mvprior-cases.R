# Lane wt-mvprior: the models and prior specifications both probes run,
# so brms 2.23.0 and frmtmb are measured on the same rows and the same
# specifications. A specification is a plain list of set_prior()
# arguments; each probe builds it with its own set_prior().

mvprior_data <- function() {
  set.seed(2309)
  n <- 80
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(rep(1:8, length.out = n)),
                  t = rep(1:10, each = 8))
  d$y1 <- 1 + d$x + rnorm(8)[d$g] + rnorm(n)
  d$y2 <- -1 + 0.5 * d$x + rnorm(8)[d$g] + rnorm(n, 0, exp(0.3 * d$z))
  d$cnt <- stats::rpois(n, exp(0.5 + 0.3 * d$x))
  d$y <- d$y1
  d$ynl <- (2 + rnorm(8, 0, 0.3)[d$g]) * exp(-0.5 * abs(d$x)) +
    rnorm(n, 0, 0.1)
  d$ax <- abs(d$x)
  # a short panel for unstr(): 20 units of 4 occasions
  d$t4 <- rep(1:4, length.out = n)
  d$g4 <- factor(rep(1:20, each = 4))
  # a three-category response and a bimodal one, for the families whose
  # location is several dpars (punch round 1, M1)
  d$cat <- factor(c("a", "b", "c")[1 + (d$x + stats::rlogis(n) > 0) +
                                      (d$z + stats::rlogis(n) > 0.5)])
  d$ym <- ifelse(stats::rbinom(n, 1, 0.5) == 1, 3 + d$x, -1 + d$x) +
    stats::rnorm(n, 0, 0.5)
  d
}

# one specification, as set_prior() arguments
sp <- function(class = "b", coef = "", group = "", resp = "", dpar = "",
               nlpar = "", prior = "normal(0, 5)", lb = NA) {
  list(prior = prior, class = class, coef = coef, group = group,
       resp = resp, dpar = dpar, nlpar = nlpar, lb = lb)
}

# `fam` names the family spelling: each probe maps it to its own
# objects. `f` is built by the probe from the same pieces.
mvprior_cases <- list(
  mv2 = list(
    what = "two-response gaussian, bf(y1 ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)), rescor FALSE",
    specs = list(
      sp("b"), sp("b", coef = "x"), sp("b", resp = "y1"),
      sp("b", coef = "x", resp = "y1"),
      sp("Intercept"), sp("Intercept", resp = "y2"),
      sp("sigma", prior = "student_t(3, 0, 2.5)"),
      sp("sigma", resp = "y1", prior = "student_t(3, 0, 2.5)"),
      sp("sd"), sp("sd", resp = "y1"),
      sp("b", prior = "", lb = 0), sp("b", resp = "y1", prior = "", lb = 0),
      sp("b", resp = "nosuch"))),
  mvsig = list(
    what = "bf(y1 ~ x, sigma ~ z) + bf(y2 ~ x), rescor FALSE",
    specs = list(
      sp("b", dpar = "sigma"), sp("b", dpar = "sigma", resp = "y1"),
      sp("Intercept", dpar = "sigma"),
      sp("Intercept", dpar = "sigma", resp = "y1"),
      sp("sigma", prior = "student_t(3, 0, 2.5)"),
      sp("sigma", resp = "y2", prior = "student_t(3, 0, 2.5)"),
      sp("b"), sp("b", resp = "y2"))),
  mixed = list(
    what = "bf(y1 ~ x) + gaussian() + bf(cnt ~ x) + poisson()",
    specs = list(
      sp("b"), sp("b", resp = "cnt"), sp("Intercept"),
      sp("Intercept", resp = "y1"),
      sp("sigma", prior = "student_t(3, 0, 2.5)"),
      sp("sigma", resp = "y1", prior = "student_t(3, 0, 2.5)"),
      sp("sigma", resp = "cnt", prior = "student_t(3, 0, 2.5)"))),
  rescor = list(
    what = "bf(mvbind(y1, y2) ~ x + (1 + x | g)) + set_rescor(TRUE)",
    specs = list(
      sp("b"), sp("b", resp = "y1"), sp("Intercept"),
      sp("Intercept", resp = "y1"),
      sp("sigma", prior = "student_t(3, 0, 2.5)"),
      sp("sigma", resp = "y2", prior = "student_t(3, 0, 2.5)"),
      sp("rescor", prior = "lkj(2)"),
      sp("rescor", resp = "y1", prior = "lkj(2)"),
      sp("cor", prior = "lkj(2)"), sp("cor", group = "g", prior = "lkj(2)"),
      sp("cor", resp = "y1", prior = "lkj(2)"))),
  mvstudent = list(
    what = "bf(y1 ~ x) + bf(y2 ~ x), student(), rescor FALSE",
    specs = list(
      sp("nu", prior = "gamma(2, 0.1)"),
      sp("nu", resp = "y1", prior = "gamma(2, 0.1)"))),
  mvar = list(
    what = "bf(y1 ~ x + ar(time = t, gr = g, cov = TRUE)) + bf(y2 ~ x), rescor FALSE",
    specs = list(
      sp("ar", prior = "normal(0, 0.5)"),
      sp("ar", resp = "y1", prior = "normal(0, 0.5)"),
      sp("ar", resp = "y2", prior = "normal(0, 0.5)"))),
  nl = list(
    what = "bf(ynl ~ a * exp(-b * ax), a ~ 1 + (1 | g), b ~ 1, nl = TRUE)",
    specs = list(
      sp("b"), sp("b", nlpar = "a"), sp("b", coef = "Intercept"),
      sp("b", coef = "Intercept", nlpar = "b"),
      sp("Intercept"),
      sp("sigma", prior = "student_t(3, 0, 2.5)"),
      sp("b", prior = "", lb = 0), sp("b", nlpar = "b", prior = "", lb = 0),
      sp("Intercept", nlpar = "a"))),
  uni = list(
    what = "univariate: bf(y ~ x + (1 | g), sigma ~ z)",
    specs = list(
      sp("b"), sp("b", coef = "x"), sp("Intercept"), sp("b", resp = "y"),
      sp("Intercept", resp = "y"), sp("b", resp = "nosuch"),
      sp("b", dpar = "sigma"), sp("Intercept", dpar = "sigma"),
      sp("sd"))),
  uni2 = list(
    what = "univariate: bf(y ~ x)",
    specs = list(
      sp("sigma", prior = "student_t(3, 0, 2.5)"),
      sp("sigma", resp = "y", prior = "student_t(3, 0, 2.5)"),
      sp("b", prior = "", lb = 0)))
)

mvprior_cases$mi <- list(
  what = "bf(y ~ mi(xm) + z) + bf(xm | mi() ~ z), rescor FALSE",
  specs = list(
    sp("b"), sp("b", resp = "y"), sp("b", coef = "mixm"),
    sp("b", coef = "mixm", resp = "y"), sp("Intercept"),
    sp("sigma", prior = "student_t(3, 0, 2.5)")))

mvprior_cases$mvac <- list(
  what = paste("bf(y1 ~ x + ma(time = t, gr = g, cov = TRUE)) +",
               "bf(y2 ~ x + cosy(time = t, gr = g)), rescor FALSE"),
  specs = list(
    sp("ma", prior = "normal(0, 0.5)"),
    sp("ma", resp = "y1", prior = "normal(0, 0.5)"),
    sp("cosy", prior = "beta(2, 2)"),
    sp("cosy", resp = "y2", prior = "beta(2, 2)")))

mvprior_cases$mvunstr <- list(
  what = "bf(y1 ~ x + unstr(time = t4, gr = g4)) + bf(y2 ~ x), rescor FALSE",
  specs = list(
    sp("cortime", prior = "lkj(2)"),
    sp("cortime", resp = "y1", prior = "lkj(2)")))

mvprior_cases$uni3 <- list(
  what = "univariate, no slope: bf(y ~ 1)",
  specs = list(sp("b"), sp("Intercept")))

mvprior_cases$cat <- list(
  what = "univariate categorical: bf(cat ~ x), categorical()",
  specs = list(
    sp("b"), sp("b", coef = "x"), sp("Intercept"),
    sp("b", dpar = "mub"), sp("b", coef = "x", dpar = "muc"),
    sp("Intercept", dpar = "mub"), sp("b", dpar = "nosuch")))

mvprior_cases$mix <- list(
  what = "mixture: bf(ym ~ x), mixture(gaussian, gaussian)",
  specs = list(
    sp("b"), sp("b", coef = "x"), sp("Intercept"),
    sp("b", dpar = "mu1"), sp("Intercept", dpar = "mu2"),
    sp("sigma1", prior = "student_t(3, 0, 2.5)"),
    sp("sigma", prior = "student_t(3, 0, 2.5)")))

mvprior_cases$mvcat <- list(
  what = "bf(y1 ~ x) + bf(cat ~ x, family = categorical()), rescor FALSE",
  specs = list(
    sp("b", resp = "cat"), sp("Intercept", resp = "cat"),
    sp("b", resp = "y1"),
    sp("b", dpar = "mub", resp = "cat"), sp("b", dpar = "mub"),
    sp("Intercept", dpar = "muc", resp = "cat"),
    sp("Intercept", dpar = "muc")))

mvprior_cases$mvnl <- list(
  what = paste("bf(y1 ~ x) + bf(ynl ~ a * exp(-b * ax), a ~ 1 + z,",
               "b ~ 1, nl = TRUE), rescor FALSE"),
  specs = list(
    sp("b", nlpar = "a"), sp("b", nlpar = "a", resp = "ynl"),
    sp("Intercept", nlpar = "a"), sp("Intercept", nlpar = "a", resp = "ynl"),
    sp("b", resp = "ynl"), sp("b", coef = "z", nlpar = "a")))

mvprior_cases$nl2 <- list(
  what = "bf(ynl ~ a * exp(-b * ax), a ~ 1 + z, b ~ 1, nl = TRUE)",
  specs = list(sp("b"), sp("b", prior = "", lb = 0),
               sp("b", coef = "z", nlpar = "a")))

mvprior_cases$catre <- list(
  what = "categorical with a group effect: bf(cat ~ x + (1 | g))",
  specs = list(sp("sd"), sp("sd", group = "g"), sp("sd", dpar = "mub"),
               sp("b", dpar = "muc")))

mvprior_cases$mixre <- list(
  what = "mixture with a group effect: bf(ym ~ x + (1 | g))",
  specs = list(sp("sd"), sp("sd", dpar = "mu1")))

# the user's decisions of 2026-09-24: no slope, univariate resp, rescor
mvprior_cases$noslope <- list(
  what = "no slope: bf(y ~ 1 + (1 | g), sigma ~ 1)",
  specs = list(sp("b"), sp("b", prior = "", lb = 0),
               sp("b", dpar = "sigma"), sp("Intercept"),
               sp("sd", resp = "y"), sp("sd")))

mvprior_cases$mvnoslope <- list(
  what = "bf(y1 ~ x) + bf(y2 ~ 1), rescor FALSE",
  specs = list(sp("b", resp = "y2"), sp("b", resp = "y1"),
               sp("Intercept", resp = "y2")))

mvprior_cases$catnoslope <- list(
  what = "bf(cat ~ 1), categorical()",
  specs = list(sp("b", dpar = "mub"), sp("Intercept", dpar = "mub")))

mvprior_cases$uniar <- list(
  what = "univariate autocorrelation: bf(y ~ x + ar(time = t, gr = g, cov = TRUE))",
  specs = list(sp("ar", prior = "normal(0, 0.5)"),
               sp("ar", resp = "y", prior = "normal(0, 0.5)")))

mvprior_cases$rescor2 <- list(
  what = "bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE)",
  specs = list(sp("rescor", prior = "lkj(2)"),
               sp("rescor", resp = "y1", prior = "lkj(2)"),
               sp("rescor", resp = "y2", prior = "lkj(2)")))

mvprior_data_mi <- function() {
  set.seed(108)
  n <- 80
  d <- data.frame(z = stats::rnorm(n))
  d$xm <- 0.5 * d$z + stats::rnorm(n)
  d$y <- 1 + 0.7 * d$xm - 0.3 * d$z + stats::rnorm(n)
  d$xm[sample.int(n, 12)] <- NA
  d
}

# The model of each case, built from whichever package `ns` names, so
# brms and frmtmb each get their own bf() and family objects. Returns
# list(formula, data, family), family NULL where the formula carries it.
mvprior_model <- function(case, ns) {
  g <- function(nm) get(nm, envir = asNamespace(ns))
  bf <- g("bf")
  rescor <- g("set_rescor")
  d <- mvprior_data()
  switch(case,
    mv2 = list(bf(y1 ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)) + rescor(FALSE),
               d, stats::gaussian()),
    mvsig = list(bf(y1 ~ x, sigma ~ z) + bf(y2 ~ x) + rescor(FALSE), d,
                 stats::gaussian()),
    mixed = list(bf(y1 ~ x, family = stats::gaussian()) +
                   bf(cnt ~ x, family = stats::poisson()) + rescor(FALSE),
                 d, NULL),
    rescor = list(bf(mvbind(y1, y2) ~ x + (1 + x | g)) + rescor(TRUE),
                  d, stats::gaussian()),
    mvstudent = list(bf(y1 ~ x) + bf(y2 ~ x) + rescor(FALSE), d,
                     g("student")()),
    mvar = list(bf(y1 ~ x + ar(time = t, gr = g, cov = TRUE)) + bf(y2 ~ x) +
                  rescor(FALSE), d, stats::gaussian()),
    nl = list(bf(ynl ~ a * exp(-b * ax), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
              d, stats::gaussian()),
    uni = list(bf(y ~ x + (1 | g), sigma ~ z), d, stats::gaussian()),
    uni2 = list(bf(y ~ x), d, stats::gaussian()),
    mvac = list(bf(y1 ~ x + ma(time = t, gr = g, cov = TRUE)) +
                  bf(y2 ~ x + cosy(time = t, gr = g)) + rescor(FALSE), d,
                stats::gaussian()),
    mvunstr = list(bf(y1 ~ x + unstr(time = t4, gr = g4)) + bf(y2 ~ x) +
                     rescor(FALSE), d, stats::gaussian()),
    uni3 = list(bf(y ~ 1), d, stats::gaussian()),
    cat = list(bf(cat ~ x), d, g("categorical")()),
    mix = list(bf(ym ~ x), d, g("mixture")(stats::gaussian(),
                                            stats::gaussian())),
    mvcat = list(bf(y1 ~ x, family = stats::gaussian()) +
                   bf(cat ~ x, family = g("categorical")()) + rescor(FALSE),
                 d, NULL),
    mvnl = list(bf(y1 ~ x) +
                  bf(ynl ~ a * exp(-b * ax), a ~ 1 + z, b ~ 1, nl = TRUE) +
                  rescor(FALSE), d, stats::gaussian()),
    catre = list(bf(cat ~ x + (1 | g)), d, g("categorical")()),
    mixre = list(bf(ym ~ x + (1 | g)), d,
                 g("mixture")(stats::gaussian(), stats::gaussian())),
    noslope = list(bf(y ~ 1 + (1 | g), sigma ~ 1), d, stats::gaussian()),
    mvnoslope = list(bf(y1 ~ x) + bf(y2 ~ 1) + rescor(FALSE), d,
                     stats::gaussian()),
    catnoslope = list(bf(cat ~ 1), d, g("categorical")()),
    uniar = list(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)), d,
                 stats::gaussian()),
    rescor2 = list(bf(mvbind(y1, y2) ~ x) + rescor(TRUE), d,
                   stats::gaussian()),
    nl2 = list(bf(ynl ~ a * exp(-b * ax), a ~ 1 + z, b ~ 1, nl = TRUE), d,
               stats::gaussian()),
    mi = list(bf(y ~ mi(xm) + z) + bf(xm | mi() ~ z) + rescor(FALSE),
              mvprior_data_mi(), stats::gaussian()))
}
