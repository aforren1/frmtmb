# Reviewer 2 cases: models and specifications, independent of the
# worker's dev/mvprior-cases.R. Seed 4417, n = 150, 4-level categorical.

r2_data <- function() {
  set.seed(4417)
  n <- 150
  d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                  g = factor(rep(1:10, length.out = n)),
                  h = factor(rep(1:5, each = 30)))
  u <- rnorm(10)[d$g]
  e1 <- 0.8 * d$x + u + rlogis(n)
  e2 <- -0.5 * d$x + d$z + rlogis(n)
  e3 <- 0.3 * d$w + rlogis(n)
  lvl <- apply(cbind(0, e1, e2, e3), 1, which.max)
  d$cat4 <- factor(c("a", "b", "c", "d")[lvl])
  d$cat2 <- factor(ifelse(d$x + rlogis(n) > 0, "yes", "no"))
  comp <- sample(1:3, n, TRUE)
  d$ym3 <- c(-4, 0, 5)[comp] + 0.7 * d$x + rnorm(n, 0, 0.6)
  d$ym <- ifelse(comp == 1, -2 + d$x, 3 + d$x) + rnorm(n, 0, 0.7)
  d$yg <- 2 + 0.5 * d$x + u + rnorm(n)
  d$ord <- factor(cut(d$x + 0.5 * d$z + rlogis(n), c(-Inf, -1, 0, 1.2, Inf),
                      labels = FALSE), ordered = TRUE)
  d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x)) * rbinom(n, 1, 0.7)
  d$pos <- ifelse(rbinom(n, 1, 0.3) == 1, 0, exp(0.5 + 0.3 * d$x +
                                                    rnorm(n, 0, 0.5)))
  d$mcnt <- rpois(n, exp(1 + 0.3 * d$x))
  d
}

sp <- function(class = "b", coef = "", group = "", resp = "", dpar = "",
               nlpar = "", prior = "normal(0, 5)", lb = NA, ub = NA) {
  list(prior = prior, class = class, coef = coef, group = group,
       resp = resp, dpar = dpar, nlpar = nlpar, lb = lb, ub = ub)
}

# name -> list(what, build(ns) -> list(formula, family), specs)
r2_cases <- list()
r2_add <- function(name, what, build, specs) {
  r2_cases[[name]] <<- list(what = what, build = build, specs = specs)
}
G <- function(ns, nm) get(nm, envir = asNamespace(ns))

r2_add("cat4", "bf(cat4 ~ x + z), categorical, 4 levels",
  function(ns) list(G(ns, "bf")(cat4 ~ x + z), G(ns, "categorical")()),
  list(sp("b"), sp("Intercept"), sp("b", coef = "x"),
       sp("b", dpar = "mub"), sp("b", dpar = "mud", coef = "z"),
       sp("Intercept", dpar = "muc"), sp("b", dpar = "mue"),
       sp("b", dpar = "mub", prior = "", lb = 0),
       sp("b", dpar = "muc", coef = "x", prior = "", ub = 1),
       sp("b", dpar = "mub", prior = "normal(0, 1)", lb = -2, ub = 2),
       sp("Intercept", coef = "", dpar = "mud", prior = "student_t(3, 0, 1)"),
       sp("b", coef = "x", dpar = "mu")))

r2_add("cat2", "bf(cat2 ~ x), categorical, 2 levels (one beyond ref)",
  function(ns) list(G(ns, "bf")(cat2 ~ x), G(ns, "categorical")()),
  list(sp("b"), sp("Intercept"), sp("b", dpar = "muyes"),
       sp("Intercept", dpar = "muyes"), sp("b", coef = "x", dpar = "muyes"),
       sp("b", dpar = "mu")))

r2_add("cat4re", "bf(cat4 ~ x + (1 | ID | g)), categorical",
  function(ns) list(G(ns, "bf")(cat4 ~ x + (1 | ID | g)),
                    G(ns, "categorical")()),
  list(sp("sd"), sp("sd", group = "g"), sp("sd", dpar = "mub"),
       sp("sd", dpar = "muc", group = "g"),
       sp("sd", dpar = "mud", group = "g", coef = "Intercept"),
       sp("cor", prior = "lkj(2)"), sp("cor", group = "g", prior = "lkj(2)"),
       sp("cor", dpar = "mub", prior = "lkj(2)"),
       sp("b", dpar = "mub"), sp("sd", group = "h")))

r2_add("cat4re2", "bf(cat4 ~ x + (1 | g)), categorical, blocks per dpar",
  function(ns) list(G(ns, "bf")(cat4 ~ x + (1 | g)), G(ns, "categorical")()),
  list(sp("sd"), sp("sd", dpar = "mub"), sp("sd", dpar = "muc", group = "g"),
       sp("cor", prior = "lkj(2)")))

r2_add("catpart", "bf(cat4 ~ x, mud ~ z), dpar formula for one category",
  function(ns) list(G(ns, "bf")(cat4 ~ x, mud ~ z), G(ns, "categorical")()),
  list(sp("b"), sp("b", dpar = "mub"), sp("b", dpar = "mud"),
       sp("b", dpar = "mud", coef = "z"), sp("b", dpar = "mud", coef = "x"),
       sp("b", dpar = "muc", coef = "x"), sp("Intercept", dpar = "mud"),
       sp("Intercept")))

r2_add("mult", "bf(mcnt | trials(tt) ~ x) placeholder",
  NULL, list())
r2_cases$mult <- NULL

r2_add("mix3", "bf(ym3 ~ x), mixture(gaussian x3)",
  function(ns) list(G(ns, "bf")(ym3 ~ x),
                    G(ns, "mixture")(stats::gaussian(), stats::gaussian(),
                                     stats::gaussian())),
  list(sp("b"), sp("Intercept"), sp("b", dpar = "mu3"),
       sp("Intercept", dpar = "mu2"), sp("b", dpar = "mu3", coef = "x"),
       sp("sigma", prior = "student_t(3, 0, 2.5)"),
       sp("sigma3", prior = "student_t(3, 0, 2.5)"),
       sp("sigma2", prior = "exponential(1)"),
       sp("theta1", prior = "normal(0, 1)"),
       sp("Intercept", dpar = "theta2"),
       sp("b", dpar = "mu4"),
       sp("b", dpar = "mu1", prior = "", lb = 0)))

r2_add("mixdiff", "bf(ym ~ x), mixture(gaussian, student)",
  function(ns) list(G(ns, "bf")(ym ~ x),
                    G(ns, "mixture")(stats::gaussian(), G(ns, "student")())),
  list(sp("b"), sp("b", dpar = "mu1"), sp("b", dpar = "mu2"),
       sp("Intercept", dpar = "mu2"), sp("nu2", prior = "gamma(2, 0.1)"),
       sp("nu", prior = "gamma(2, 0.1)"),
       sp("sigma1", prior = "student_t(3, 0, 2.5)")))

r2_add("mixsig", "bf(ym ~ x, sigma2 ~ z, theta2 ~ w), mixture(gaussian x2)",
  function(ns) list(G(ns, "bf")(ym ~ x, sigma2 ~ z, theta2 ~ w),
                    G(ns, "mixture")(stats::gaussian(), stats::gaussian())),
  list(sp("b"), sp("b", dpar = "sigma2"), sp("Intercept", dpar = "sigma2"),
       sp("b", dpar = "theta2"), sp("b", dpar = "theta2", coef = "w"),
       sp("sigma1", prior = "student_t(3, 0, 2.5)"),
       sp("Intercept", dpar = "mu1"), sp("b", dpar = "sigma1")))

r2_add("mixre", "bf(ym ~ x + (1 | g)), mixture(gaussian x2)",
  function(ns) list(G(ns, "bf")(ym ~ x + (1 | g)),
                    G(ns, "mixture")(stats::gaussian(), stats::gaussian())),
  list(sp("sd"), sp("sd", dpar = "mu1"), sp("sd", dpar = "mu2", group = "g"),
       sp("sd", group = "g")))

r2_add("mvcat", "bf(yg ~ x) + bf(cat4 ~ x + (1 | g), categorical)",
  function(ns) list(G(ns, "bf")(yg ~ x, family = stats::gaussian()) +
                      G(ns, "bf")(cat4 ~ x + (1 | g),
                                  family = G(ns, "categorical")()) +
                      G(ns, "set_rescor")(FALSE), NULL),
  list(sp("b", resp = "cat4"), sp("b", resp = "yg"),
       sp("b", dpar = "mub", resp = "cat4"),
       sp("b", dpar = "mub"), sp("Intercept", dpar = "mud", resp = "cat4"),
       sp("sd", resp = "cat4"), sp("sd", resp = "cat4", dpar = "muc"),
       sp("sd", dpar = "muc"), sp("sigma", resp = "yg",
                                  prior = "student_t(3, 0, 2.5)"),
       sp("b", resp = "cat4", dpar = "mub", coef = "x"),
       sp("b", resp = "yg", dpar = "mub")))

r2_add("cumul", "bf(ord ~ x + z + (1 | g)), cumulative",
  function(ns) list(G(ns, "bf")(ord ~ x + z + (1 | g)), G(ns, "cumulative")()),
  list(sp("b"), sp("b", coef = "x"), sp("Intercept"),
       sp("Intercept", coef = "2"), sp("sd"), sp("b", dpar = "mu")))

r2_add("sratio", "bf(ord ~ x), sratio", function(ns)
  list(G(ns, "bf")(ord ~ x), G(ns, "sratio")()),
  list(sp("b"), sp("Intercept")))
r2_add("acat", "bf(ord ~ x), acat", function(ns)
  list(G(ns, "bf")(ord ~ x), G(ns, "acat")()),
  list(sp("b"), sp("Intercept")))
r2_add("cratio", "bf(ord ~ x, disc ~ z), cratio", function(ns)
  list(G(ns, "bf")(ord ~ x, disc ~ z), G(ns, "cratio")()),
  list(sp("b"), sp("Intercept"), sp("b", dpar = "disc")))

r2_add("zip", "bf(cnt ~ x + (1 | g), zi ~ z), zero_inflated_poisson",
  function(ns) list(G(ns, "bf")(cnt ~ x + (1 | g), zi ~ z),
                    G(ns, "zero_inflated_poisson")()),
  list(sp("b"), sp("b", coef = "x"), sp("Intercept"), sp("sd"),
       sp("b", dpar = "zi"), sp("Intercept", dpar = "zi")))

r2_add("hln", "bf(pos ~ x), hurdle_lognormal",
  function(ns) list(G(ns, "bf")(pos ~ x), G(ns, "hurdle_lognormal")()),
  list(sp("b"), sp("Intercept"), sp("hu", prior = "beta(1, 1)"),
       sp("sigma", prior = "student_t(3, 0, 2.5)")))

# item 3: refusals must not over-reach
r2_add("slopes_in_sigma", "bf(y ~ 1, sigma ~ x): slopes only in sigma",
  function(ns) list(G(ns, "bf")(yg ~ 1, sigma ~ x), stats::gaussian()),
  list(sp("b"), sp("b", dpar = "sigma"), sp("b", dpar = "sigma", coef = "x"),
       sp("Intercept"), sp("b", coef = "x")))

r2_add("slopes_in_nl", "bf(yg ~ a + b * x, a ~ 1, b ~ 1 + z, nl = TRUE)",
  function(ns) list(G(ns, "bf")(yg ~ a + b * x, a ~ 1, b ~ 1 + z, nl = TRUE),
                    stats::gaussian()),
  list(sp("b"), sp("b", nlpar = "b"), sp("b", nlpar = "a"),
       sp("b", nlpar = "b", coef = "z"), sp("b", nlpar = "a", coef = "Intercept")))

r2_add("icpt_only", "bf(yg ~ 1): coef = on intercept-only",
  function(ns) list(G(ns, "bf")(yg ~ 1), stats::gaussian()),
  list(sp("b", coef = "Intercept"), sp("b", coef = "x"), sp("b"),
       sp("Intercept")))

r2_add("mvdrop", "bf(yg ~ x) + bf(ym ~ 1 + (1 | g)): mv, second no slope",
  function(ns) list(G(ns, "bf")(yg ~ x) + G(ns, "bf")(ym ~ 1 + (1 | g)) +
                      G(ns, "set_rescor")(FALSE), stats::gaussian()),
  list(sp("b", resp = "yg"), sp("b", resp = "ym"),
       sp("Intercept", resp = "ym"), sp("sd", resp = "ym"),
       sp("sd", resp = "yg"), sp("b", resp = "nope"),
       sp("sigma", resp = "ym", prior = "student_t(3, 0, 2.5)"),
       sp("b", coef = "x", resp = "yg")))
