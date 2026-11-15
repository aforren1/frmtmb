# Lane gpfield: the data and spellings of the rows this lane probes, by tag.
# sweep:<body>:<fam>:<n>:<seed>[:k]  lane gpstop's review sweep, mk() of
#   dev/gpstop-fx2-sweep.R (set.seed(81100 + 1000 fam + 100 seed + n + 10
#   body)).
# l96:<design>:<curve>:<n>:<seed>  lane fitedge's lane96 set, mk() of
#   dev/gpstop-gpfresh.R (set.seed(18000 + 100 seed + n + curve)).
# ng:<case>:<curve>:<n>:<seed>  lane arstep's 150-row set, dev/gpstop-ng.R
#   (set.seed(74000 + 31 seed + n)).
# gate:<approx_a|approx_b>  the review's gate rows (dev/review-fitedge-
#   designs.R rv_mk("E" or "L", 2, 120, 2)).
gf_bodies <- c("exp_a", "exp_b", "logis_a", "mm_V", "mm_K", "pow_a", "pow_b")
gf_sweep_spell <- function(body, k = NULL) {
  g <- if (is.null(k)) quote(gp(x)) else bquote(gp(x, k = .(k)))
  switch(body,
    exp_a = bf(y ~ a * exp(b * z), as.formula(bquote(a ~ .(g))), b ~ 1,
               nl = TRUE),
    exp_b = bf(y ~ a * exp(b * z), a ~ 1, as.formula(bquote(b ~ .(g))),
               nl = TRUE),
    logis_a = bf(y ~ a / (1 + exp(-b * z)), as.formula(bquote(a ~ .(g))),
                 b ~ 1, nl = TRUE),
    mm_V = bf(y ~ V * w / (K + w), as.formula(bquote(V ~ .(g))), K ~ 1,
              nl = TRUE),
    mm_K = bf(y ~ V * w / (K + w), V ~ 1, as.formula(bquote(K ~ .(g))),
              nl = TRUE),
    pow_a = bf(y ~ a * w^b, as.formula(bquote(a ~ .(g))), b ~ 1, nl = TRUE),
    pow_b = bf(y ~ a * w^b, a ~ 1, as.formula(bquote(b ~ .(g))), nl = TRUE))
}
gf_sweep_mk <- function(body, fam, n, seed) {
  set.seed(81100 + 1000 * match(fam, c("gaussian", "poisson", "bernoulli",
                                       "student")) +
             100 * seed + n + 10 * match(body, gf_bodies))
  x <- sort(round(runif(n, 0, 10), 3))
  z <- round(rnorm(n), 3)
  w <- round(runif(n, 0.2, 5), 3)
  f <- if (seed %% 2) sin(x) else 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  eta <- switch(body,
    exp_a = (2 + fs) * exp(0.4 * z),
    exp_b = 2 * exp((0.3 + 0.2 * fs) * z),
    logis_a = (3 + fs) / (1 + exp(-1.5 * z)),
    mm_V = (3 + fs) * w / (1 + w),
    mm_K = 3 * w / (1 + 0.6 * fs + w),
    pow_a = (2 + fs) * w^0.5,
    pow_b = w^(0.5 + 0.3 * fs))
  y <- switch(fam,
    gaussian = eta + 0.1 * rnorm(n),
    student = eta + 0.1 * rt(n, 3),
    poisson = rpois(n, exp(0.6 * eta)),
    bernoulli = rbinom(n, 1, plogis(1.5 * (eta - median(eta)) / sd(eta))))
  data.frame(x, z, w, y)
}
gf_famf <- function(fam) {
  switch(fam, gaussian = gaussian(), poisson = poisson(),
         student = student(), bernoulli = bernoulli())
}
gf_l96_curves <- c("sin", "wig", "exp", "step")
gf_l96_mk <- function(seed, n, ci, body) {
  set.seed(18000 + 100 * seed + n + ci)
  x <- sort(runif(n, 0, 10)); z <- rnorm(n); w <- runif(n, 1, 2)
  f <- switch(gf_l96_curves[ci], sin = sin(x), wig = 2 * sin(3 * x) + x,
              exp = exp(0.3 * x), step = sign(sin(1.5 * x)))
  y <- switch(body,
    A = f + 0.3 * z + 0.1 * rnorm(n),
    B = rpois(n, exp(1 + 0.5 * f / max(abs(f)) + 0.3 * z)),
    C = exp(0.5 + 0.5 * f / max(abs(f))) + 0.3 * z + 0.1 * rnorm(n),
    D = 1 + f * w + 0.1 * rnorm(n))
  data.frame(x, z, w, y)
}
gf_l96_spell <- list(
  A = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
           fam = gaussian()),
  B = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
           fam = poisson()),
  C = list(nl = bf(y ~ exp(a) + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
           fam = gaussian()),
  D = list(nl = bf(y ~ b + a * w, a ~ gp(x), b ~ 1, nl = TRUE),
           fam = gaussian()))
gf_ng_curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x,
                     wig = function(x) sin(7 * x) + 0.2 * x,
                     two = function(x) sin(x) + 0.4 * sin(10 * x),
                     smooth = function(x) sin(x),
                     step = function(x) 2 * (x > 5) - 1)
gf_ng_mk <- function(case, cn, n, seed) {
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  eta <- gf_ng_curves[[cn]](x)
  switch(case,
    pois_low = list(d = data.frame(x = x, y = rpois(n, exp(eta / 2))),
                    fml = y ~ gp(x), fam = poisson()),
    pois_mid = list(d = data.frame(x = x, y = rpois(n, exp(1.5 + eta / 2))),
                    fml = y ~ gp(x), fam = poisson()),
    bern = list(d = data.frame(x = x, y = rbinom(n, 1,
                                                 plogis(eta - mean(eta)))),
                fml = y ~ gp(x), fam = bernoulli()),
    binom5 = list(d = data.frame(x = x, tr = 5,
                                 y = rbinom(n, 5, plogis(eta - mean(eta)))),
                  fml = y | trials(tr) ~ gp(x), fam = binomial()),
    negbin = list(d = data.frame(x = x, y = rnbinom(n, size = 20,
                                                    mu = exp(4 + eta))),
                  fml = y ~ gp(x), fam = negbinomial()),
    gamma = list(d = data.frame(x = x, y = rgamma(n, 400, 400 / exp(eta))),
                 fml = y ~ gp(x), fam = Gamma(link = "log")),
    lognorm = list(d = data.frame(x = x,
                                  y = exp(eta + rnorm(n, sd = 0.02))),
                   fml = y ~ gp(x), fam = lognormal()),
    student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)),
                   fml = y ~ gp(x), fam = student()),
    gsigma = list(d = data.frame(x = x, z = rnorm(n),
                                 y = eta + 0.02 * rnorm(n)),
                  fml = bf(y ~ gp(x), sigma ~ z), fam = gaussian()),
    gar = list(d = data.frame(x = x, t = seq_len(n),
                              y = eta + as.numeric(arima.sim(list(ar = 0.4),
                                                             n, sd = 0.03))),
               fml = y ~ gp(x) + ar(time = t), fam = gaussian()))
}
gf_rv_curve <- function(x, ci) {
  switch(ci, sin(x), 2 * sin(3 * x) + x, exp(0.3 * x))
}
gf_gate_mk <- function(design) {
  di <- match(design, c("E", "F", "G", "H", "I", "J", "K", "L"))
  n <- 120L; seed <- 2; ci <- 2
  set.seed(52000 + 100 * seed + n + 10 * di + ci)
  x <- sort(round(runif(n, 0, 10), 3))
  z <- round(rnorm(n), 3)
  g <- factor(sample(1:8, n, TRUE))
  fby <- factor(sample(c("p", "q"), n, TRUE))
  f <- gf_rv_curve(x, ci)
  fs <- f / max(abs(f))
  ge <- rnorm(8, 0, 0.5)[g]
  y <- switch(design,
    E = (2 + fs) * exp(0.4 * z) + 0.1 * rnorm(n),
    L = 2 * exp((0.3 + 0.2 * fs) * z) + 0.1 * rnorm(n))
  data.frame(x, z, g, fby, y)
}
# list(fml, data, family) for a tag
gf_row <- function(tag) {
  p <- strsplit(tag, ":")[[1]]
  switch(p[1],
    sweep = {
      k <- if (length(p) > 5) as.integer(p[6])
      list(fml = gf_sweep_spell(p[2], k),
           data = gf_sweep_mk(p[2], p[3], as.integer(p[4]), as.integer(p[5])),
           family = gf_famf(p[3]))
    },
    l96 = {
      sp <- gf_l96_spell[[p[2]]]
      list(fml = sp$nl,
           data = gf_l96_mk(as.integer(p[5]), as.integer(p[4]),
                            match(p[3], gf_l96_curves), p[2]),
           family = sp$fam)
    },
    ng = {
      m <- gf_ng_mk(p[2], p[3], as.integer(p[4]), as.integer(p[5]))
      list(fml = m$fml, data = m$d, family = m$fam)
    },
    gate = {
      k <- if (length(p) > 2) as.integer(p[3]) else 20L
      if (p[2] == "approx_a") {
        list(fml = bf(y ~ a * exp(b * z), as.formula(bquote(a ~ gp(x, k = .(k)))),
                      b ~ 1, nl = TRUE),
             data = gf_gate_mk("E"), family = gaussian())
      } else {
        list(fml = bf(y ~ a * exp(b * z), a ~ 1,
                      as.formula(bquote(b ~ gp(x, k = .(k)))), nl = TRUE),
             data = gf_gate_mk("L"), family = gaussian())
      }
    })
}
