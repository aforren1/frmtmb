# Review of lane fitedge, item 1: fresh designs for gp() in a nonlinear body.
# Sourced by review-fitedge-gpfresh.R and review-fitedge-gpprobe.R.
# Seeds: set.seed(52000 + 100 * seed + n + design index).
rv_curve <- function(x, ci) {
  switch(ci, sin(x), 2 * sin(3 * x) + x, exp(0.3 * x))
}
rv_designs <- c("E", "F", "G", "H", "I", "J", "K", "L")
rv_mk <- function(design, seed, n, ci) {
  di <- match(design, rv_designs)
  set.seed(52000 + 100 * seed + n + 10 * di + ci)
  x <- sort(round(runif(n, 0, 10), 3))
  z <- round(rnorm(n), 3)
  g <- factor(sample(1:8, n, TRUE))
  fby <- factor(sample(c("p", "q"), n, TRUE))
  f <- rv_curve(x, ci)
  fs <- f / max(abs(f))
  ge <- rnorm(8, 0, 0.5)[g]
  y <- switch(design,
    # a * exp(b z): a carries the field
    E = (2 + fs) * exp(0.4 * z) + 0.1 * rnorm(n),
    # logistic in z, the field in the asymptote
    F = (3 + fs) / (1 + exp(-1.5 * (z - 0.2))) + 0.1 * rnorm(n),
    # a + b z^2
    G = f + 0.5 * z^2 + 0.1 * rnorm(n),
    # field and a group effect in the same nlpar
    H = f + ge + 0.3 * z + 0.1 * rnorm(n),
    # gp(x, by = fby)
    I = ifelse(fby == "p", f, 0.5 * f) + 0.3 * z + 0.1 * rnorm(n),
    # bernoulli
    J = rbinom(n, 1, plogis(1.5 * fs + 0.5 * z)),
    # Gamma(log)
    K = rgamma(n, 5, 5 / exp(0.5 + 0.5 * fs + 0.3 * z)),
    # the field in the exponent's coefficient: J_b = a z exp(b z)
    L = 2 * exp((0.3 + 0.2 * fs) * z) + 0.1 * rnorm(n))
  data.frame(x, z, g, fby, y)
}
rv_spell <- list(
  E = list(nl = bf(y ~ a * exp(b * z), a ~ gp(x), b ~ 1, nl = TRUE),
           lin = NULL, fam = gaussian()),
  F = list(nl = bf(y ~ a / (1 + exp(-b * (z - c))), a ~ gp(x), b ~ 1, c ~ 1,
                   nl = TRUE), lin = NULL, fam = gaussian()),
  G = list(nl = bf(y ~ a + b * z^2, a ~ gp(x), b ~ 1, nl = TRUE),
           lin = y ~ I(z^2) + gp(x), fam = gaussian()),
  H = list(nl = bf(y ~ a + b * z, a ~ gp(x) + (1 | g), b ~ 1, nl = TRUE),
           lin = y ~ z + gp(x) + (1 | g), fam = gaussian()),
  I = list(nl = bf(y ~ a + b * z, a ~ gp(x, by = fby), b ~ 1, nl = TRUE),
           lin = y ~ z + gp(x, by = fby), fam = gaussian()),
  J = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
           lin = y ~ z + gp(x), fam = bernoulli()),
  K = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
           lin = y ~ z + gp(x), fam = Gamma(link = "log")),
  L = list(nl = bf(y ~ a * exp(b * z), a ~ 1, b ~ gp(x), nl = TRUE),
           lin = NULL, fam = gaussian()))
