# Lane bigx: model shapes beyond the grid, each against the same model
# on the centered column (an exact reparametrization): an ordinal fit,
# a shifted column in sigma's design, an interaction with a factor, two
# shifted columns, a smooth beside a shifted column, a year column, a
# bound on the shifted coefficient (the optimizer keeps its natural
# coordinates), and raw polynomials of degree 6 with (1 | g), against
# lme4's maximum.
#   BIGX_ARM=before|after Rscript dev/bigx-extra.R
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
suppressPackageStartupMessages(library(lme4))
ll <- function(x) if (x$err) NA_real_ else as.numeric(logLik(x$f))
co <- function(x, nm) {
  if (x$err) return(NA_real_)
  b <- fixef(x$f)
  if (is.matrix(b)) b <- b[, "Estimate"]
  if (nm %in% names(b)) b[[nm]] else NA_real_
}
sec <- function(x, nm) {
  if (x$err) return(NA_real_)
  v <- suppressWarnings(tryCatch(vcov(x$f), error = function(e) NULL))
  if (is.null(v) || !nm %in% rownames(v)) return(NA_real_)
  sqrt(v[nm, nm])
}
say <- function(lbl, xw, xc, nw, nc) {
  cat(sprintf("%-34s dll %+.3e slope z %+9.4f SE ratio %.5f code %s warn %d %s\n",
              lbl, ll(xw) - ll(xc), (co(xw, nw) - co(xc, nc)) / sec(xc, nc),
              sec(xw, nw) / sec(xc, nc),
              if (xw$err) "ERR" else format(xw$f$opt$convergence),
              length(xw$w), if (length(xw$w)) substr(xw$w[1], 1, 50) else ""))
}
for (s in 1:3) {
  cat("== data seed", s, "\n")
  # ordinal
  for (m in c(1e4, 1e5, 1e6)) {
    set.seed(s)
    G <- 15L
    g <- factor(rep(seq_len(G), each = 10))
    z <- rnorm(150)
    lat <- 0.7 * z + rnorm(G, 0, 0.7)[g] + rlogis(150)
    y <- ordered(as.integer(cut(lat, c(-Inf, -1, 0.5, Inf))))
    d <- data.frame(y, w = 5 * z + m, wc = 5 * z, g)
    say(sprintf("cumulative (1|g) 5z+%g", m),
        bigx_fit(bf(y ~ w + (1 | g)), family = cumulative(), data = d),
        bigx_fit(bf(y ~ wc + (1 | g)), family = cumulative(), data = d),
        "w", "wc")
  }
  # sigma's design
  for (m in c(1e5, 1e6)) {
    set.seed(s)
    z <- rnorm(200)
    d <- data.frame(y = 1 + 0.5 * z + rnorm(200) * exp(-0.5 + 0.3 * z),
                    w = z + m, z = z)
    say(sprintf("gaussian sigma ~ w, +%g", m),
        bigx_fit(bf(y ~ w, sigma ~ w), family = gaussian(), data = d),
        bigx_fit(bf(y ~ z, sigma ~ z), family = gaussian(), data = d),
        "w", "z")
  }
  # interaction with a factor
  for (m in c(1e6, 1e7)) {
    set.seed(s)
    G <- 15L
    g <- factor(rep(seq_len(G), length.out = 200))
    z <- rnorm(200)
    f <- factor(sample(c("a", "b"), 200, TRUE))
    y <- rpois(200, exp(0.3 + 0.4 * z - 0.3 * z * (f == "b") +
                          rnorm(G, 0, 0.5)[g]))
    d <- data.frame(y, w = z + m, z, f, g)
    say(sprintf("poisson w*f (1|g), +%g", m),
        bigx_fit(bf(y ~ w * f + (1 | g)), family = poisson(), data = d),
        bigx_fit(bf(y ~ z * f + (1 | g)), family = poisson(), data = d),
        "w:fb", "z:fb")
  }
  # two shifted columns
  set.seed(s)
  z1 <- rnorm(200)
  z2 <- rnorm(200)
  d <- data.frame(y = rpois(200, exp(0.5 + 0.3 * z1 - 0.2 * z2)),
                  w1 = z1 + 1e6, w2 = 3 * z2 - 5e5, z1, z2)
  say("poisson w1 + w2 (1e6, -5e5)",
      bigx_fit(bf(y ~ w1 + w2), family = poisson(), data = d),
      bigx_fit(bf(y ~ z1 + z2), family = poisson(), data = d), "w2", "z2")
  # a smooth beside a shifted column
  set.seed(s)
  t <- runif(200)
  z <- rnorm(200)
  d <- data.frame(y = sin(2 * pi * t) + 0.4 * z + rnorm(200, 0, 0.3),
                  t, w = z + 1e6, z)
  say("gaussian s(t) + (x + 1e6)",
      bigx_fit(bf(y ~ s(t) + w), family = gaussian(), data = d),
      bigx_fit(bf(y ~ s(t) + z), family = gaussian(), data = d), "w", "z")
  # a year column (dev/reviews/2026-11-02-impsamp.md)
  set.seed(s)
  G <- 15L
  g <- factor(rep(seq_len(G), each = 10))
  zz <- rnorm(150)
  year <- round(2000 + 5 * zz)
  d <- data.frame(year, yc = year - 2000, g,
                  y = rpois(150, exp(0.3 + 0.08 * (year - 2000) +
                                       rnorm(G, 0, 0.5)[g])),
                  yb = rbinom(150, 1, plogis(0.08 * (year - 2000) +
                                               rnorm(G, 0, 0.5)[g])))
  say("poisson year (1|g)",
      bigx_fit(bf(y ~ year + (1 | g)), family = poisson(), data = d),
      bigx_fit(bf(y ~ yc + (1 | g)), family = poisson(), data = d),
      "year", "yc")
  say("bernoulli year (1|g)",
      bigx_fit(bf(yb ~ year + (1 | g)), family = bernoulli(), data = d),
      bigx_fit(bf(yb ~ yc + (1 | g)), family = bernoulli(), data = d),
      "year", "yc")
  # a bound on the shifted coefficient: natural coordinates
  d <- bigx_data("poisson", 1e5, s)
  pr <- set_prior("normal(0, 10)", class = "b", lb = 0)
  say("poisson (1|g) +1e5, b lb = 0",
      bigx_fit(bf(y ~ w + (1 | g)), family = poisson(), data = d,
               prior = pr),
      bigx_fit(bf(y ~ z + (1 | g)), family = poisson(), data = d,
               prior = pr), "w", "z")
  # raw polynomial of degree 6 on [1, 2] with (1 | g), against lme4
  set.seed(s)
  G <- 20L
  g <- factor(rep(seq_len(G), each = 10))
  x <- runif(200, 1, 2)
  d <- data.frame(x, g, y = sin(3 * x) + rnorm(G, 0, 0.5)[g] +
                    rnorm(200, 0, 0.2))
  fw <- bigx_fit(bf(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) +
                      (1 | g)), family = gaussian(), data = d)
  lm4 <- lmer(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + (1 | g),
              data = d, REML = FALSE)
  cat(sprintf("%-34s dll vs lmer %+.3e code %s warn %d %s\n",
              "gaussian raw poly 6 (1|g)", ll(fw) - as.numeric(logLik(lm4)),
              format(fw$f$opt$convergence), length(fw$w),
              if (length(fw$w)) substr(fw$w[1], 1, 50) else ""))
}
