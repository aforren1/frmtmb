# Designs NOT on the lane's grid. Each returns data, the model call, and a
# hand-written group log-likelihood `ll(par, d, u)`: u is a numeric vector of
# node values (scalar RE) or a 2-col matrix (2-d), and the result is the sum
# over the group's rows at each node. `sdpar` maps par to the RE sd or Sigma.
# The densities are base R's d*() functions, never the package's.
rows_ll <- function(lp_fun, dens) {
  function(par, d, u) {
    if (is.matrix(u)) m <- nrow(u) else m <- length(u)
    out <- numeric(m)
    for (i in seq_len(m)) {
      ui <- if (is.matrix(u)) u[i, ] else u[i]
      out[i] <- sum(dens(par, d, lp_fun(par, d, ui)))
    }
    out
  }
}
d_common <- function(sd1 = 0.8) NULL

design_list <- list(
  # unbalanced groups, 1 to 40 rows, Bernoulli
  unbal_bern = function(seed) {
    set.seed(seed)
    ns <- c(1, 1, 1, 2, 2, 3, 5, 8, 13, 21, 34, 40, 1, 2, 3, 1, 2, 5, 8, 2)
    g <- factor(rep(seq_along(ns), ns))
    x <- rnorm(length(g))
    u <- rnorm(length(ns), 0, 1.2)
    y <- rbinom(length(g), 1, plogis(-0.3 + 0.7 * x + u[g]))
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = bernoulli(), draws = 500L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dbinom(d$y, 1, plogis(par[1] + par[2] * d$x + ui), log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  },
  # large RE variance
  bigsd_bern = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:30, each = 5))
    x <- rnorm(150)
    u <- rnorm(30, 0, 4)
    y <- rbinom(150, 1, plogis(0.5 + 0.5 * x + u[g]))
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = bernoulli(), draws = 500L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dbinom(d$y, 1, plogis(par[1] + par[2] * d$x + ui), log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  },
  # near separation: most groups are all 0 or all 1 and x is strong
  nearsep = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:25, each = 6))
    x <- rep(c(-1, -1, -1, 1, 1, 1), 25) + rnorm(150, 0, 0.1)
    u <- rnorm(25, 0, 2)
    y <- rbinom(150, 1, plogis(4 * x + u[g]))
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = bernoulli(), draws = 500L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dbinom(d$y, 1, plogis(par[1] + par[2] * d$x + ui), log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  },
  # small counts
  smallpois = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:30, each = 3))
    x <- rnorm(90)
    u <- rnorm(30, 0, 1)
    y <- rpois(90, exp(-1.5 + 0.3 * x + u[g]))
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = poisson(), draws = 500L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dpois(d$y, exp(par[1] + par[2] * d$x + ui), log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  },
  # offset with a wide range
  offpois = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 6))
    x <- rnorm(90)
    ex <- exp(rnorm(90, 0, 1.5))
    u <- rnorm(15, 0, 0.7)
    y <- rpois(90, ex * exp(0.2 + 0.4 * x + u[g]))
    list(data = data.frame(y, x, g, ex), formula = bf(y ~ x + offset(log(ex)) + (1 | g)),
         family = poisson(), draws = 200L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dpois(d$y, d$ex * exp(par[1] + par[2] * d$x + ui), log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  },
  # case weights, integer and not, on a binomial
  wtbin = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 6))
    x <- rnorm(90)
    w <- sample(c(0.5, 1, 3), 90, TRUE)
    u <- rnorm(15, 0, 1)
    y <- rbinom(90, 1, plogis(0.2 + 0.8 * x + u[g]))
    list(data = data.frame(y, x, g, w), formula = bf(y | weights(w) ~ x + (1 | g)),
         family = bernoulli(), draws = 300L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(d$w * dbinom(d$y, 1, plogis(par[1] + par[2] * d$x + ui),
                              log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  },
  # non-canonical: probit and cloglog
  probit = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:20, each = 5))
    x <- rnorm(100)
    u <- rnorm(20, 0, 1)
    y <- rbinom(100, 1, pnorm(0.2 + 0.5 * x + u[g]))
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = bernoulli(link = "probit"), draws = 300L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dbinom(d$y, 1, pnorm(par[1] + par[2] * d$x + ui), log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  },
  cloglog = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:20, each = 5))
    x <- rnorm(100)
    u <- rnorm(20, 0, 1)
    y <- rbinom(100, 1, 1 - exp(-exp(-0.5 + 0.5 * x + u[g])))
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = bernoulli(link = "cloglog"), draws = 300L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             pr <- 1 - exp(-exp(par[1] + par[2] * d$x + ui))
             sum(dbinom(d$y, 1, pr, log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  },
  # gaussian mean with a log link: non-canonical and non-linear in u
  gaussloglink = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 5))
    x <- rnorm(75)
    u <- rnorm(15, 0, 0.6)
    y <- rnorm(75, exp(0.5 + 0.3 * x + u[g]), 0.5)
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = gaussian(link = "log"), draws = 300L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dnorm(d$y, exp(par[1] + par[2] * d$x + ui), exp(par[3]),
                       log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[4]))
  },
  # multi-dpar: sigma has its own regression, RE in mu only
  distreg = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 6))
    x <- rnorm(90)
    z <- rnorm(90)
    u <- rnorm(15, 0, 0.7)
    y <- rgamma(90, shape = exp(0.5 + 0.5 * z),
                rate = exp(0.5 + 0.5 * z) / exp(0.3 + 0.4 * x + u[g]))
    list(data = data.frame(y, x, z, g),
         formula = bf(y ~ x + (1 | g), shape ~ z),
         family = Gamma(link = "log"), draws = 300L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sh <- exp(par[3] + par[4] * d$z)
             mu <- exp(par[1] + par[2] * d$x + ui)
             sum(dgamma(d$y, shape = sh, rate = sh / mu, log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[5]))
  },
  # RE only in the dispersion
  re_sigma = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:20, each = 8))
    x <- rnorm(160)
    u <- rnorm(20, 0, 0.5)
    y <- rnorm(160, 0.2 + 0.5 * x, exp(-0.2 + u[g]))
    list(data = data.frame(y, x, g),
         formula = bf(y ~ x, sigma ~ 1 + (1 | g)),
         family = gaussian(), draws = 300L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dnorm(d$y, par[1] + par[2] * d$x, exp(par[3] + ui), log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[4]))
  },
  # small groups in the dispersion: the mode is near a saddle
  re_sigma_small = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:30, each = 2))
    x <- rnorm(60)
    u <- rnorm(30, 0, 1)
    y <- rnorm(60, 0.2 + 0.5 * x, exp(-0.2 + u[g]))
    list(data = data.frame(y, x, g),
         formula = bf(y ~ x, sigma ~ 1 + (1 | g)),
         family = gaussian(), draws = 500L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dnorm(d$y, par[1] + par[2] * d$x, exp(par[3] + ui), log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[4]))
  },
  # constant dpar: nu is mapped; a student RE fit, non-log-concave in u
  t_nu5 = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 6))
    x <- rnorm(90)
    u <- rnorm(15, 0, 0.7)
    y <- 0.3 + 0.5 * x + u[g] + 0.5 * rt(90, 5)
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g), nu = 5),
         family = student(), draws = 300L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             s <- exp(par[3])
             sum(dt((d$y - par[1] - par[2] * d$x - ui) / s, 5, log = TRUE) -
                   log(s))
           }, 0)
         },
         sdof = function(par) exp(par[4]))
  },
  # student with heavy outliers in small groups: multimodal integrands, the
  # PSD cut and the rate's denominator are both exposed
  t_outlier = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 3))
    x <- rnorm(45)
    u <- rnorm(15, 0, 0.7)
    y <- 0.3 + 0.5 * x + u[g] + 0.3 * rt(45, 2)
    y[seq(1, 45, by = 9)] <- y[seq(1, 45, by = 9)] + 6
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g), nu = 3),
         family = student(), draws = 500L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             s <- exp(par[3])
             sum(dt((d$y - par[1] - par[2] * d$x - ui) / s, 3, log = TRUE) -
                   log(s))
           }, 0)
         },
         sdof = function(par) exp(par[4]))
  },
  # fixed sigma (a mapped constant), gaussian: exact at any draw count
  gauss_sigma_fixed = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 4))
    x <- rnorm(60)
    u <- rnorm(15, 0, 0.7)
    y <- 0.3 + 0.5 * x + u[g] + rnorm(60, 0, 0.5)
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g), sigma = 0.5),
         family = gaussian(), draws = 64L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dnorm(d$y, par[1] + par[2] * d$x + ui, 0.5, log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  },
  # gaussian with se(): residual sd is known per row and unequal
  gauss_se = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:12, each = 4))
    x <- rnorm(48)
    s <- runif(48, 0.1, 1.5)
    u <- rnorm(12, 0, 0.7)
    y <- 0.3 + 0.5 * x + u[g] + rnorm(48, 0, s)
    list(data = data.frame(y, x, g, s), formula = bf(y | se(s) ~ x + (1 | g)),
         family = gaussian(), draws = 64L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dnorm(d$y, par[1] + par[2] * d$x + ui, d$s, log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[3]))
  },
  # negative binomial, shape estimated, overdispersed small groups
  nb = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 4))
    x <- rnorm(60)
    u <- rnorm(15, 0, 1)
    y <- rnbinom(60, size = 1, mu = exp(1 + 0.4 * x + u[g]))
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = negbinomial(), draws = 300L,
         ll = function(par, d, u) {
           vapply(u, function(ui) {
             sum(dnbinom(d$y, size = exp(par[3]),
                         mu = exp(par[1] + par[2] * d$x + ui), log = TRUE))
           }, 0)
         },
         sdof = function(par) exp(par[4]))
  }
)
