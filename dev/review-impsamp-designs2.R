# More designs not on the lane's grid: non-log-concave families, ordinal,
# censoring, truncation, the bounded families, and scales that break the fixed
# step of imp_anchor_slopes().
lv <- function(f) function(par, d, u) {
  vapply(u, function(ui) sum(f(par, d, ui)), 0)
}
design_list <- list(
  # zero inflation: the likelihood is not log-concave in u at y = 0
  zip = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:20, each = 6))
    x <- rnorm(120)
    u <- rnorm(20, 0, 1)
    mu <- exp(1 + 0.3 * x + u[g])
    y <- ifelse(runif(120) < 0.3, 0L, rpois(120, mu))
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = zero_inflated_poisson(), draws = 500L,
         ll = lv(function(par, d, ui) {
           zi <- plogis(par[3])
           m <- exp(par[1] + par[2] * d$x + ui)
           ifelse(d$y == 0, log(zi + (1 - zi) * exp(-m)),
                  log(1 - zi) + dpois(d$y, m, log = TRUE))
         }),
         sdof = function(par) exp(par[4]))
  },
  # a hurdle: the RE's curvature is the truncated Poisson's
  hurdle = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:20, each = 6))
    x <- rnorm(120)
    u <- rnorm(20, 0, 0.8)
    mu <- exp(0.5 + 0.3 * x + u[g])
    y <- rpois(120, mu)
    z <- runif(120) < 0.3
    y[z] <- 0L
    y[!z & y == 0] <- 1L
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = hurdle_poisson(), draws = 500L,
         ll = lv(function(par, d, ui) {
           hu <- plogis(par[3])
           m <- exp(par[1] + par[2] * d$x + ui)
           ifelse(d$y == 0, log(hu),
                  log(1 - hu) + dpois(d$y, m, log = TRUE) -
                    log(1 - exp(-m)))
         }),
         sdof = function(par) exp(par[4]))
  },
  # beta
  betafam = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:20, each = 5))
    x <- rnorm(100)
    u <- rnorm(20, 0, 0.7)
    mu <- plogis(0.2 + 0.4 * x + u[g])
    phi <- 6
    y <- pmin(pmax(rbeta(100, mu * phi, (1 - mu) * phi), 1e-4), 1 - 1e-4)
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = Beta(), draws = 300L,
         ll = lv(function(par, d, ui) {
           m <- plogis(par[1] + par[2] * d$x + ui)
           p <- exp(par[3])
           dbeta(d$y, m * p, (1 - m) * p, log = TRUE)
         }),
         sdof = function(par) exp(par[4]))
  },
  # censored and truncated rows
  cens_gauss = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 5))
    x <- rnorm(75)
    u <- rnorm(15, 0, 0.8)
    y0 <- 0.5 + 0.5 * x + u[g] + rnorm(75, 0, 0.5)
    cen <- ifelse(y0 > 1, "right", "none")
    y <- pmin(y0, 1)
    list(data = data.frame(y, x, g, cen), formula = bf(y | cens(cen) ~ x + (1 | g)),
         family = gaussian(), draws = 400L,
         ll = lv(function(par, d, ui) {
           m <- par[1] + par[2] * d$x + ui
           s <- exp(par[3])
           ifelse(d$cen == "right",
                  pnorm(d$y, m, s, lower.tail = FALSE, log.p = TRUE),
                  dnorm(d$y, m, s, log = TRUE))
         }),
         sdof = function(par) exp(par[4]))
  },
  trunc_pois = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:20, each = 5))
    x <- rnorm(100)
    u <- rnorm(20, 0, 0.7)
    y <- integer(100)
    for (i in 1:100) {
      repeat {
        v <- rpois(1, exp(0.3 + 0.4 * x[i] + u[g[i]]))
        if (v >= 1) break
      }
      y[i] <- v
    }
    list(data = data.frame(y, x, g), formula = bf(y | trunc(lb = 1) ~ x + (1 | g)),
         family = poisson(), draws = 300L,
         ll = lv(function(par, d, ui) {
           m <- exp(par[1] + par[2] * d$x + ui)
           dpois(d$y, m, log = TRUE) - log(1 - exp(-m))
         }),
         sdof = function(par) exp(par[3]))
  },
  # an ordinal response: three thresholds, a logit latent
  ordinal = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:25, each = 5))
    x <- rnorm(125)
    u <- rnorm(25, 0, 1)
    lat <- 0.7 * x + u[g] + rlogis(125)
    y <- as.integer(cut(lat, c(-Inf, -1, 0.5, 1.5, Inf)))
    list(data = data.frame(y = ordered(y), yi = y, x, g),
         formula = bf(y ~ x + (1 | g)), family = cumulative(), draws = 400L,
         ordinal = TRUE,
         ll = NULL, sdof = NULL)
  },
  # x at a scale of 1e5: h = 1e-4 moves the predictor by 10 at the slope's
  # coordinate, so imp_anchor_slopes' central difference is far from local
  scale1e5 = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 10))
    x <- rnorm(150)
    u <- rnorm(15, 0, 0.8)
    y <- rbinom(150, 1, plogis(0.2 + 0.7 * x + u[g]))
    xb <- x * 1e5
    list(data = data.frame(y, xb, g), formula = bf(y ~ xb + (1 | g)),
         family = bernoulli(), draws = 300L,
         ll = lv(function(par, d, ui) {
           dbinom(d$y, 1, plogis(par[1] + par[2] * d$xb + ui), log = TRUE)
         }),
         sdof = function(par) exp(par[3]))
  },
  pois1e5 = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 10))
    x <- rnorm(150)
    u <- rnorm(15, 0, 0.6)
    y <- rpois(150, exp(0.5 + 0.5 * x + u[g]))
    xb <- x * 1e5
    list(data = data.frame(y, xb, g), formula = bf(y ~ xb + (1 | g)),
         family = poisson(), draws = 300L,
         ll = lv(function(par, d, ui) {
           dpois(d$y, exp(par[1] + par[2] * d$xb + ui), log = TRUE)
         }),
         sdof = function(par) exp(par[3]))
  },
  # uncentered x around 2000: the intercept is far from 1 and the slope small
  year = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:15, each = 10))
    x <- 2000 + rnorm(150, 0, 5)
    u <- rnorm(15, 0, 0.8)
    y <- rbinom(150, 1, plogis(0.2 + 0.1 * (x - 2000) + u[g]))
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = bernoulli(), draws = 300L,
         ll = lv(function(par, d, ui) {
           dbinom(d$y, 1, plogis(par[1] + par[2] * d$x + ui), log = TRUE)
         }),
         sdof = function(par) exp(par[3]))
  },
  # a slope RE at x scale 1e-3: theta's coordinate is large on the log scale
  # and H[2,2] is small against the intercept's
  tiny_x_slope = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:30, each = 2))
    x <- rnorm(60)
    u <- rnorm(30, 0, 1.5)
    y <- rbinom(60, 1, plogis(u[g]))
    xs <- x * 1e-3
    list(data = data.frame(y, xs, g), formula = bf(y ~ xs + (1 | g)),
         family = bernoulli(), draws = 500L,
         ll = lv(function(par, d, ui) {
           dbinom(d$y, 1, plogis(par[1] + par[2] * d$xs + ui), log = TRUE)
         }),
         sdof = function(par) exp(par[3]))
  },
  # one huge group among small ones
  huge_group = function(seed) {
    set.seed(seed)
    ns <- c(2000, rep(3, 15))
    g <- factor(rep(seq_along(ns), ns))
    x <- rnorm(length(g))
    u <- rnorm(length(ns), 0, 1)
    y <- rbinom(length(g), 1, plogis(0.3 + 0.5 * x + u[g]))
    list(data = data.frame(y, x, g), formula = bf(y ~ x + (1 | g)),
         family = bernoulli(), draws = 200L,
         ll = lv(function(par, d, ui) {
           dbinom(d$y, 1, plogis(par[1] + par[2] * d$x + ui), log = TRUE)
         }),
         sdof = function(par) exp(par[3]))
  },
  # binomial with varying trials and a zero-weight group
  zero_wt = function(seed) {
    set.seed(seed)
    g <- factor(rep(1:12, each = 6))
    x <- rnorm(72)
    u <- rnorm(12, 0, 1)
    w <- rep(1, 72)
    w[g == 3] <- 0
    y <- rbinom(72, 1, plogis(0.2 + 0.6 * x + u[g]))
    list(data = data.frame(y, x, g, w), formula = bf(y | weights(w) ~ x + (1 | g)),
         family = bernoulli(), draws = 300L,
         ll = lv(function(par, d, ui) {
           d$w * dbinom(d$y, 1, plogis(par[1] + par[2] * d$x + ui), log = TRUE)
         }),
         sdof = function(par) exp(par[3]))
  }
)
# the ordinal one needs thresholds, whose chart is read from the fit's names
design_list$ordinal <- NULL
