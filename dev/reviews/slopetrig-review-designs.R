# Review of lane slopetrig: out-of-sample designs for the ratio trigger.
# Each entry: data(seed) -> data.frame; form; family; ref(d) -> the same
# model on a column moved by an affine map the model absorbs (centered,
# standardized; within the nonzero rows where the zero rows carry their
# own indicator in the random term), or NULL where no exact twin exists
# (then "short" is read against the best arm only).
# Sourced by dev/reviews/slopetrig-review-oos.R and -probe.R.
std <- function(x) (x - mean(x)) / sd(x)
base_g <- function(seed, ng = 20, per = 15) {
  set.seed(seed)
  n <- ng * per
  list(n = n, g = factor(rep(seq_len(ng), each = per)), z = rnorm(n),
       b0 = rnorm(ng, 0, 0.5), b1 = rnorm(ng, 0, 0.5))
}
designs <- list(
  binom_year = list(
    data = function(seed) {
      b <- base_g(seed)
      eta <- -0.5 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z
      data.frame(g = b$g, x = 2000 + 3 * b$z,
                 y = rbinom(b$n, 10, plogis(eta)), n = 10L)
    },
    form = y | trials(n) ~ 1 + (1 + x | g), family = binomial(),
    ref = function(d) { d$x <- std(d$x); d }),
  binom_dose = list(
    data = function(seed) {
      b <- base_g(seed)
      eta <- -0.5 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z
      data.frame(g = b$g, x = 5 + 0.1 * b$z,
                 y = rbinom(b$n, 10, plogis(eta)), n = 10L)
    },
    form = y | trials(n) ~ x + (1 + x | g), family = binomial(),
    ref = function(d) { d$x <- std(d$x); d }),
  nb_year = list(
    data = function(seed) {
      b <- base_g(seed)
      eta <- 1 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z
      data.frame(g = b$g, x = 2000 + 3 * b$z,
                 y = rnbinom(b$n, mu = exp(eta), size = 5))
    },
    form = y ~ 1 + (1 + x | g), family = negbinomial(),
    ref = function(d) { d$x <- std(d$x); d }),
  nb_dose = list(
    data = function(seed) {
      b <- base_g(seed)
      eta <- 1 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z
      data.frame(g = b$g, x = 5 + 0.1 * b$z,
                 y = rnbinom(b$n, mu = exp(eta), size = 5))
    },
    form = y ~ x + (1 + x | g), family = negbinomial(),
    ref = function(d) { d$x <- std(d$x); d }),
  # x1 far out, x2 just under the threshold (ratio about 1.94)
  two_slopes = list(
    data = function(seed) {
      b <- base_g(seed)
      z2 <- rnorm(b$n)
      b2 <- rnorm(20, 0, 0.5)
      eta <- 0.5 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z + b2[b$g] * z2
      data.frame(g = b$g, x1 = 2000 + 3 * b$z, x2 = 1 + 0.6 * z2,
                 y = eta + rnorm(b$n))
    },
    form = y ~ 1 + (1 + x1 + x2 | g), family = gaussian(),
    ref = function(d) { d$x1 <- std(d$x1); d$x2 <- std(d$x2); d }),
  # x2 alone far out (5 + 0.1 z), x1 centered
  two_slopes_b = list(
    data = function(seed) {
      b <- base_g(seed)
      z2 <- rnorm(b$n)
      b2 <- rnorm(20, 0, 0.5)
      eta <- 0.5 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z + b2[b$g] * z2
      data.frame(g = b$g, x1 = b$z, x2 = 5 + 0.1 * z2,
                 y = eta + rnorm(b$n))
    },
    form = y ~ 1 + (1 + x1 + x2 | g), family = gaussian(),
    ref = function(d) { d$x1 <- std(d$x1); d$x2 <- std(d$x2); d }),
  sigma_year = list(
    data = function(seed) {
      b <- base_g(seed, 30, 20)
      s0 <- rnorm(30, 0, 0.3)
      s1 <- rnorm(30, 0, 0.3)
      ls <- s0[b$g] + (0.2 + s1[b$g]) * b$z
      data.frame(g = b$g, x = 2000 + 3 * b$z,
                 y = 1 + b$b0[b$g] + rnorm(b$n, 0, exp(ls)))
    },
    form = bf(y ~ 1 + (1 | g), sigma ~ 1 + (1 + x | g)), family = gaussian(),
    ref = function(d) { d$x <- std(d$x); d }),
  nested_year = list(
    data = function(seed) {
      b <- base_g(seed)
      g1 <- factor(rep(1:10, each = 30))
      g2 <- factor(rep(1:3, each = 10, times = 10))
      c0 <- rnorm(10, 0, 0.5)
      c1 <- rnorm(10, 0, 0.5)
      g12 <- interaction(g1, g2)
      eta <- 0.5 + c0[g1] + c1[g1] * b$z + b$b0[as.integer(g12) %% 20 + 1] +
        b$b1[as.integer(g12) %% 20 + 1] * b$z + 0.3 * b$z
      data.frame(g1 = g1, g2 = g2, x = 2000 + 3 * b$z, y = eta + rnorm(b$n))
    },
    form = y ~ 1 + (1 + x | g1 / g2), family = gaussian(),
    ref = function(d) { d$x <- std(d$x); d }),
  mm_year = list(
    data = function(seed) {
      b <- base_g(seed)
      g1 <- sample(1:20, b$n, TRUE)
      g2 <- sample(1:20, b$n, TRUE)
      eta <- 0.5 + 0.5 * (b$b0[g1] + b$b0[g2]) +
        (0.3 + 0.5 * (b$b1[g1] + b$b1[g2])) * b$z
      data.frame(g1 = factor(g1, 1:20), g2 = factor(g2, 1:20),
                 x = 2000 + 3 * b$z, y = eta + rnorm(b$n))
    },
    form = y ~ 1 + (1 + x | mm(g1, g2)), family = gaussian(),
    ref = function(d) { d$x <- std(d$x); d }),
  # 5 + 0.1 z, with 3 data-entry outliers at 50: the ratio over every
  # row falls to about 1.6 while the bulk is the 5033 column
  outlier_dose = list(
    data = function(seed) {
      b <- base_g(seed)
      eta <- 0.5 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z
      x <- 5 + 0.1 * b$z
      x[c(7, 150, 290)] <- 50
      data.frame(g = b$g, x = x, y = eta + rnorm(b$n))
    },
    form = y ~ x + (1 + x | g), family = gaussian(),
    ref = function(d) { d$x <- std(d$x); d }),
  # a calendar year with 30 percent of rows coded 0 (ratio about 1.83)
  year_zero30 = list(
    data = function(seed) {
      b <- base_g(seed)
      eta <- 0.5 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z
      x <- 2000 + 3 * b$z
      x[seq_len(b$n) %% 10 < 3] <- 0
      data.frame(g = b$g, x = x, y = eta + rnorm(b$n))
    },
    form = y ~ 1 + (1 + x | g), family = gaussian(),
    ref = function(d) { d$x <- std(d$x); d }),
  # a factor-by-numeric column: year on treated rows, 0 on control rows
  # (ratio about 1.41 over every row)
  fac_year = list(
    data = function(seed) {
      b <- base_g(seed)
      trt <- rep(c(1, 0), length.out = b$n)
      bt <- rnorm(20, 0, 0.5)
      eta <- 0.5 + b$b0[b$g] + trt * (0.2 + bt[b$g]) +
        trt * (0.3 + b$b1[b$g]) * b$z
      data.frame(g = b$g, trt = trt, x = trt * (2000 + 3 * b$z),
                 y = eta + rnorm(b$n))
    },
    form = y ~ trt + x + (1 + trt + x | g), family = gaussian(),
    ref = function(d) {
      on <- d$trt == 1
      d$x[on] <- std(d$x[on])
      d
    }),
  # the same with a dose of 5 + 0.1 z on treated rows
  fac_dose = list(
    data = function(seed) {
      b <- base_g(seed)
      trt <- rep(c(1, 0), length.out = b$n)
      bt <- rnorm(20, 0, 0.5)
      eta <- 0.5 + b$b0[b$g] + trt * (0.2 + bt[b$g]) +
        trt * (0.3 + b$b1[b$g]) * b$z
      data.frame(g = b$g, trt = trt, x = trt * (5 + 0.1 * b$z),
                 y = eta + rnorm(b$n))
    },
    form = y ~ trt + x + (1 + trt + x | g), family = gaussian(),
    ref = function(d) {
      on <- d$trt == 1
      d$x[on] <- std(d$x[on])
      d
    }),
  grby_year = list(
    data = function(seed) {
      b <- base_g(seed)
      h <- factor(ifelse(as.integer(b$g) <= 10, "a", "b"))
      eta <- 0.5 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z
      data.frame(g = b$g, h = h, x = 2000 + 3 * b$z, y = eta + rnorm(b$n))
    },
    form = y ~ 1 + (1 + x | gr(g, by = h)), family = gaussian(),
    ref = function(d) { d$x <- std(d$x); d }),
  ar_year = list(
    data = function(seed) {
      b <- base_g(seed)
      t <- rep(1:15, 20)
      e <- unlist(lapply(1:20, function(i)
        as.numeric(stats::arima.sim(list(ar = 0.5), 15))))
      eta <- 0.5 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z
      data.frame(g = b$g, t = t, x = 2000 + 3 * b$z, y = eta + e)
    },
    form = y ~ 1 + (1 + x | g) + ar(time = t, gr = g), family = gaussian(),
    ref = function(d) { d$x <- std(d$x); d }),
  cs_year = list(
    data = function(seed) {
      b <- base_g(seed)
      eta <- 0.5 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z
      data.frame(g = b$g, x = 2000 + 3 * b$z, y = eta + rnorm(b$n))
    },
    form = y ~ 1 + cs(1 + x | g), family = gaussian(),
    ref = function(d) { d$x <- std(d$x); d }),
  # a multivariate model whose first response is fit on half the rows;
  # x is 0 on the rows that response does not use
  subset_year = list(
    data = function(seed) {
      b <- base_g(seed)
      s1 <- rep(c(TRUE, FALSE), length.out = b$n)
      eta <- 0.5 + b$b0[b$g] + (0.3 + b$b1[b$g]) * b$z
      data.frame(g = b$g, s1 = s1, x = ifelse(s1, 2000 + 3 * b$z, 0),
                 y1 = eta + rnorm(b$n), y2 = rnorm(b$n) + b$b0[b$g])
    },
    form = bf(y1 | subset(s1) ~ 1 + (1 + x | g)) + bf(y2 ~ 1 + (1 | g)) +
      set_rescor(FALSE),
    family = gaussian(),
    ref = function(d) {
      d$x[d$s1] <- std(d$x[d$s1])
      d
    })
)
