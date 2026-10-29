## Review rescorpp: data sets of the review's own fits.
rv_data <- function(w) {
  switch(w,
    mix = {
      set.seed(31); n <- 150
      x <- rnorm(n)
      data.frame(x = x, y1 = 1 + x + rnorm(n),
                 cnt = rpois(n, exp(0.5 + 0.3 * x)))
    },
    us = {
      set.seed(34); n <- 100
      x <- rnorm(n)
      E <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .7, .7, 1), 2))
      data.frame(x = x, y_a = 1 + x + E[, 1], y_b = 2 * E[, 2])
    },
    ord = {
      set.seed(33); n <- 150
      x <- rnorm(n)
      data.frame(x = x, o1 = cut(x + rnorm(n), c(-Inf, -0.5, 0.5, Inf),
                                 labels = FALSE), y2 = 0.5 * x + rnorm(n))
    },
    re = {
      set.seed(32); n <- 240; g <- rep(1:12, each = 20)
      x <- rnorm(n)
      u <- MASS::mvrnorm(12, c(0, 0), matrix(c(1, .9, .9, 1), 2))
      e <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .6, .6, 1), 2))
      data.frame(x = x, g = g, y1 = 1 + x + u[g, 1] + e[, 1],
                 y2 = -1 + 0.5 * x + u[g, 2] + e[, 2])
    })
}
rv_brms_form <- function(w) {
  switch(w,
    mix = brms::mvbf(brms::bf(y1 ~ x), brms::bf(cnt ~ x, family = poisson()),
                     rescor = FALSE),
    us = brms::mvbf(brms::bf(y_a ~ x), brms::bf(y_b ~ x), rescor = TRUE),
    ord = brms::mvbf(brms::bf(o1 ~ x, family = brms::cumulative()),
                     brms::bf(y2 ~ x, family = gaussian()), rescor = FALSE),
    re = brms::mvbf(brms::bf(y1 ~ x + (1 | p | g)),
                    brms::bf(y2 ~ x + (1 | p | g)), rescor = TRUE))
}
rv_frm_form <- function(w) {
  switch(w,
    mix = mvbf(bf(y1 ~ x), bf(cnt ~ x, family = poisson()), rescor = FALSE),
    us = mvbf(bf(y_a ~ x), bf(y_b ~ x), rescor = TRUE),
    ord = mvbf(bf(o1 ~ x, family = cumulative()),
               bf(y2 ~ x, family = gaussian()), rescor = FALSE),
    re = mvbf(bf(y1 ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g)),
              rescor = TRUE))
}
