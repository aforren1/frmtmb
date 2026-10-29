## Lane rescorpp: the data sets every comparison script reads.
rescorpp_data <- function(which) {
  switch(which,
    g3 = {
      # three responses, constant sigma: resp subsetting of a correlated set
      set.seed(11); n <- 200
      x <- rnorm(n)
      R <- matrix(c(1, .8, .5, .8, 1, .3, .5, .3, 1), 3)
      S <- diag(c(1, 2, 0.5)) %*% R %*% diag(c(1, 2, 0.5))
      e <- MASS::mvrnorm(n, c(0, 0, 0), S)
      data.frame(x = x, y1 = 1 + x + e[, 1], y2 = -1 + 0.5 * x + e[, 2],
                 y3 = 0.3 * x + e[, 3])
    },
    gs = {
      # two responses, sigma of y1 predicted: Sigma varies by row
      set.seed(12); n <- 200
      x <- rnorm(n)
      s1 <- exp(0.5 * x)
      z <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .7, .7, 1), 2))
      data.frame(x = x, y1 = 1 + x + s1 * z[, 1], y2 = -1 + 0.5 * x + z[, 2])
    },
    st = {
      # two responses, multivariate t with nu = 5
      set.seed(13); n <- 200
      x <- rnorm(n)
      z <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .6, .6, 1), 2))
      w <- sqrt(rchisq(n, 5) / 5)
      data.frame(x = x, y1 = 1 + x + z[, 1] / w,
                 y2 = -1 + 0.5 * x + 1.5 * z[, 2] / w)
    })
}
rescorpp_forms <- list(
  g3 = quote(brms::mvbf(brms::bf(y1 ~ x), brms::bf(y2 ~ x), brms::bf(y3 ~ x)) +
               brms::set_rescor(TRUE)),
  gs = quote(brms::mvbf(brms::bf(y1 ~ x, sigma ~ x), brms::bf(y2 ~ x)) +
               brms::set_rescor(TRUE)),
  st = quote(brms::mvbf(brms::bf(y1 ~ x), brms::bf(y2 ~ x)) +
               brms::set_rescor(TRUE)),
  g2f = quote(brms::mvbf(brms::bf(y1 ~ x), brms::bf(y2 ~ x)) +
                brms::set_rescor(FALSE)))
rescorpp_frm_forms <- list(
  g3 = quote(mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ x)) + set_rescor(TRUE)),
  gs = quote(mvbf(bf(y1 ~ x, sigma ~ x), bf(y2 ~ x)) + set_rescor(TRUE)),
  st = quote(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE)),
  g2f = quote(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(FALSE)))
rescorpp_family <- function(which) {
  if (identical(which, "st")) "student" else "gaussian"
}
rescorpp_dname <- function(which) if (identical(which, "g2f")) "g3" else which
