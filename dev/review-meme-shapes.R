# Review of lane meme: the extended shapes, shared by the brms and the
# frmtmb arms. Sourced; defines rv_data() and rv_shapes(bf, mvbf, ...).
rv_data <- function(seed = 2028, n = 40) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                  sx = 0.3, sz = runif(n, 0.2, 0.4), sw = 0.35,
                  g = factor(rep(1:8, each = 5)), h = factor(rep(1:5, 8)))
  d$xg <- rnorm(8)[d$g]
  d$wg <- rnorm(8)[d$g]
  d$sg <- 0.25
  d$y <- 1 + 0.5 * d$x - 0.3 * d$z + rnorm(n)
  d$y2 <- 0.5 * d$x + rnorm(n)
  d$v <- 0.4 * d$x + rnorm(n)
  d$v[c(3, 17, 29)] <- NA
  d$sy <- 0.2
  d
}

# each shape: formula builder, a user prior set (as a list of rows
# class/coef/group/prior) and the rows that need a correlated group
rv_shapes <- function(bf, mvbf, set_mecor, set_rescor, mi_ok = TRUE) {
  s <- list()
  s$mixed3 <- list(
    f = bf(y ~ me(xg, sg, gr = g) + me(z, sz) + me(wg, sg, gr = g)),
    pr = list(c("meanme", "", "", "normal(0, 2)"),
              c("meanme", "mez", "", "normal(1, 0.5)"),
              c("sdme", "", "", "exponential(1)"),
              c("corme", "", "g", "lkj(2)")))
  s$mixed4 <- list(
    f = bf(y ~ me(xg, sg, gr = g) + me(z, sz) + me(wg, sg, gr = g) +
             me(x, sx)),
    pr = list(c("meanme", "", "", "normal(0, 2)"),
              c("sdme", "", "", "exponential(1)"),
              c("sdme", "mewg", "", "exponential(2)"),
              c("corme", "", "", "lkj(3)"),
              c("corme", "", "g", "lkj(2)")))
  s$inter <- list(
    f = bf(y ~ me(x, sx) * z + me(w, sw)),
    pr = list(c("meanme", "", "", "normal(0, 2)"),
              c("sdme", "", "", "exponential(1)"),
              c("corme", "", "", "lkj(2)")))
  s$dpar <- list(
    f = bf(y ~ me(x, sx), sigma ~ me(z, sz)),
    pr = list(c("meanme", "", "", "normal(0, 2)"),
              c("sdme", "", "", "exponential(1)"),
              c("corme", "", "", "lkj(2)")))
  s$mv <- list(
    f = mvbf(bf(y ~ me(x, sx)), bf(y2 ~ me(x, sx) + me(z, sz))) +
      set_rescor(TRUE),
    pr = list(c("meanme", "", "", "normal(0, 2)"),
              c("sdme", "", "", "exponential(1)"),
              c("corme", "", "", "lkj(2)")))
  s$nosdx <- list(f = bf(y ~ me(x)), pr = list())
  if (mi_ok) {
    s$mi <- list(
      f = bf(y ~ me(x, sx) + mi(v)) + bf(v | mi() ~ z) + set_rescor(FALSE),
      pr = list(c("meanme", "", "", "normal(0, 2)"),
                c("sdme", "", "", "exponential(1)")))
    s$misdy <- list(
      f = bf(y | mi(sy) ~ me(x, sx) + me(z, sz)),
      pr = list(c("meanme", "", "", "normal(0, 2)"),
                c("sdme", "", "", "exponential(1)"),
                c("corme", "", "", "lkj(2)")))
  }
  s
}
