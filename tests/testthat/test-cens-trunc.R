test_that("right-censored gaussian matches survival::survreg", {
  skip_if_not_installed("survival")
  set.seed(101)
  n <- 400
  x <- rnorm(n)
  ystar <- 1 + 0.8 * x + rnorm(n, 0, 1.2)
  cpoint <- 2
  y <- pmin(ystar, cpoint)
  cen <- as.numeric(ystar > cpoint)   # 1 = right-censored
  dd <- data.frame(y = y, x = x, cen = cen)

  fit <- frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd)
  ref <- survival::survreg(survival::Surv(y, 1 - cen) ~ x,
                           data = dd, dist = "gaussian")
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-5)
  expect_vector_equal(fixef(fit)$mu, unname(coef(ref)), tol = 1e-4)
  expect_lt(abs(exp(fixef(fit)$sigma[[1]]) - ref$scale), 1e-3)
})

test_that("mixed left/right censoring matches a hand-rolled reference", {
  set.seed(102)
  n <- 400
  x <- rnorm(n)
  ystar <- 0.5 * x + rnorm(n)
  lo <- -1; hi <- 1.5
  y <- pmin(pmax(ystar, lo), hi)
  cen <- ifelse(ystar < lo, -1, ifelse(ystar > hi, 1, 0))
  dd <- data.frame(y = y, x = x, cen = cen)
  fit <- frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd)

  nll_ref <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    mu <- p$b[1] + p$b[2] * x
    s <- exp(p$ls)
    ll <- RTMB::dnorm(y, mu, s, log = TRUE)
    Fv <- RTMB::pnorm((y - mu) / s)
    ir <- which(cen == 1); il <- which(cen == -1)
    ll[ir] <- log(1 - Fv[ir])
    ll[il] <- log(Fv[il])
    -sum(ll)
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), ls = 0),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
})

test_that("truncation matches a hand-rolled reference", {
  set.seed(103)
  n <- 2000
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n)
  keep <- y > 0
  dd <- data.frame(y = y[keep], x = x[keep])

  fit <- frm(bf(y | trunc(lb = 0) ~ x) + gaussian(), data = dd)
  yv <- dd$y; xv <- dd$x
  nll_ref <- function(p) {
    mu <- p$b[1] + p$b[2] * xv
    s <- exp(p$ls)
    -sum(RTMB::dnorm(yv, mu, s, log = TRUE) -
           log(1 - RTMB::pnorm((0 - mu) / s)))
  }
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), ls = 0),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
  # truncation-corrected fit recovers the latent coefficients
  expect_vector_equal(fixef(fit)$mu, c(1, 0.5), tol = 0.15)
})

test_that("brms character and factor censoring codes fit identically", {
  set.seed(104)
  n <- 300
  x <- rnorm(n)
  ystar <- 0.5 * x + rnorm(n)
  lo <- -1; hi <- 1.5
  y <- pmin(pmax(ystar, lo), hi)
  cen <- ifelse(ystar < lo, -1, ifelse(ystar > hi, 1, 0))
  lab <- c("left", "none", "right")[cen + 2]
  dd <- data.frame(y = y, x = x, cen = cen, lab = lab,
                   flab = factor(lab),
                   # prefix matching, case-insensitively
                   pre = c("LEFT", "n", "r")[cen + 2])

  fit_num <- frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd)
  for (nm in c("lab", "flab", "pre")) {
    fo <- stats::as.formula(paste0("y | cens(", nm, ") ~ x"))
    fit <- frm(bf(fo) + gaussian(), data = dd)
    expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(fit_num))),
              1e-8)
    expect_vector_equal(fixef(fit)$mu, fixef(fit_num)$mu, tol = 1e-8)
    expect_vector_equal(fixef(fit)$sigma, fixef(fit_num)$sigma, tol = 1e-8)
  }
})

test_that("string interval censoring matches its numeric code", {
  set.seed(105)
  n <- 200
  x <- rnorm(n)
  y <- 0.4 * x + rnorm(n)
  iv <- y > 0.5
  dd <- data.frame(y = y, x = x, y2 = y + 1,
                   cs = ifelse(iv, "interval", "none"),
                   cn = ifelse(iv, 2, 0))
  fit_s <- frm(bf(y | cens(cs, y2) ~ x) + gaussian(), data = dd)
  fit_n <- frm(bf(y | cens(cn, y2) ~ x) + gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(fit_s)) - as.numeric(logLik(fit_n))), 1e-8)
  expect_vector_equal(fixef(fit_s)$mu, fixef(fit_n)$mu, tol = 1e-8)
})

test_that("cens/trunc validation", {
  dd <- data.frame(y = rpois(50, 3), x = rnorm(50), cen = 0)
  expect_error(frm(bf(y | cens(cen) ~ x) + poisson(), data = dd),
               "discrete")
  dd2 <- data.frame(y = rnorm(50), x = rnorm(50), cen = 3)
  expect_error(frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd2),
               "codes")
  dd3 <- data.frame(y = rnorm(50), x = rnorm(50), cen = 2)
  expect_error(frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd3),
               "needs upper bounds")
  expect_error(frm(bf(y | trunc(0) ~ x) + gaussian(),
                   data = NULL, dry_run = "spec"),
               "named bounds")
})

test_that("unusable censoring labels error informatively", {
  dd <- data.frame(y = rnorm(50), x = rnorm(50),
                   bad = rep(c("right", "sideways"), 25),
                   amb = "",
                   iv = "interval")
  expect_error(frm(bf(y | cens(bad) ~ x) + gaussian(), data = dd),
               "sideways")
  expect_error(frm(bf(y | cens(bad) ~ x) + gaussian(), data = dd),
               "interval")
  # an empty label prefix-matches every spelling, so it is not a code
  expect_error(frm(bf(y | cens(amb) ~ x) + gaussian(), data = dd),
               "cannot decode")
  # "interval" still needs the upper bounds, same as the numeric code
  expect_error(frm(bf(y | cens(iv) ~ x) + gaussian(), data = dd),
               "needs upper bounds")
  # factor labels go through the same decoder
  dd$fbad <- factor(dd$bad)
  expect_error(frm(bf(y | cens(fbad) ~ x) + gaussian(), data = dd),
               "sideways")
})
