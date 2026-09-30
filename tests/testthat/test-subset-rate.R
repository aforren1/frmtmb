# brms's addition terms subset(), index(), rate() and cat(), and the
# mi(x, idx = ) predictor. The agreement with brms's own compiled
# programs is in test-brms-likelihood.R (rows 24 and 25, gated); this
# file holds the identities that need no Stan and every refusal.

eps <- .Machine$double.eps

# a joint objective and a sum of separate ones at one parameter point
# built from the joint model's own vector: the pieces share no
# parameter, so the sum is an IDENTITY, not a measurement
split_sum_at <- function(mv, parts, pm) {
  lm_ <- mv$obj$env$parList(pm)
  tot <- 0
  for (r in names(parts)) {
    f <- parts[[r]]
    pl <- f$obj$env$parList(f$opt$par)
    for (nm in names(pl)) {
      src <- names(mv$frame$par_template[[nm]])
      pl[[nm]] <- lm_[[nm]][match(paste0(r, "_",
                                         names(f$frame$par_template[[nm]])),
                                  src)]
    }
    tot <- tot + f$obj$fn(unlist(pl))
  }
  tot
}

subset_data <- function() {
  set.seed(25)
  n <- 120
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(rep(rep(letters[1:6], each = 2), n / 12)),
                  s1 = rep(c(TRUE, FALSE), n / 2),
                  s2 = c(rep(TRUE, 90), rep(FALSE, 30)))
  d$y1 <- 1 + d$x + rnorm(6, 0, 0.7)[d$g] + rnorm(n)
  d$y2 <- rpois(n, exp(0.5 - 0.3 * d$z))
  # outside y2's rows its variables may be missing, as in brms
  d$y2[!d$s2] <- NA
  d$z[!d$s2] <- NA
  d
}

# ------------------------------------------------------------------ rate

rate_data <- function() {
  set.seed(24)
  n <- 200
  d <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4))
  d$y <- rpois(n, exp(0.2 + 0.5 * d$x) * d$time)
  d$yn <- rnbinom(n, mu = exp(0.2 + 0.5 * d$x) * d$time,
                  size = 2 * d$time)
  d
}

test_that("rate() on poisson is the offset(log(denom)) model", {
  d <- rate_data()
  f1 <- frm(y | rate(time) ~ x, data = d, family = poisson())
  f2 <- frm(y ~ x + offset(log(time)), data = d, family = poisson())
  # IDENTITY: under the log link rate() adds log(denom) to the linear
  # predictor, in the same order as the offset, so the objectives agree
  # bit for bit at any point
  set.seed(1)
  p <- f1$opt$par + rnorm(length(f1$opt$par), 0, 0.1)
  expect_identical(f1$obj$fn(p), f2$obj$fn(p))
  expect_identical(as.numeric(f1$obj$gr(p)), as.numeric(f2$obj$gr(p)))
  # the expected count carries the exposure; the mu dpar does not, as
  # brms's posterior_epred() and fitted(dpar = "mu")
  ep <- fitted(f1)[, "Estimate"]
  mu <- fitted(f1, dpar = "mu")[, "Estimate"]
  expect_lt(max(abs(ep - mu * d$time)), 64 * eps * max(ep))
  expect_lt(max(abs(frm_linpred(f1, type = "link") -
                      log(mu))), 64 * eps * max(abs(log(mu))))
  # on newdata the exposure is newdata's own, and it is required
  nd <- d[1:4, ]
  nd$time <- c(1, 2, 3, 4)
  en <- fitted(f1, newdata = nd)[, "Estimate"]
  mn <- fitted(f1, newdata = nd, dpar = "mu")[, "Estimate"]
  expect_lt(max(abs(en - mn * nd$time)), 64 * eps * max(en))
  expect_error(fitted(f1, newdata = d[1:4, "x", drop = FALSE]),
               "rate\\(time\\).*newdata has no column time")
})

test_that("rate() scales the negative binomial's mean and shape", {
  d <- rate_data()
  f <- frm(yn | rate(time) ~ x, data = d, family = negbinomial())
  # brms's neg_binomial_2_log(eta + log(d), shape * d), written out
  est <- f$estimates
  eta <- as.vector(f$frame$linpreds[[1]]$X %*% est$beta)
  shape <- exp(est$betad[[1]])
  ref <- sum(stats::dnbinom(d$yn, mu = exp(eta) * d$time,
                            size = shape * d$time, log = TRUE))
  expect_lt(abs(as.numeric(logLik(f)) - ref), 1e3 * eps * abs(ref))
  # the variance the pearson residual divides by is mu d (1 + mu / shape)
  mu <- exp(eta)
  r <- residuals(f, type = "pearson")[, "Estimate"]
  v <- mu * d$time + (mu * d$time)^2 / (shape * d$time)
  expect_lt(max(abs(r - (d$yn - mu * d$time) / sqrt(v))),
            1e3 * eps * max(abs(r)))
})

test_that("rate() on geometric is the negative binomial of size denom", {
  d <- rate_data()
  f <- frm(y | rate(time) ~ x, data = d, family = geometric())
  eta <- as.vector(f$frame$linpreds[[1]]$X %*% f$estimates$beta)
  ref <- sum(stats::dnbinom(d$y, mu = exp(eta) * d$time, size = d$time,
                            log = TRUE))
  expect_lt(abs(as.numeric(logLik(f)) - ref), 1e3 * eps * abs(ref))
})

test_that("rate() is refused where brms refuses it", {
  d <- rate_data()
  # the families brms takes it for, and no others
  expect_error(frm(y | rate(time) ~ x, data = d, family = gaussian(),
                   dry_run = "frame"),
               "`rate\\(\\)` is not one this family reads")
  expect_error(frm(y | rate(time) ~ x, data = d,
                   family = zero_inflated_poisson(), dry_run = "frame"),
               "`rate\\(\\)` is not one this family reads")
  # brms's "Rate denomiators should be positive"
  d$time[3] <- 0
  expect_error(frm(y | rate(time) ~ x, data = d, family = poisson(),
                   dry_run = "frame"),
               "rate denominators should be positive")
  # and the case where the guard's condition is absent: a positive
  # denominator of every size fits
  d$time[3] <- 1e-3
  expect_s3_class(frm(y | rate(time) ~ x, data = d, family = poisson(),
                      dry_run = "frame"), "frmtmb_frame")
})

test_that("predict() and simulate() draw around mu times the exposure", {
  d <- rate_data()
  f <- frm(y | rate(time) ~ x, data = d, family = poisson())
  # a large exposure moves the draws by its factor
  nd <- d[1:2, ]
  nd$time <- c(1, 1000)
  set.seed(11)
  pr <- predict(f, newdata = nd, ndraws = 400)
  mu <- fitted(f, newdata = nd, dpar = "mu")[, "Estimate"]
  expect_gt(pr[2, "Estimate"] / mu[2], 900)
  expect_lt(pr[1, "Estimate"] / mu[1], 3)
  s <- simulate(f, nsim = 200, seed = 3)
  expect_gt(cor(rowMeans(s), fitted(f)[, "Estimate"]), 0.95)
})

test_that("newdata's exposure must be positive, as the data's must", {
  d <- rate_data()
  f <- frm(y | rate(time) ~ x, data = d, family = poisson())
  nd <- d[1:3, ]
  # brms's "Rate denomiators should be positive." on newdata too: a zero
  # or negative exposure would give a non-positive expected count
  for (bad in c(0, -2)) {
    nd$time[2] <- bad
    expect_error(fitted(f, newdata = nd),
                 "rate\\(time\\) on newdata: rate denominators should")
    expect_error(predict(f, newdata = nd, ndraws = 10),
                 "rate\\(time\\) on newdata: rate denominators should")
  }
  # the case where the guard's condition is absent: a small positive
  # exposure predicts, and the prediction is mu times it
  nd$time[2] <- 1e-3
  e <- fitted(f, newdata = nd)[, "Estimate"]
  m <- fitted(f, newdata = nd, dpar = "mu")[, "Estimate"]
  expect_lt(max(abs(e - m * nd$time)), 64 * eps * max(e))
})

# ------------------------------------------------------------------- cat

test_that("cat(x) is brms's deprecated thres(x - 1)", {
  set.seed(4)
  dc <- data.frame(s = sample(1:5, 60, TRUE), x = rnorm(60))
  fc <- allow_warnings(
    frm(s | cat(6) ~ x, data = dc, family = cumulative(),
        dry_run = "frame"),
    "Addition argument 'cat' is deprecated",
    require = "Addition argument 'cat' is deprecated")
  ft <- frm(s | thres(5) ~ x, data = dc, family = cumulative(),
            dry_run = "frame")
  expect_identical(fc$par_template$tau_raw, ft$par_template$tau_raw)
  expect_length(fc$par_template$tau_raw, 5L)
  # brms's resp_cat(x) takes its argument by name too
  fx <- allow_warnings(
    frm(s | cat(x = 6) ~ x, data = dc, family = cumulative(),
        dry_run = "frame"),
    "Addition argument 'cat' is deprecated")
  expect_identical(fx$par_template$tau_raw, ft$par_template$tau_raw)
  expect_error(frm(s | cat(k = 6) ~ x, data = dc, family = cumulative(),
                   dry_run = "frame"),
               "cat\\(\\) takes the number of categories")
  # the two spellings set one datum, which brms lets cat() win silently,
  # in either order
  allow_warnings(
    expect_error(frm(s | cat(6) + thres(5) ~ x, data = dc,
                     family = cumulative(), dry_run = "frame"),
                 "number of thresholds is set twice"),
    "Addition argument 'cat' is deprecated")
  expect_error(frm(s | thres(5) + cat(6) ~ x, data = dc,
                   family = cumulative(), dry_run = "frame"),
               "number of thresholds is set twice")
  # two thres() keep their own refusal, which cat() does not touch
  expect_error(frm(s | thres(5) + thres(gr = x) ~ x, data = dc,
                   family = cumulative(), dry_run = "frame"),
               "Duplicated addition term")
  # an ordinal family's term, like thres()
  allow_warnings(
    expect_error(frm(s | cat(6) ~ x, data = dc, family = gaussian(),
                     dry_run = "frame"),
                 "not a valid addition term for family 'gaussian'"),
    "Addition argument 'cat' is deprecated")
})

# ---------------------------------------------------------------- subset

test_that("subset() fits each response on its own rows", {
  d <- subset_data()
  fr <- frm(bf(y1 | subset(s1) ~ x + (1 | g)) + gaussian() +
              bf(y2 | subset(s2) ~ z) + poisson(), data = d,
            dry_run = "frame")
  # brms's N_y1, N_y2, and the design rows; nobs() is the whole frame,
  # as brms's
  expect_identical(length(fr$y$y1), sum(d$s1))
  expect_identical(length(fr$y$y2), sum(d$s2))
  expect_identical(nrow(fr$linpreds[["y2.mu"]]$X), sum(d$s2))
  expect_identical(nrow(fr$linpreds[["y1.mu"]]$Z), sum(d$s1))
  expect_identical(fr$n_obs, nrow(d))
  expect_identical(fr$y$y2, as.numeric(d$y2[d$s2]))
})

test_that("a subset() model is the sum of the separate fits", {
  d <- subset_data()
  mv <- frm(bf(y1 | subset(s1) ~ x) + gaussian() +
              bf(y2 | subset(s2) ~ z) + poisson(), data = d)
  u <- list(y1 = frm(y1 ~ x, data = d[d$s1, ], family = gaussian()),
            y2 = frm(y2 ~ z, data = d[d$s2, ], family = poisson()))
  set.seed(2)
  pm <- mv$opt$par + rnorm(length(mv$opt$par), 0, 0.1)
  tot <- split_sum_at(mv, u, pm)
  expect_lt(abs(mv$obj$fn(pm) - tot), 64 * eps * abs(tot))
  # brms's nobs(): the data's rows, and a response's own rows with resp
  expect_identical(nobs(mv), nrow(d))
  expect_identical(nobs(mv, resp = "y1"), sum(d$s1))
  expect_identical(nobs(mv, resp = "y2"), sum(d$s2))
  expect_error(nobs(mv, resp = "y3"), "y3")
})

test_that("an NA drops its row only where a response uses it", {
  d <- subset_data()
  # z is NA on the 30 rows outside y2's subset, and nothing else uses z:
  # every row stays, so y1 keeps all of its own
  fr <- frm(bf(y1 | subset(s1) ~ x) + gaussian() +
              bf(y2 | subset(s2) ~ z) + poisson(), data = d,
            dry_run = "frame")
  expect_identical(fr$n_obs, nrow(d))
  # the case where the rule does NOT apply: an NA on a row inside y2's
  # subset drops that row, as brms's na_omit() drops it
  d$z[2] <- NA
  fr2 <- suppressMessages(
    frm(bf(y1 | subset(s1) ~ x) + gaussian() +
          bf(y2 | subset(s2) ~ z) + poisson(), data = d,
        dry_run = "frame"))
  expect_identical(fr2$n_obs, nrow(d) - 1L)
  expect_identical(length(fr2$y$y2), sum(d$s2) - 1L)
  # a variable both responses use: exempt on a row that BOTH leave out
  # (row 120 is outside s1 and s2), dropped on a row that one of them
  # uses (row 1 is inside both)
  d <- subset_data()
  d$x[c(1, nrow(d))] <- NA
  fr3 <- suppressMessages(
    frm(bf(y1 | subset(s1) ~ x) + gaussian() +
          bf(y2 | subset(s2) ~ x) + poisson(), data = d,
        dry_run = "frame"))
  expect_identical(fr3$n_obs, nrow(d) - 1L)
  expect_identical(length(fr3$y$y1), sum(d$s1) - 1L)
})

test_that("subset() in a univariate model is a row filter", {
  d <- subset_data()
  fu <- frm(y1 | subset(s1) ~ x, data = d, family = gaussian())
  fv <- frm(y1 ~ x, data = d[d$s1, ], family = gaussian())
  expect_identical(fu$obj$fn(fv$opt$par), fv$obj$fn(fv$opt$par))
  # brms's nobs() counts the data's rows, nrow(model.frame(fit)) there,
  # while the likelihood sums over the subset's
  expect_identical(nobs(fu), nrow(d))
  expect_identical(nobs(fu, resp = "y1"), nrow(d))
  expect_identical(attr(logLik(fu), "nobs"), sum(d$s1))
  # simulate() still draws one value per fitted row
  expect_identical(nrow(simulate(fu, nsim = 2, seed = 1)), sum(d$s1))
  # on newdata, the rows where the subset is TRUE, as brms keeps them
  expect_identical(nrow(fitted(fu, newdata = d[1:6, ])), 3L)
})

test_that("post-fit methods take one response at a time, as brms asks", {
  d <- subset_data()
  mv <- frm(bf(y1 | subset(s1) ~ x) + gaussian() +
              bf(y2 | subset(s2) ~ z) + poisson(), data = d)
  msg <- "argument 'resp' must be a single variable name"
  expect_error(fitted(mv), msg)
  expect_error(predict(mv), msg)
  expect_error(frm_linpred(mv), msg)
  expect_error(fitted(mv, resp = c("y1", "y2")), msg)
  expect_identical(dim(fitted(mv, resp = "y1")), c(sum(d$s1), 4L))
  expect_identical(dim(fitted(mv, resp = "y2")), c(sum(d$s2), 4L))
  expect_identical(dim(predict(mv, resp = "y2", ndraws = 50)),
                   c(sum(d$s2), 4L))
  # newdata: the response's own subset picks the rows, and only that
  # response's subset variable is needed
  nd <- d[1:8, ]
  expect_identical(nrow(fitted(mv, newdata = nd, resp = "y1")),
                   sum(nd$s1))
  nd$s1 <- NULL
  expect_identical(nrow(fitted(mv, newdata = nd, resp = "y2")),
                   sum(nd$s2))
  expect_error(fitted(mv, newdata = nd, resp = "y1"),
               "newdata has no column s1")
  # the fit's character rule holds on newdata too
  nd$s1 <- as.character(d$s1[1:8])
  expect_error(fitted(mv, newdata = nd, resp = "y1"),
               "on newdata: a character subset variable is refused")
  # the in-sample fitted values are y1's own rows times its coefficients
  b <- fixef(mv)[c("y1_Intercept", "y1_x"), "Estimate"]
  e1 <- fitted(mv, resp = "y1")[, "Estimate"]
  expect_lt(max(abs(e1 - (b[[1]] + b[[2]] * d$x[d$s1]))),
            64 * eps * max(abs(e1)))
})

test_that("subset() is refused where it cannot be followed", {
  d <- subset_data()
  fr <- function(f, ...) frm(f, data = d, family = gaussian(),
                              dry_run = "frame", ...)
  expect_error(fr(bf(y1 | subset(s1) ~ x) + bf(y2 | subset(s2) ~ z) +
                    set_rescor(TRUE)),
               "cannot be combined with rescor = TRUE")
  expect_error(fr(bf(y1 | subset(s1) ~ x) + bf(y2 ~ x),
                  na.action = na.exclude),
               "cannot be combined with na.action = na.exclude")
  d$sm <- 0.5
  expect_error(fr(bf(y1 | subset(s1) ~ me(x, sm)) + bf(y2 ~ x)),
               "'subset' is not supported when using 'me' terms")
  d$s3 <- d$s1
  d$s3[5] <- NA
  expect_error(fr(bf(y1 | subset(s3) ~ x) + bf(y2 ~ x)),
               "subset variables may not contain NAs")
  # brms refuses a character subset ("invalid argument type"), which
  # as.logical() would read; a factor of the same values is accepted
  d$s5 <- as.character(d$s1)
  expect_error(fr(bf(y1 | subset(s5) ~ x) + bf(y2 ~ x)),
               "a character subset variable is refused")
  d$s6 <- factor(d$s5)
  expect_s3_class(fr(bf(y1 | subset(s6) ~ x) + bf(y2 | subset(s2) ~ z)),
                  "frmtmb_frame")
  d$s4 <- FALSE
  expect_error(fr(bf(y1 | subset(s4) ~ x) + bf(y2 | subset(s2) ~ z)),
               "leaves it no rows")
  # |ID|-linked terms whose rows carry different levels
  d$g2 <- factor(ifelse(d$s2, as.character(d$g), "h"))
  expect_error(fr(bf(y1 | subset(s1) ~ x + (1 | p | g2)) +
                    bf(y2 | subset(s2) ~ z + (1 | p | g2))),
               "A response with subset\\(\\) has the levels")
  # the case where the guard's condition is absent: the same |ID| terms
  # on rows that carry the same levels are accepted
  expect_s3_class(fr(bf(y1 | subset(s1) ~ x + (1 | p | g)) +
                       bf(y2 | subset(s2) ~ z + (1 | p | g))),
                  "frmtmb_frame")
})

# ------------------------------------------------------ index / mi(idx)

mi_idx_data <- function() {
  set.seed(26)
  n <- 120
  dm <- data.frame(g1 = sample(seq(1, n - 1, 2), n, TRUE), g2 = seq_len(n),
                   s = rep(c(TRUE, FALSE), n / 2), w = rnorm(n))
  dm$x <- rnorm(n)
  dm$y <- 1 + 0.5 * dm$x[match(dm$g1, dm$g2)] + 0.3 * dm$w +
    rnorm(n, sd = 0.5)
  dm
}

test_that("mi(x, idx = ) reads the rows of x that index() names", {
  dm <- mi_idx_data()
  dm$x[c(3, 7, 8)] <- NA
  fr <- frm(bf(y ~ mi(x, idx = g1) + w) +
              bf(x | mi() + index(g2) + subset(s) ~ 1),
            data = dm, family = gaussian(), dry_run = "frame")
  # brms's idxl: match(idx, index(x) on x's rows)
  xs <- dm$g2[dm$s]
  expect_identical(fr$linpreds[["y.mu"]]$mi[[1]]$idxl, match(dm$g1, xs))
  # x's latent values are its NAs inside the subset; row 8 is outside
  expect_identical(as.integer(fr$mi_map$x$rows), which(is.na(dm$x[dm$s])))
  # brms's coefficient name, rename("mi(x, idx = g1)")
  expect_true("mixidxEQg1" %in% colnames(fr$linpreds[["y.mu"]]$X))
})

test_that("with x observed, mi(x, idx = ) is the regression on x[idx]", {
  dm <- mi_idx_data()
  mv <- frm(bf(y ~ mi(x, idx = g1) + w) +
              bf(x | mi() + index(g2) + subset(s) ~ 1),
            data = dm, family = gaussian())
  dm$xm <- dm$x[match(dm$g1, dm$g2)]
  fy <- frm(y ~ xm + w, data = dm, family = gaussian())
  fx <- frm(x ~ 1, data = dm[dm$s, ], family = gaussian())
  # IDENTITY at a shared point: with nothing to impute the joint
  # density factorizes into the two regressions
  set.seed(5)
  pm <- mv$opt$par + rnorm(length(mv$opt$par), 0, 0.1)
  lm_ <- mv$obj$env$parList(pm)
  bn <- names(mv$frame$par_template$beta)
  py <- fy$obj$env$parList(fy$opt$par)
  py$beta <- lm_$beta[match(c("y_(Intercept)", "y_mixidxEQg1", "y_w"), bn)]
  py$betad <- lm_$betad[1]
  px <- fx$obj$env$parList(fx$opt$par)
  px$beta <- lm_$beta[match("x_(Intercept)", bn)]
  px$betad <- lm_$betad[2]
  tot <- fy$obj$fn(unlist(py)) + fx$obj$fn(unlist(px))
  expect_lt(abs(mv$obj$fn(pm) - tot), 64 * eps * abs(tot))
  # on newdata the rows are matched within newdata, as brms matches
  # them: the whole data again gives the in-sample values
  e_in <- fitted(mv, resp = "y")[, "Estimate"]
  e_nd <- fitted(mv, newdata = dm, resp = "y")[, "Estimate"]
  expect_lt(max(abs(e_in - e_nd)), 64 * eps * max(abs(e_in)))
  # and ten rows hold too few of x's rows to match every idx, which
  # brms refuses with the same words
  expect_error(fitted(mv, newdata = dm[1:10, ], resp = "y"),
               "Could not match all indices in response 'x'")
})

test_that("mi(x, idx = ) without a subset is mi(x) when idx names the row", {
  dm <- mi_idx_data()
  dm$x[c(3, 7)] <- NA
  a <- frm(bf(y ~ mi(x) + w) + bf(x | mi() ~ 1), data = dm,
           family = gaussian())
  b <- frm(bf(y ~ mi(x, idx = g2) + w) + bf(x | mi() + index(g2) ~ 1),
           data = dm, family = gaussian())
  set.seed(6)
  p <- a$opt$par + rnorm(length(a$opt$par), 0, 0.1)
  expect_identical(a$obj$fn(p), b$obj$fn(p))
})

test_that("mi() across responses with different rows is refused as brms", {
  dm <- mi_idx_data()
  fr <- function(f) frm(f, data = dm, family = gaussian(),
                        dry_run = "frame")
  expect_error(fr(bf(y ~ mi(x, idx = g1)) + bf(x | mi() + subset(s) ~ 1)),
               "Response 'x' needs to have an 'index' addition term")
  expect_error(fr(bf(y ~ mi(x)) + bf(x | mi() + subset(s) + index(g2) ~ 1)),
               "mi() terms of subsetted variables require the 'idx'",
               fixed = TRUE)
  expect_error(fr(bf(y | mi() ~ mi(x, idx = g1)) +
                    bf(x | mi() + subset(s) + index(g2) ~ mi(y))),
               "mi() terms in subsetted formulas require the 'idx'",
               fixed = TRUE)
  dm$g1[2] <- 2
  expect_error(fr(bf(y ~ mi(x, idx = g1)) +
                    bf(x | mi() + index(g2) + subset(s) ~ 1)),
               "Could not match all indices in response 'x'")
  dm <- mi_idx_data()
  dm$g2[3] <- 1
  expect_error(fr(bf(y ~ mi(x, idx = g1)) +
                    bf(x | mi() + index(g2) + subset(s) ~ 1)),
               "Index of response 'x' contains duplicated values")
  expect_error(fr(bf(y ~ mi(x, idx = g1 + 1)) + bf(x | mi() ~ 1)),
               "idx takes one variable name")
})
