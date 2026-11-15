# Response coding, refits and the refusals they lifted. Each test was seen
# to fail on frmtmb 0.82.0 (dev/refits-findings.md has the output).

test_that("a character ordinal response is refused in brms's words", {
  set.seed(2501)
  d <- data.frame(x = stats::rnorm(60))
  lab <- c("none", "mild", "moderate", "severe")
  idx <- pmin(4, pmax(1, round(2.5 + d$x + stats::rnorm(60, 0, 0.7))))
  d$y <- lab[idx]
  # brms 2.23.0 refuses labels and digit text alike with this sentence;
  # 0.82.0 stopped on labels with "missing value where TRUE/FALSE
  # needed" and fitted the digits
  d$yn <- as.character(idx)
  for (fam in list(cumulative(), sratio())) {
    for (v in c("y", "yn")) {
      expect_error(
        frm(bf(stats::as.formula(paste(v, "~ x"))), family = fam, data = d),
        paste0("Family '", fam$family, "' requires either positive ",
               "integers or ordered factors as responses. '", v,
               "' is a character vector"),
        fixed = TRUE, class = "frmtmb_error")
    }
  }
  # the absent case: the same labels as an ordered factor, and the codes
  # as integers, fit the same model
  d$yo <- factor(d$y, levels = lab, ordered = TRUE)
  d$yi <- as.integer(idx)
  fo <- frm(bf(yo ~ x), family = cumulative(), data = d)
  fi <- frm(bf(yi ~ x), family = cumulative(), data = d)
  expect_identical(unname(fixef(fo)[, 1]), unname(fixef(fi)[, 1]))
  expect_length(fixef(fo)[, 1], 4L)
})

test_that("fitted(newdata = ) of a multivariate ordinal fit takes one row", {
  set.seed(77)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$yo <- findInterval(d$x + stats::rlogis(n), c(-0.5, 0.8)) + 1L
  d$yo2 <- findInterval(d$z + stats::rlogis(n), c(-0.3, 0.6)) + 1L
  mv <- frm(bf(yo ~ x) + bf(yo2 ~ z), family = sratio(), data = d)
  nd <- data.frame(x = c(0.2, 1), z = c(-0.4, 0))
  # 0.82.0: "values must be length 1, but FUN(X[[1]]) result is length 0"
  one <- fitted(mv, newdata = nd[1, ])
  two <- fitted(mv, newdata = nd)
  expect_identical(dim(one), c(1L, 4L, 6L))
  expect_identical(dimnames(one)[2:3], dimnames(two)[2:3])
  expect_identical(one[1, , ], two[1, , ])
})

test_that("y ~ . - w refits on data without w, as in brms", {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(80), z = stats::rnorm(80),
                  w = stats::rnorm(80))
  d$y <- d$x + stats::rnorm(80)
  f <- frm(y ~ . - w, data = d)
  expect_false("w" %in% all.vars(getCall(f)$formula))
  d5 <- d
  d5$w <- NULL
  # 0.82.0: "The model uses `w`, which is not a column of `data`"
  expect_identical(fixef(update(f, newdata = d5)), fixef(f))
  expect_identical(rownames(fixef(f)), c("Intercept", "x", "z"))
  # the absent case: a removed column the formula reads elsewhere stays
  f2 <- frm(y ~ . - w + w:x, data = d)
  expect_true("w" %in% all.vars(getCall(f2)$formula))
  expect_error(update(f2, newdata = d5), "The model uses `w`",
               fixed = TRUE, class = "frmtmb_error")
})

test_that("refit() takes a bernoulli response's own values", {
  set.seed(9)
  d <- data.frame(x = stats::rnorm(100))
  d$yb <- ifelse(stats::runif(100) < stats::plogis(d$x), -1, -2)
  fb <- frm(yb ~ x, family = bernoulli(), data = d)
  codes <- as.numeric(d$yb == -1)
  # 0.82.0: "recode newresp to 0 and 1"
  expect_identical(fixef(refit(fb, newresp = d$yb)),
                   fixef(refit(fb, newresp = codes)))
  d$yf <- factor(ifelse(d$yb == -1, "yes", "no"))
  ff <- frm(yf ~ x, family = bernoulli(), data = d)
  expect_identical(fixef(refit(ff, newresp = d$yf)),
                   fixef(refit(ff, newresp = codes)))
  expect_error(refit(fb, newresp = replace(d$yb, 1, 5)),
               "'5', which is neither", fixed = TRUE,
               class = "frmtmb_error")
  expect_error(refit(fb, newresp = replace(d$yb, 1, 0)),
               "mixes the codes with the values", fixed = TRUE,
               class = "frmtmb_error")
})

test_that("refit() reads a number one way on a fit coded 1 and 2", {
  set.seed(9)
  n <- 40
  x <- stats::rnorm(n)
  s <- stats::rbinom(n, 1, stats::plogis(0.3 + x))
  f12 <- frm(y ~ x, data = data.frame(x = x, y = s + 1),
             family = bernoulli())
  f01 <- frm(y ~ x, data = data.frame(x = x, y = s), family = bernoulli())
  # residuals(newdata = ), as brms, read 1 as the value 1, a failure
  r12 <- residuals(f12, newdata = data.frame(x = x, y = 1),
                   type = "response")[, 1]
  r01 <- residuals(f01, newdata = data.frame(x = x, y = 0),
                   type = "response")[, 1]
  expect_identical(r12, r01)
  # the values, then the codes, then a marked vector of either
  expect_identical(fixef(refit(f12, s + 1)), fixef(refit(f01, s)))
  expect_identical(fixef(refit(f12, s)), fixef(refit(f01, s)))
  expect_identical(fixef(refit(f12, factor(s + 1, levels = 1:2))),
                   fixef(refit(f01, s)))
  expect_identical(fixef(refit(f12, as.logical(s))), fixef(refit(f01, s)))
  # the frmtmb 0.82.0 lane build read rep(1, n) as codes, a success, and
  # c(2, 1, 1, ...) as values, so 1 was a failure there
  expect_error(refit(f12, rep(1, n)), "reads two ways on this fit",
               fixed = TRUE, class = "frmtmb_error")
  # the bootstrap's draws are codes, whatever the fit's values are
  b12 <- frm_bootstrap(f12, FUN = function(f) fixef(f)[, 1], nsim = 6,
                       seed = 5)
  b01 <- frm_bootstrap(f01, FUN = function(f) fixef(f)[, 1], nsim = 6,
                       seed = 5)
  expect_identical(b12$t, b01$t)
})

test_that("update() records a delta as a formula a reader can read", {
  set.seed(7)
  d <- data.frame(x = stats::rnorm(60), x3 = stats::rnorm(60),
                  g = factor(rep(1:10, 6)))
  d$z <- d$x + stats::rnorm(10)[d$g] + stats::rnorm(60)
  fk <- frm(z ~ x + x3 + (1 | g), data = d)
  u <- update(fk, ~ . - x3)
  # 0.82.0 recorded the deparsed bf() object, 3125 characters here
  expect_identical(deparse1(getCall(u)$formula), "z ~ x + (1 | g)")
  expect_identical(fixef(eval(getCall(u))), fixef(u))
  fs <- frm(bf(z ~ x + x3, sigma ~ x), data = d)
  us <- update(fs, ~ . - x3)
  expect_identical(getCall(us)$formula[[1L]], quote(frmtmb::bf))
  expect_identical(fixef(eval(getCall(us))), fixef(us))
  # a family given to frm() stays in its argument
  fp <- frm(z ~ x + x3, family = student(), data = d)
  up <- update(fp, ~ . - x3)
  expect_identical(deparse1(getCall(up)$formula), "z ~ x")
  expect_identical(getCall(up)$family, quote(student()))
})

test_that("confint() of a saved ar(cov = TRUE) fit answers on the mean", {
  set.seed(4)
  d <- data.frame(t = rep(1:12, 5), h = factor(rep(1:5, each = 12)),
                  x = stats::rnorm(60))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(60)
  fit <- frm(y ~ x + ar(time = t, gr = h, cov = TRUE), data = d)
  old <- fit
  # what frmtmb 0.76.0 stored
  old$frame$autocor[[1L]]$sigma_is <- NULL
  # 0.82.0 refused this one too
  expect_identical(confint(old, parm = "x"), confint(fit, parm = "x"))
  msg <- "frmtmb 0.76.0 or earlier"
  expect_error(confint(old), msg, class = "frmtmb_error")
  expect_error(confint(old, parm = "sigma_(Intercept)"), msg,
               class = "frmtmb_error")
})

test_that("vcov_cluster() reads an outside cluster only from the fit's data", {
  set.seed(1)
  dd <- data.frame(cl = factor(rep(1:10, each = 12)), x = stats::rnorm(120))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(120)
  fit <- frm(y ~ x, data = dd)
  ref <- vcov_cluster(fit, dd$cl)
  expect_identical(vcov_cluster(fit, ~ cl), ref)
  # rows reordered since the fit, row names kept: 0.82.0 read the
  # cluster by position, so every row took its neighbor's
  d0 <- dd
  dd <- dd[c(2:120, 1L), ]
  expect_identical(vcov_cluster(fit, ~ cl), ref)
  # rows dropped since the fit
  dd <- d0[-1L, ]
  expect_error(vcov_cluster(fit, ~ cl), "no longer hold the data",
               fixed = TRUE, class = "frmtmb_error")
  rm(dd)
  expect_error(vcov_cluster(fit, ~ cl), "which are gone", fixed = TRUE,
               class = "frmtmb_error")
})

test_that("vcov_cluster() reads the data's column before the workspace's", {
  set.seed(1)
  dd <- data.frame(cl = factor(rep(1:10, each = 12)), x = stats::rnorm(120))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(120)
  fit <- frm(y ~ x, data = dd)
  # an object of the same name and length where the formula is written:
  # 0.82.0 read it in place of dd$cl, as sandwich::vcovCL() does not
  cl <- factor(rep(1:2, length.out = 120))
  expect_identical(vcov_cluster(fit, ~ cl), vcov_cluster(fit, dd$cl))
  # data that are gone may have held cl: the workspace's is not it
  # (0.82.0 and the first fix round read it)
  mk <- function() {
    d_local <- dd
    frm(y ~ x, data = d_local)
  }
  fl <- mk()
  expect_error(vcov_cluster(fl, ~ cl), "which are gone", fixed = TRUE,
               class = "frmtmb_error")
})

test_that("vcov_cluster() on a refit of a matrix response reads its data", {
  set.seed(5)
  n <- 200
  x <- stats::rnorm(n)
  p <- cbind(a = 1, b = exp(0.3 + 0.4 * x), c = exp(-0.2 * x))
  p <- p / rowSums(p)
  Y <- t(vapply(seq_len(n), function(i) stats::rmultinom(1, 12, p[i, ]),
                numeric(3)))
  colnames(Y) <- c("a", "b", "c")
  dd <- data.frame(x = x, N = 12L, cl = factor(rep(1:20, each = 10)))
  dd$Y <- Y
  fm <- frm(Y | trials(N) ~ x, data = dd, family = multinomial())
  set.seed(6)
  rf <- refit(fm, Y[sample(n), ])
  # the refit's response is not its call's data's; the first fix round
  # refused the call's data for that reason
  expect_identical(vcov_cluster(rf, ~ cl), vcov_cluster(rf, dd$cl))
})

test_that("influence() of a refit deletes rows of the refit's response", {
  set.seed(8)
  d <- data.frame(x = stats::rnorm(30), y = stats::rnorm(30))
  fu <- frm(y ~ x, data = d)
  ynew <- stats::rnorm(30, 2)
  ru <- refit(fu, ynew)
  d2 <- d
  d2$y <- ynew
  iu <- influence(ru, force = TRUE)$fixed
  i2 <- influence(frm(y ~ x, data = d2), force = TRUE)$fixed
  i0 <- influence(fu, force = TRUE)$fixed
  # 0.82.0 deleted rows of the original response: the table was the
  # original fit's, 1.85 away from the refit's
  expect_lt(max(abs(iu - i2)) / max(abs(i2 - i0)), 1e-4)
})

test_that("update() takes na.action and se from the fit", {
  set.seed(2)
  d <- data.frame(x = stats::rnorm(50))
  d$y <- 1 + d$x + stats::rnorm(50)
  d$x[3] <- NA
  mk <- function(dat) {
    my_na <- stats::na.exclude
    s <- FALSE
    frm(y ~ x, data = dat, na.action = my_na, se = s)
  }
  fit <- suppressMessages(mk(d))
  # 0.82.0: "object 'my_na' not found", and an `s` in the caller's
  # workspace became the se flag
  s <- "not a flag"
  u <- update(fit)
  expect_identical(fixef(u), fixef(fit))
  expect_identical(NROW(residuals(u)), 50L)
  u2 <- suppressMessages(update(fit, newdata = d))
  expect_identical(NROW(residuals(u2)), 50L)
})

test_that("marginaleffects predicts a fit with a matrix response", {
  skip_if_not_installed("marginaleffects")
  set.seed(5)
  d <- data.frame(x = stats::rnorm(100), N = 10L)
  pr <- cbind(1, exp(0.3 + 0.5 * d$x), exp(-0.2 * d$x))
  pr <- pr / rowSums(pr)
  d$Y <- t(sapply(1:100, function(i) stats::rmultinom(1, 10, pr[i, ])))
  colnames(d$Y) <- c("a", "b", "c")
  fm <- frm(Y | trials(N) ~ x, family = multinomial(), data = d)
  # 0.82.0: "Column 3 ['Y'] is length 300 but column 1 is length 100"
  p <- marginaleffects::predictions(fm)
  ft <- fitted(fm)
  expect_identical(nrow(p), 300L)
  expect_equal(p$estimate[p$group == "b"], unname(ft[, 1, 2]),
               tolerance = 1e-10)
  # a refit's data hold its own response
  Y2 <- d$Y[c(2:100, 1L), ]
  r2 <- refit(fm, newresp = Y2)
  got <- as.matrix(attr(r2, "marginaleffects_modeldata")[, c("Y.a", "Y.b",
                                                               "Y.c")])
  expect_true(all(unname(got) == unname(Y2)))
  # insight still reads the matrix response
  if (requireNamespace("insight", quietly = TRUE)) {
    expect_identical(dim(insight::get_response(fm)), c(100L, 3L))
  }
  # the absent case: a fit without a matrix column carries nothing
  dv <- data.frame(x = d$x, y = d$x + stats::rnorm(100))
  expect_null(attr(frm(y ~ x, data = dv), "marginaleffects_modeldata"))
})

test_that("a recover_data() refusal keeps its class when data are given", {
  skip_if_not_installed("emmeans")
  set.seed(6)
  d <- data.frame(x = stats::rnorm(50), f = factor(rep(c("a", "b"), 25)))
  d$y <- 1 + d$x + stats::rnorm(50)
  fit <- frm(y ~ x + f, data = d)
  # 0.82.0 raised a simpleError on both routes
  expect_error(emmeans::ref_grid(fit, dpar = "nope", data = d),
               "'nope' is not part of the model", class = "frmtmb_error")
  # without data, emmeans's try() leaves the text alone
  expect_error(emmeans::ref_grid(fit, dpar = "nope"),
               "'nope' is not part of the model")
})

test_that("an influence() deletion keeps the fit's smooth basis", {
  set.seed(8)
  d <- data.frame(x = stats::runif(60, 0, 3), z = stats::rnorm(60),
                  blk = factor(rep(1:5, 12)))
  d$y <- sin(d$x) + 0.3 * d$z + stats::rnorm(60, 0, 0.3)
  fl <- frm(y ~ z + s(x, k = 5), data = d)
  inf <- influence(fl, groups = "blk", data = d)
  # the deletion of block 1 with the fit's basis written out: the wiggly
  # columns at the kept rows as one iid block with one sd, the
  # null-space column as a fixed effect
  si <- fl$frame$linpreds[[1L]]$smooths[[1L]]
  keep <- d[d$blk != "1", ]
  M <- frmtmb:::smooth_basis_at(si, keep)
  nr <- sum(si$nr)
  for (k in seq_len(nr)) keep[[paste0("r", k)]] <- M[, k]
  keep$fx <- M[, nr + 1L]
  keep$id <- factor(1)
  fo <- stats::as.formula(paste0("y ~ z + fx + homdiag(0 + ",
                                 paste0("r", seq_len(nr), collapse = " + "),
                                 " | id)"))
  ref <- fixef(frm(fo, data = keep))[, 1]
  got <- inf$fixed["1", c("z", "s(x).fx1")]
  moved <- abs(got - inf$fixed_full[c("z", "s(x).fx1")])
  # the deletion moves each coefficient; the two fits agree far closer
  # than that. 0.82.0 rebuilt the basis on the kept rows
  expect_lt(max(abs(got - ref[c("z", "fx")]) / moved), 1e-3)
})

test_that("an uncorrelated factor over gr(cov = ) is brms's model", {
  set.seed(11)
  ng <- 12
  A <- 0.5 ^ abs(outer(1:ng, 1:ng, "-"))
  dimnames(A) <- list(paste0("g", 1:ng), paste0("g", 1:ng))
  d <- data.frame(g = factor(rep(paste0("g", 1:ng), each = 10),
                             levels = paste0("g", 1:ng)),
                  f = factor(rep(c("a", "b", "c"), 40)))
  L <- t(chol(A))
  r <- cbind(0.8 * L %*% stats::rnorm(ng), 0.6 * L %*% stats::rnorm(ng),
             1.2 * L %*% stats::rnorm(ng))
  gi <- as.integer(d$g)
  d$y <- 1 + r[cbind(gi, as.integer(d$f))] + stats::rnorm(120, 0, 0.5)
  d$fb <- as.numeric(d$f == "b")
  d$fc <- as.numeric(d$f == "c")
  d$fa <- as.numeric(d$f == "a")
  # 0.82.0: "it is not available yet"
  f1 <- frm(y ~ 1 + (0 + f || gr(g, cov = A)), data = d,
            data2 = list(A = A))
  f2 <- frm(y ~ 1 + (0 + fa + fb + fc || gr(g, cov = A)), data = d,
            data2 = list(A = A))
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-10)
  expect_identical(grep("^sd_", variables(f1), value = TRUE),
                   c("sd_g__fa", "sd_g__fb", "sd_g__fc"))
  # the gaussian marginal density at the fit's estimates, written out:
  # y ~ N(b0, sum_k sd_k^2 Z_k A Z_k' + sigma^2 I)
  vc <- VarCorr(f1)
  sds <- vc$g$sd[, "Estimate"]
  V <- diag(vc$residual__$sd[, "Estimate"]^2, 120)
  for (k in 1:3) {
    Z <- matrix(0, 120, ng)
    Z[cbind(1:120, gi)] <- as.integer(d$f) == k
    V <- V + sds[k]^2 * Z %*% A %*% t(Z)
  }
  R <- chol(V)
  e <- backsolve(R, d$y - fixef(f1)[1, 1], transpose = TRUE)
  exact <- -sum(log(diag(R))) - 0.5 * sum(e^2) - 60 * log(2 * pi)
  expect_equal(as.numeric(logLik(f1)), exact, tolerance = 1e-9)
  # an intercept beside the factor: brms codes Intercept, fb, fc
  h1 <- frm(y ~ 1 + (1 + f || gr(g, cov = A)), data = d,
            data2 = list(A = A))
  h2 <- frm(y ~ 1 + (1 + fb + fc || gr(g, cov = A)), data = d,
            data2 = list(A = A))
  expect_equal(as.numeric(logLik(h1)), as.numeric(logLik(h2)),
               tolerance = 1e-10)
  expect_identical(grep("^sd_", variables(h1), value = TRUE),
                   c("sd_g__Intercept", "sd_g__fb", "sd_g__fc"))
  # new data rebuild each column's design
  nd <- data.frame(g = factor(c("g1", "g5"), levels = levels(d$g)),
                   f = factor(c("b", "c"), levels = c("a", "b", "c")),
                   fa = 0, fb = c(1, 0), fc = c(0, 1))
  expect_equal(fitted(f1, newdata = nd)[, 1], fitted(f2, newdata = nd)[, 1],
               tolerance = 1e-10)
  # and a known precision is the same model
  g1 <- frm(y ~ 1 + (0 + f || gr(g, prec = Q)), data = d,
            data2 = list(Q = solve(A)))
  expect_equal(as.numeric(logLik(g1)), as.numeric(logLik(f1)),
               tolerance = 1e-8)
})
