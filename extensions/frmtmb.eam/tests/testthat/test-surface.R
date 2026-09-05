## The rest of the frmtmb surface on a fitted wiener model. A family is
## not finished when it fits; it is finished when the methods a user
## reaches for next either work or refuse for a stated reason.

ddm_fit <- function(n = 400, seed = 404) {
  set.seed(seed)
  cond <- rep(c(0, 1), each = n / 2)
  dat <- ddm_simulate(n, mu = 0.4 + 0.9 * cond, bs = 1.4, ndt = 0.28)
  dat$cond <- cond
  list(dat = dat,
       fit = frm(bf(rt | vint(upper) ~ cond, bias = 0.5),
                 family = wiener(), data = dat))
}

test_that("summary, fixef and logLik work", {
  skip_if_not_installed("RWiener")
  f <- ddm_fit()$fit
  expect_s3_class(f, "frmtmb_fit")
  expect_no_error(summary(f))
  fx <- fixef(f)
  expect_setequal(names(fx), c("mu", "bs", "ndt", "bias"))
  expect_true(is.finite(as.numeric(logLik(f))))
  expect_true(is.finite(AIC(f)))
})

test_that("fitted and predict return conditional mean response times", {
  skip_if_not_installed("RWiener")
  o <- ddm_fit()
  f <- o$fit
  ft <- fitted(f)
  expect_length(ft, nrow(o$dat))
  expect_true(all(is.finite(ft)))
  # a mean response time is a non-decision time plus a positive
  # decision time, so it exceeds ndt everywhere
  e <- unlist(fixef(f))
  ndt_hat <- min(o$dat$rt) / (1 + exp(-e[["ndt.(Intercept)"]]))
  expect_true(all(ft > ndt_hat))
  # and it is in the right ballpark for the data it was fitted to
  expect_equal(mean(ft), mean(o$dat$rt), tolerance = 0.1)

  expect_no_error(predict(f, type = "link"))
  expect_equal(predict(f, type = "response"), ft, tolerance = 1e-8)
})

test_that("predict on newdata requires the decision indicator", {
  skip_if_not_installed("RWiener")
  f <- ddm_fit()$fit
  nd <- data.frame(cond = c(0, 1), upper = c(1, 1))
  p <- predict(f, newdata = nd, type = "response")
  expect_length(p, 2L)
  expect_true(all(is.finite(p)))
  # omitting it is refused by frmtmb, not silently defaulted
  expect_error(predict(f, newdata = data.frame(cond = c(0, 1)),
                       type = "response"),
               "could not be evaluated on newdata")

  # At bias = 0.5, as above, the two boundaries share a mean, so the
  # boundary must be varied on a fit with a FREE bias to show that it
  # reaches the prediction at all.
  set.seed(99)
  d2 <- ddm_simulate(900, mu = 0.7, bs = 1.4, ndt = 0.25, bias = 0.35)
  f2 <- frm(bf(rt | vint(upper) ~ 1), family = wiener(), data = d2)
  p1 <- predict(f2, newdata = data.frame(upper = 1), type = "response")
  p0 <- predict(f2, newdata = data.frame(upper = 0), type = "response")
  expect_gt(abs(p1 - p0), 0.05)
  # and each tracks the observed mean at its own boundary
  expect_equal(p1, mean(d2$rt[d2$upper == 1]), tolerance = 0.05)
  expect_equal(p0, mean(d2$rt[d2$upper == 0]), tolerance = 0.05)
})

test_that("all four dpars recover when none is held fixed", {
  skip_if_not_installed("RWiener")
  set.seed(99)
  d <- ddm_simulate(900, mu = 0.7, bs = 1.4, ndt = 0.25, bias = 0.35)
  f <- frm(bf(rt | vint(upper) ~ 1), family = wiener(), data = d)
  e <- unlist(fixef(f))
  expect_equal(e[["mu.(Intercept)"]], 0.7, tolerance = 0.2)
  expect_equal(exp(e[["bs.(Intercept)"]]), 1.4, tolerance = 0.1)
  expect_equal(1 / (1 + exp(-e[["bias.(Intercept)"]])), 0.35,
               tolerance = 0.05)
  expect_equal(min(d$rt) / (1 + exp(-e[["ndt.(Intercept)"]])), 0.25,
               tolerance = 0.03)
})

test_that("residuals: response works, pearson and deviance refuse", {
  skip_if_not_installed("RWiener")
  o <- ddm_fit()
  r <- residuals(o$fit, type = "response")
  expect_length(r, nrow(o$dat))
  expect_equal(r, o$dat$rt - fitted(o$fit), tolerance = 1e-8)
  # declared omissions, and the refusals name the reason
  expect_error(residuals(o$fit, type = "pearson"), "no variance function")
  expect_error(residuals(o$fit, type = "deviance"), "unit deviance")
})

test_that("simulate draws response times at each row's own boundary", {
  skip_if_not_installed("RWiener")
  o <- ddm_fit()
  set.seed(11)
  s <- as.matrix(simulate(o$fit, nsim = 20))
  expect_equal(dim(s), c(nrow(o$dat), 20L))
  expect_true(all(is.finite(s)))
  e <- unlist(fixef(o$fit))
  ndt_hat <- min(o$dat$rt) / (1 + exp(-e[["ndt.(Intercept)"]]))
  expect_true(all(s > ndt_hat))
  # the draws sit around the observed response times
  expect_equal(mean(s), mean(o$dat$rt), tolerance = 0.08)
})

test_that("par_template and set_prior reach the wiener dpars", {
  skip_if_not_installed("RWiener")
  o <- ddm_fit()
  tmpl <- par_template(o$fit)
  expect_true(is.list(tmpl))
  # the dpar intercepts live in betad, under their dpar-prefixed names
  nms <- names(unlist(tmpl))
  expect_true(any(grepl("ndt_", nms, fixed = TRUE)))
  expect_true(any(grepl("bs_", nms, fixed = TRUE)))
  expect_true(any(grepl("bias_", nms, fixed = TRUE)))

  gp <- get_prior(bf(rt | vint(upper) ~ cond, bias = 0.5),
                  family = wiener(), data = o$dat)
  expect_true(is.data.frame(gp))
  # mu is the main formula, so its rows carry an empty dpar; the other
  # dpars are named. bias is fixed at 0.5, so it has no targetable slot.
  expect_true(all(c("bs", "ndt") %in% gp$dpar))
  expect_false("bias" %in% gp$dpar)
  expect_true(any(gp$class == "b" & gp$coef == "cond"))

  # a prior on a dpar changes the answer it is set on
  f2 <- frm(bf(rt | vint(upper) ~ cond, bias = 0.5), family = wiener(),
            data = o$dat,
            prior = set_prior("normal(0, 0.1)", class = "b", dpar = "mu"))
  # toward zero, and not past it. Bounding only the magnitude admits a
  # sign flip, which is what a natural-scale placement for a dpar prior
  # produces here: the density on exp(coef) pushes the coefficient
  # negative and |shrunk| < |unpenalized| stays true while the estimate
  # changes sign. The ratio bounds both at once.
  shrink <- unlist(fixef(f2))[["mu.cond"]] /
    unlist(fixef(o$fit))[["mu.cond"]]
  expect_gt(shrink, 0)
  expect_lt(shrink, 1)
})

test_that("weights() and a random effect both work", {
  skip_if_not_installed("RWiener")
  o <- ddm_fit()
  d <- o$dat
  d$wt <- 1
  fw <- frm(bf(rt | vint(upper) + weights(wt) ~ cond, bias = 0.5),
            family = wiener(), data = d)
  expect_equal(as.numeric(logLik(fw)), as.numeric(logLik(o$fit)),
               tolerance = 1e-6)

  d$subj <- factor(rep(1:20, length.out = nrow(d)))
  fr <- frm(bf(rt | vint(upper) ~ cond + (1 | subj), bias = 0.5),
            family = wiener(), data = d)
  expect_s3_class(fr, "frmtmb_fit")
  expect_true(is.finite(as.numeric(logLik(fr))))
  expect_true(nrow(VarCorr(fr)$subj$sd) >= 1 ||
                length(ranef(fr)) >= 1)
})

test_that("cens() and trunc() are refused for want of a CDF", {
  skip_if_not_installed("RWiener")
  o <- ddm_fit()
  d <- o$dat
  d$cc <- 0
  expect_error(frm(bf(rt | vint(upper) + trunc(lb = 0.3) ~ cond,
                      bias = 0.5), family = wiener(), data = d),
               "need a family with a CDF")
  expect_error(frm(bf(rt | vint(upper) + cens(cc) ~ cond, bias = 0.5),
                   family = wiener(), data = d),
               "need a family with a CDF")
})

test_that("dec() is the spelling now, and vint() is the same model", {
  skip_if_not_installed("RWiener")
  # This assertion used to read the other way round: it pinned core's
  # refusal of `dec()`, because frmtmb's addition terms were a closed
  # set and the only route for the indicator was vint(). frmtmb 0.49.0
  # added frmtmb_register_aterm(), this package registers `dec` when it
  # loads, and the refusal it pinned no longer exists. The pin moves to
  # the property that replaced it: the two spellings are one model.
  o <- ddm_fit()
  f_dec <- frm(bf(rt | dec(upper) ~ cond, bias = 0.5),
               family = wiener(), data = o$dat)
  expect_equal(as.numeric(logLik(f_dec)), as.numeric(logLik(o$fit)),
               tolerance = 1e-10)
  expect_equal(unlist(fixef(f_dec)), unlist(fixef(o$fit)),
               tolerance = 1e-6)
})
