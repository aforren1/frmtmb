# update() with a complete formula keeps the stored parameter formulas,
# as brms's update.brmsfit() does through update.brmsformula(), which
# pools object$pforms with the new ones. Before, a complete formula
# replaced the whole bf(), and bf(y ~ a + b, nl = TRUE) was refused at
# bf() for having no parameter formulas (ledger row brmsfit-methods:955;
# brms's own output in dev/formrobust-log/brms-update.txt).

up_data <- local({
  set.seed(41)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- 2 + 1.5 * d$x + stats::rnorm(n, 0, 0.5)
  d
})

test_that("a new nonlinear body keeps the fit's parameter formulas", {
  f0 <- frm(bf(y ~ a + b * x, a ~ 1 + z, b ~ 1, nl = TRUE), data = up_data)
  up <- NULL
  expect_message(
    up <- update(f0, formula. = bf(y ~ b * x + a, nl = TRUE)),
    "will completely replace the original formula in non-linear models",
    fixed = TRUE)
  # the same model written with the terms swapped: the same fit
  expect_identical(sort(rownames(fixef(up))), sort(rownames(fixef(f0))))
  expect_lt(abs(as.numeric(logLik(up) - logLik(f0))),
            1e-10 * abs(as.numeric(logLik(f0))))
  expect_identical(deparse1(up$bform$pforms[["a"]]), "a ~ 1 + z")
  # a formula the update gives replaces the stored one of that name
  msgs <- character()
  up2 <- withCallingHandlers(
    update(f0, formula. = bf(y ~ a + b * x, a ~ 1, nl = TRUE)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("Replacing initial definitions of parameters a",
                        msgs, fixed = TRUE)))
  expect_identical(sort(rownames(fixef(up2))), c("a_Intercept",
                                                 "b_Intercept"))
})

test_that("a plain formula keeps a linear model's dpar formulas", {
  f0 <- frm(bf(y ~ x, sigma ~ z), data = up_data)
  up <- update(f0, y ~ x + z)
  # the call holds the pooled formula as the bf() call that builds it,
  # not the evaluated object with its family's closures
  expect_identical(deparse1(getCall(up)$formula),
                   paste0("frmtmb::bf(y ~ x + z, sigma ~ z, ",
                          "family = stats::gaussian(link = \"identity\"))"))
  expect_identical(logLik(eval(getCall(up))), logLik(up))
  expect_true(all(c("sigma_Intercept", "sigma_z", "z") %in%
                    rownames(fixef(up))))
  # the guard-absent case: a model with nothing to keep is refitted from
  # the formula as written, and its call keeps that formula
  f1 <- frm(bf(y ~ x), data = up_data)
  u1 <- update(f1, y ~ x + z)
  expect_identical(deparse1(getCall(u1)$formula), "y ~ x + z")
})

test_that("bf(nl = TRUE) without parameter formulas waits for frm()", {
  f <- bf(y ~ a + b * x, nl = TRUE)
  expect_s3_class(f, "frmtmb_formula")
  expect_error(frm(f, data = up_data),
               "needs at least one nonlinear-parameter formula",
               fixed = TRUE)
  # lf() added afterwards completes it, as in brms
  fit <- frm(bf(y ~ a + b * x, nl = TRUE) + lf(a ~ 1) + lf(b ~ 1),
             data = up_data)
  expect_identical(rownames(fixef(fit)), c("a_Intercept", "b_Intercept"))
})

test_that("a complete formula on a multivariate fit is refused, as brms", {
  fm <- frm(bf(y ~ x, sigma ~ z) + bf(z ~ x), data = up_data)
  expect_error(update(fm, y ~ x + z),
               "Updating formulas of multivariate models is not yet",
               fixed = TRUE)
  expect_error(update(fm, mvbf(bf(y ~ x), bf(z ~ x))),
               "Updating formulas of multivariate models is not yet",
               fixed = TRUE)
  # the guard-absent case: the data can still be updated
  up <- update(fm, newdata = up_data[1:60, ])
  expect_identical(nobs(up), 60L)
})

test_that("a pooled update keeps a family's options that are not links", {
  # the round-1 call wrote huber(k = 3) as huber(link = , link_sigma = )
  # and the update refitted at the default k = 1.345 without a word
  # (dev/formrobust-rev2-famcall.R, logLik -231.496 against -237.750)
  set.seed(95)
  d <- data.frame(x = stats::rnorm(120), z = stats::rnorm(120))
  d$y <- 0.5 * d$x + stats::rt(120, 3)
  f0 <- frm(bf(y ~ x, sigma ~ z, family = huber(k = 3)), data = d)
  up <- update(f0, y ~ x + z)
  ref <- frm(bf(y ~ x + z, sigma ~ z, family = huber(k = 3)), data = d)
  expect_identical(logLik(up), logLik(ref))
  # the family a constructor call cannot say is stored as the object
  for (fam in list(huber(k = 3), whittle(tapers = 4), cox(df = 6))) {
    expect_s3_class(family_call_of(as_frmtmb_family(fam)), "frmtmb_family")
  }
  # the guard-absent case: a constructor of links alone is written as
  # its call, which a reader can read
  cl <- family_call_of(as_frmtmb_family(student(link_sigma = "identity")))
  expect_true(is.call(cl))
  expect_identical(deparse1(cl), paste0("frmtmb::student(link = ",
                                        "\"identity\", link_sigma = ",
                                        "\"identity\", link_nu = \"logm1\")"))
})
