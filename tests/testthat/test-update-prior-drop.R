# update() drops a stored prior the updated model cannot take, as
# brms's update() does (lane surface, 2026-10-06, vigport defect 4).
#
# Seen to fail on 0.68.1 (rellib-r6): the update of brms_overview's
# fit2, which removes the last correlation, stopped with "No
# random-effect correlations match class=cor" (dev/surface-repros.R).
# brms 2.23.0's update.brmsfit() marks the stored prior
# `allow_invalid_prior` and drops the lkj(2) row silently
# (dev/vigport-rev-brms-update.R, dev/surface-brms-src2.R).

kidney_like <- function() {
  set.seed(3)
  n_p <- 30
  d <- data.frame(patient = factor(rep(seq_len(n_p), each = 2)),
                  age = rep(stats::runif(n_p, 20, 70), each = 2),
                  sex = factor(rep(sample(c("m", "f"), n_p, TRUE),
                                   each = 2)))
  u <- stats::rnorm(n_p, 0, 0.5)[d$patient]
  d$time <- exp(3 + 0.01 * d$age + u + stats::rnorm(2 * n_p, 0, 0.6))
  d
}

pr3 <- function() {
  c(set_prior("normal(0,5)", class = "b"),
    set_prior("cauchy(0,2)", class = "sd"),
    set_prior("lkj(2)", class = "cor"))
}

test_that("update() drops a cor prior once the last correlation is gone", {
  d <- kidney_like()
  f1 <- allow_warnings(
    frm(time ~ age + sex + (1 + age | patient), data = d,
        family = lognormal(), prior = pr3()),
    c("Standard errors are not available", "singular"))
  msg <- character()
  u <- withCallingHandlers(
    allow_warnings(update(f1, formula. = ~ . - (1 + age | patient) +
                            (1 | patient)),
                   c("Standard errors are not available", "singular")),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_s3_class(u, "frmtmb_fit")
  # said, where brms is silent
  expect_true(any(grepl("update() drops the prior specification",
                        msg, fixed = TRUE)))
  expect_true(any(grepl("lkj(2) on cor", msg, fixed = TRUE)))
  # the two that still apply are kept, and only those
  kept <- as.data.frame(u$prior)
  expect_setequal(kept$class, c("b", "sd"))
  direct <- allow_warnings(
    frm(time ~ age + sex + (1 | patient), data = d, family = lognormal(),
        prior = c(set_prior("normal(0,5)", class = "b"),
                  set_prior("cauchy(0,2)", class = "sd"))),
    c("Standard errors are not available", "singular"))
  expect_identical(as.numeric(logLik(u)), as.numeric(logLik(direct)))
})

test_that("a prior that still applies is kept without a message", {
  # the absent case: nothing to drop, nothing said
  d <- kidney_like()
  f1 <- allow_warnings(
    frm(time ~ age + sex + (1 + age | patient), data = d,
        family = lognormal(), prior = pr3()),
    c("Standard errors are not available", "singular"))
  expect_silent(u <- allow_warnings(update(f1, formula. = ~ . - sex),
                                    c("Standard errors are not available",
                                      "singular")))
  expect_equal(nrow(as.data.frame(u$prior)), 3L)
})

test_that("a prior given to update() itself is still checked", {
  d <- kidney_like()
  f1 <- allow_warnings(
    frm(time ~ age + sex + (1 + age | patient), data = d,
        family = lognormal(), prior = pr3()),
    c("Standard errors are not available", "singular"))
  expect_error(update(f1, formula. = ~ . - (1 + age | patient) +
                        (1 | patient), prior = pr3()),
               "No random-effect correlations match")
  # and frm() refuses the same prior written directly, as brms does
  expect_error(frm(time ~ age + (1 | patient), data = d,
                   family = lognormal(),
                   prior = set_prior("lkj(2)", class = "cor")),
               "No random-effect correlations match")
})
