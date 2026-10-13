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

# update(prior = ) merges the new prior into the stored one, as brms's
# update.brmsfit() does (dev/moparity-findings.md, item 3). Seen to
# fail on 0.69.0 (rellib-r7): the new prior replaced the stored one, so
# a stored normal(0, 0.05) on b_x was lost (b_x 0.236 to 0.766 in the
# surface review's construction). brms 2.23.0 keeps it
# (dev/moparity-brms-update.R: stored rows normal(0, 0.05) on b x and
# normal(0, 10) on Intercept, then exponential(2) on sigma added).
merge_data <- function() {
  set.seed(11)
  dd <- data.frame(x = rnorm(60), z = rnorm(60))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
  dd
}
merge_pr0 <- function() {
  c(set_prior("normal(0, 0.05)", class = "b", coef = "x"),
    set_prior("normal(0, 10)", class = "Intercept"))
}

test_that("update(prior = ) keeps the stored rows of other slots", {
  dd <- merge_data()
  f0 <- frm(y ~ x + z, family = gaussian(), data = dd, prior = merge_pr0())
  sig <- set_prior("exponential(2)", class = "sigma")
  u <- update(f0, prior = sig)
  direct <- frm(y ~ x + z, family = gaussian(), data = dd,
                prior = c(sig, merge_pr0()))
  expect_identical(fixef(u), fixef(direct))
  expect_setequal(as.data.frame(u$prior)$class,
                  c("b", "Intercept", "sigma"))
  # the inverse: replacing the prior, as 0.69.0 did, moves b_x
  replaced <- frm(y ~ x + z, family = gaussian(), data = dd, prior = sig)
  expect_gt(abs(fixef(replaced)["x", 1]), 5 * abs(fixef(u)["x", 1]))
})

test_that("a new specification wins its slot, and NULL keeps the prior", {
  dd <- merge_data()
  f0 <- frm(y ~ x + z, family = gaussian(), data = dd, prior = merge_pr0())
  bx <- set_prior("normal(0, 1)", class = "b", coef = "x")
  u <- update(f0, prior = bx)
  direct <- frm(y ~ x + z, family = gaussian(), data = dd,
                prior = c(bx, set_prior("normal(0, 10)",
                                        class = "Intercept")))
  expect_identical(fixef(u), fixef(direct))
  expect_equal(nrow(as.data.frame(u$prior)), 2L)
  u0 <- update(f0, prior = NULL)
  expect_identical(fixef(u0), fixef(f0))
})

test_that("a stored row is dropped and a new row is checked", {
  dd <- merge_data()
  f0 <- frm(y ~ x + z, family = gaussian(), data = dd, prior = merge_pr0())
  sig <- set_prior("exponential(2)", class = "sigma")
  msg <- character()
  u <- withCallingHandlers(update(f0, formula. = ~ . - x, prior = sig),
                           message = function(m) {
                             msg <<- c(msg, conditionMessage(m))
                             invokeRestart("muffleMessage")
                           })
  expect_true(any(grepl("normal(0, 0.05) on b_x", msg, fixed = TRUE)))
  expect_setequal(as.data.frame(u$prior)$class, c("Intercept", "sigma"))
  # brms drops a new row that matches nothing; frmtmb refuses it, as
  # frm() refuses it written directly
  expect_error(update(f0, prior = set_prior("normal(0, 1)", class = "b",
                                            coef = "nosuch")),
               "nosuch")
})

# Every spelling frm() takes merges as set_prior() does (punch round 1,
# B1). Seen to fail on the lane build before the fix: a brmsprior and a
# default_prior() table went to frm() untranslated and replaced the
# stored prior, so b_x moved from 0.081 to 0.456 in the review's
# construction (dev/moparity-rev-update.R, dev/moparity-b1-probe.R).
test_that("update(prior = ) merges a brms prior() object", {
  skip_if_not_installed("brms")
  dd <- merge_data()
  f0 <- frm(y ~ x + z, family = gaussian(), data = dd, prior = merge_pr0())
  bp <- brms::prior(exponential(2), class = sigma)
  u <- update(f0, prior = bp)
  direct <- frm(y ~ x + z, family = gaussian(), data = dd,
                prior = c(set_prior("exponential(2)", class = "sigma"),
                          merge_pr0()))
  expect_identical(fixef(u), fixef(direct))
  expect_equal(nrow(as.data.frame(u$prior)), 3L)
  # a fit made with a brmsprior merges the same way
  fb <- frm(y ~ x + z, family = gaussian(), data = dd,
            prior = brms::prior(normal(0, 0.05), class = b, coef = x) +
              brms::prior(normal(0, 10), class = Intercept))
  ub <- update(fb, prior = bp)
  expect_identical(fixef(ub), fixef(direct))
  # an empty brms row takes its slot, as in brms
  ue <- update(f0, prior = brms::prior_string("", class = "b", coef = "x"))
  flat_x <- frm(y ~ x + z, family = gaussian(), data = dd,
                prior = set_prior("normal(0, 10)", class = "Intercept"))
  expect_identical(fixef(ue), fixef(flat_x))
})

# Every row of a prior table takes its slot, an empty row too (punch
# round 2, B2). brms 2.23.0 binds the new rows ahead of the stored ones
# and drops the stored row of each slot a new row names, so "" returns
# the slot to its default: a one-row table b x "" gives b_x 0.625 in
# brms, 0.0511 on the lane before the fix, which kept the stored
# normal(0, 0.05) (dev/moparity-rev2-brms-update.R,
# dev/moparity-rev2-update.R). Seen to fail on the lane build of round 1.
test_that("a prior table's rows take their slots, empty rows too", {
  dd <- merge_data()
  f0 <- frm(y ~ x + z, family = gaussian(), data = dd, prior = merge_pr0())
  tbl <- default_prior(y ~ x + z, data = dd, family = gaussian())
  # one edited row merges as set_prior() does
  sig <- tbl[tbl$class == "sigma", ]
  sig$prior <- "exponential(2)"
  u <- update(f0, prior = sig)
  direct <- frm(y ~ x + z, family = gaussian(), data = dd,
                prior = c(set_prior("exponential(2)", class = "sigma"),
                          merge_pr0()))
  expect_identical(fixef(u), fixef(direct))
  # one empty row returns that slot, and only that slot, to flat
  bx <- tbl[tbl$class == "b" & tbl$coef == "x", ]
  bx$prior <- ""
  u1 <- update(f0, prior = bx)
  flat_x <- frm(y ~ x + z, family = gaussian(), data = dd,
                prior = set_prior("normal(0, 10)", class = "Intercept"))
  expect_identical(fixef(u1), fixef(flat_x))
  # the whole table with the sigma row edited resets every slot it lists
  t2 <- tbl
  t2$prior[t2$class == "sigma"] <- "exponential(2)"
  u2 <- update(f0, prior = t2)
  only_sig <- frm(y ~ x + z, family = gaussian(), data = dd,
                  prior = set_prior("exponential(2)", class = "sigma"))
  expect_identical(fixef(u2), fixef(only_sig))
  # a prior with no rows keeps the stored prior, as NULL does
  u0 <- update(f0, prior = tbl[0, ])
  expect_identical(fixef(u0), fixef(f0))
  expect_equal(nrow(as.data.frame(u0$prior)), 2L)
})

test_that("a named list merges by name and cannot mix with set_prior()", {
  dd <- merge_data()
  f0 <- frm(y ~ x + z, family = gaussian(), data = dd,
            prior = list(beta = prior_normal(0, 0.05)))
  u <- update(f0, prior = list(betad = prior_normal(0, 1)))
  direct <- frm(y ~ x + z, family = gaussian(), data = dd,
                prior = list(betad = prior_normal(0, 1),
                             beta = prior_normal(0, 0.05)))
  expect_identical(fixef(u), fixef(direct))
  expect_error(update(f0, prior = set_prior("exponential(2)",
                                            class = "sigma")),
               "one spelling")
})
