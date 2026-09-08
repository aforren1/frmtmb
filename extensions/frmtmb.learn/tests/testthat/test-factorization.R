## The two factorization slots, and what admitting them turned on.
##
## The families declare `loglik_row` and `loglik_group` beside `loglik`,
## all three off one call to the recursion. What has to be true for that
## to be honest rather than convenient is checked here: the three agree,
## they agree on the STACKED design the importance correction calls them
## with, and the block's `group` is the partition they claim.

ln_small <- function(ns = 6L, nt = 30L, seed = 61L, sd_u = 0) {
  d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                       seed = seed)
  set.seed(seed)
  al <- stats::plogis(stats::qlogis(0.4) +
                        stats::rnorm(ns, 0, sd_u)[as.integer(d$id)])
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = al, tau = 3), seed = seed)[[1L]]$choice
  d
}

# the two parameters on their natural scales, for the longhand
# reference below
ln_nat_alpha_tau <- function(fit) {
  b <- unlist(frmtmb::fixef(fit))
  fam <- frmtmb::single_response(fit, "a fit")$family
  c(alpha = fam$links$alpha$linkinv(b[["alpha.(Intercept)"]]),
    tau = fam$links$tau$linkinv(b[["tau.(Intercept)"]]))
}

# the three slots at the estimates, off the tape, the way the core
# calls them
ln_slots <- function(fit, resp = "choice") {
  st <- frmtmb::single_response(fit, "a fit")$family$structure
  blk <- frmtmb::frame_block_of(fit$frame, resp)
  list(st = st, blk = blk,
       dp = lapply(frmtmb::eval_dpars(fit)[[resp]], as.numeric),
       av = fit$frame$aterm_values[[resp]],
       y = as.numeric(fit$frame$y[[resp]]))
}

test_that("the total, the per-subject values and the per-row values agree", {
  d <- ln_small()
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  s <- ln_slots(fit)
  tot <- s$st$loglik(s$y, s$dp, s$av, 1, s$blk, NULL)
  grp <- s$st$loglik_group(s$y, s$dp, s$av, 1, s$blk, NULL)
  row <- s$st$loglik_row(s$y, s$dp, s$av, 1, s$blk, NULL)
  expect_equal(length(grp), nlevels(factor(d$id)))
  expect_equal(length(row), nrow(d))
  expect_equal(sum(grp), tot, tolerance = 1e-12)
  expect_equal(sum(row), tot, tolerance = 1e-12)
  # with no random effects the conditional total IS the marginal one
  expect_equal(tot, as.numeric(stats::logLik(fit)), tolerance = 1e-9)
  # and the per-row values are the trace's own per-trial factors
  expect_equal(as.numeric(row), log(frm_value_trace(fit)$p),
               tolerance = 1e-12)
})

test_that("a subject's value is the sum of its own rows and no others", {
  d <- ln_small()
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  s <- ln_slots(fit)
  grp <- as.numeric(s$st$loglik_group(s$y, s$dp, s$av, 1, s$blk, NULL))
  row <- as.numeric(s$st$loglik_row(s$y, s$dp, s$av, 1, s$blk, NULL))
  g <- s$blk[["group"]]
  expect_true(is.factor(g))
  expect_equal(length(g), nrow(d))
  # the order the core fixes: the factor's own levels
  expect_equal(unname(vapply(split(row, g), sum, 0)), grp,
               tolerance = 1e-12)
})

test_that("the slots honor the stacking the correction calls them with", {
  d <- ln_small()
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  s <- ln_slots(fit)
  grp <- as.numeric(s$st$loglik_group(s$y, s$dp, s$av, 1, s$blk, NULL))
  row <- as.numeric(s$st$loglik_row(s$y, s$dp, s$av, 1, s$blk, NULL))
  n <- s$blk[["n"]]
  nrep <- 3L
  ridx <- rep.int(seq_len(n), nrep)
  st <- function(v) if (length(v) > 1L) v[ridx] else v
  ys <- s$y[ridx]
  dps <- lapply(s$dp, st)
  avs <- lapply(s$av, st)
  g2 <- as.numeric(s$st$loglik_group(ys, dps, avs, 1, s$blk, NULL))
  r2 <- as.numeric(s$st$loglik_row(ys, dps, avs, 1, s$blk, NULL))
  # the same design repeated gives the same values repeated, in the
  # replicate-major order the protocol fixes
  expect_equal(g2, rep(grp, nrep), tolerance = 1e-12)
  expect_equal(r2, rep(row, nrep), tolerance = 1e-12)
  # and a replicate whose parameters differ must differ, or the check
  # above would pass on a family that ignored the stacking entirely
  dps2 <- dps
  dps2$alpha <- dps$alpha * rep(c(1, 1.5, 2), each = n)
  g3 <- as.numeric(s$st$loglik_group(ys, dps2, avs, 1, s$blk, NULL))
  ng <- length(grp)
  expect_equal(g3[seq_len(ng)], grp, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(g3[seq_len(ng)], g3[ng + seq_len(ng)])))
})

test_that("a nominal choice has a saturated value and a density has none", {
  d <- ln_small()
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  s <- ln_slots(fit)
  row <- s$st$loglik_row(s$y, s$dp, s$av, 1, s$blk, NULL)
  # a saturated fit puts probability one on the option that was taken
  expect_equal(attr(row, "saturated"), rep(0, nrow(d)))

  skip_if_not_installed("RWiener")
  dr <- frm_task_design("bandit2arm", n_subject = 5L, n_trial = 30L,
                        seed = 62L)
  dr <- frm_task_simulate(
    rlddm(subject = id, trial = trial), dr,
    pars = list(alpha = 0.4, drift = 3, bs = 1.6, ndt = 0.2, bias = 0.5),
    seed = 62)[[1L]]
  fr <- frmtmb::frm(
    frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1,
               bs ~ 1, ndt ~ 1, bias ~ 1),
    family = rlddm(subject = id, trial = trial), data = dr)
  sr <- ln_slots(fr, "rt")
  rr <- sr$st$loglik_row(sr$y, sr$dp, sr$av, 1, sr$blk, NULL)
  # a density has no saturated value: its supremum over the parameters
  # at a fixed response time is unbounded
  expect_null(attr(rr, "saturated"))
  expect_equal(sum(rr), as.numeric(stats::logLik(fr)), tolerance = 1e-9)
})

test_that("the importance correction runs on a learn fit", {
  skip_on_cran()
  d <- ln_small(ns = 8L, nt = 40L, seed = 63L, sd_u = 0.6)
  f <- frmtmb::bf(choice | reward(pay1, pay2) ~ 1 + (1 | id), tau ~ 1)
  fam <- bandit2arm_delta(subject = id, trial = trial)
  # WHAT THIS ASSERTS AND WHAT IT DOES NOT. The seam: the family's
  # per-group values reach the correction, and `imp_verify()` inside it
  # compares them against the plain objective per group and in total at
  # the first freeze, so a fit returning at all is that check having
  # passed. NOT the answer: eight subjects is far too few for the
  # correction to converge, and it says so rather than being quiet about
  # it. The iteration takes the SAME step for all five rounds at 24, 50,
  # 100 and 200 draws alike, which is measured in
  # dev/learn2-findings.md; the warning is asserted here so a future
  # change that made it silent would fail rather than pass.
  #
  # frmtmb >= 0.55.0 says which kind of capped this is. The step here is
  # a property of the draws: at 50 draws these eight subjects walk at
  # 0.94961, the same number a core Bernoulli fit of eight groups walks
  # at with the same seed and draw count.
  expect_warning(
    fit <- frmtmb::frm(f, family = fam, data = d, importance = 24),
    "moved by the same amount")
  expect_s3_class(fit, "frmtmb_fit")
  expect_equal(fit$importance$draws, 24)
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  expect_equal(length(fit$importance$ess), 8L)
  expect_true(fit$importance$capped)
  # the moves are equal to a part in ten thousand of themselves, which
  # is what "the same step" means and what the warning keys on
  mv <- fit$importance$moves
  expect_lt((max(mv) - min(mv)) / mean(mv), 1e-3)
})

test_that("importance refuses a grouping that is not the family's unit", {
  skip_on_cran()
  d <- ln_small(ns = 6L, nt = 20L, seed = 64L)
  # a second factor that cuts across subjects: the family's units and
  # the correction's grouping levels are then different partitions, and
  # summing one into the other means nothing
  d$block <- factor(rep_len(c("a", "b"), nrow(d)))
  expect_error(
    frmtmb::frm(
      frmtmb::bf(choice | reward(pay1, pay2) ~ 1 + (1 | block), tau ~ 1),
      family = bandit2arm_delta(subject = id, trial = trial), data = d,
      importance = 8),
    "same partition of the rows")
})

test_that("the compatibility table says importance works, with the numbers", {
  x <- frmtmb::frm_compat("bandit2arm_delta", "importance")
  expect_identical(x$status, "works")
  expect_match(x$note, "loglik_group")
  expect_match(x$note, "learn2-findings")
  # and the deviance refusal now names the sign rather than the slot
  d <- ln_small(ns = 4L, nt = 12L, seed = 65L)
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  expect_error(stats::residuals(fit, type = "deviance"),
               "SIGN rather than the magnitude")
})

# An independent longhand reference for the per-trial factors: one
# subject and one trial at a time, no mask, no vectorization. The two
# tests below compare rows against it in the two layouts the fixtures
# above do not have.
ln_ref_rows <- function(d, alpha, tau) {
  out <- rep(NA_real_, nrow(d))
  for (s in unique(d$id)) {
    r <- which(d$id == s)
    r <- r[order(d$trial[r])]
    q <- c(0, 0)
    for (i in r) {
      k <- d$choice[i]
      out[i] <- tau * q[k] - log(sum(exp(tau * q)))
      pay <- if (k == 1) d$pay1[i] else d$pay2[i]
      q[k] <- q[k] + alpha * (pay - q[k])
    }
  }
  out
}

test_that("the per-row values are right on a PADDED block", {
  # unequal trial counts, so the block really is padded and the `keep`
  # mask in ln_loglik_row() is load-bearing. Every fixture above has
  # equal counts, so none of them reaches it.
  len <- c(11L, 20L, 24L, 20L, 18L, 30L)
  d <- frm_task_design("bandit2arm", n_subject = 6L, n_trial = 30L,
                       seed = 66L)
  keep <- unlist(lapply(seq_along(len), function(s) {
    which(as.integer(d$id) == s)[seq_len(len[[s]])]
  }))
  d <- d[sort(keep), ]
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = 0.4, tau = 3), seed = 66)[[1L]]$choice
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  s <- ln_slots(fit)
  row <- as.numeric(s$st$loglik_row(s$y, s$dp, s$av, 1, s$blk, NULL))
  p <- ln_nat_alpha_tau(fit)
  expect_equal(row, ln_ref_rows(d, p[["alpha"]], p[["tau"]]),
               tolerance = 1e-12)
  # a padded cell repeats its subject's FIRST row number, so a scatter
  # that wrote the pads would leave those rows at a masked zero
  expect_false(any(row == 0))
  expect_equal(sum(row), as.numeric(stats::logLik(fit)), tolerance = 1e-9)
  expect_equal(unname(table(s$blk[["group"]])), unname(table(d$id)))
})

test_that("the per-subject values follow levels(group), not order of appearance", {
  d <- frm_task_design("bandit2arm", n_subject = 6L, n_trial = 25L,
                       seed = 67L)
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = 0.4, tau = 3), seed = 67)[[1L]]$choice
  # levels reversed AND the rows shuffled, so order of appearance and
  # level order disagree. The core fixes the order of what
  # loglik_group returns as the factor's levels; a family that returned
  # them in order of appearance would be aligned to the wrong groups
  # and nothing about the total would show it.
  d$id <- factor(d$id, levels = rev(levels(factor(d$id))))
  set.seed(67)
  d <- d[sample(nrow(d)), ]
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  s <- ln_slots(fit)
  row <- as.numeric(s$st$loglik_row(s$y, s$dp, s$av, 1, s$blk, NULL))
  grp <- as.numeric(s$st$loglik_group(s$y, s$dp, s$av, 1, s$blk, NULL))
  p <- ln_nat_alpha_tau(fit)
  expect_equal(row, ln_ref_rows(d, p[["alpha"]], p[["tau"]]),
               tolerance = 1e-12)
  by_level <- unname(vapply(split(row, s$blk[["group"]]), sum, 0))
  expect_equal(grp, by_level, tolerance = 1e-12)
  # and the check is not vacuous: the two orderings really differ
  by_appearance <- unname(vapply(
    split(row, factor(as.character(d$id), levels = unique(as.character(d$id)))),
    sum, 0))
  expect_false(isTRUE(all.equal(by_level, by_appearance)))
})
