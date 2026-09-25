## session =: every family's value store starts again at each session
## boundary, and the subject stays the unit a random effect and the
## importance correction group on.
##
## Every value below is computed inside tryCatch(), so on a build
## without `session =` each assertion FAILS rather than the block
## stopping at its first error with the rest unrun.

ss_try <- function(expr) tryCatch(expr, error = function(e) NULL)

ss_data <- function(ns = 6L, nt = 30L, seed = 71L) {
  d1 <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                        seed = seed)
  d2 <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                        seed = seed + 1L)
  d1$session <- "a"
  d2$session <- "b"
  d <- rbind(d1, d2)
  d$idsess <- interaction(d$id, d$session, lex.order = TRUE)
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = idsess, trial = trial), d,
    pars = list(alpha = 0.4, tau = 3), seed = seed)[[1L]]$choice
  d
}

ss_f <- frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1)

test_that("two sessions fit as two single-session fits sharing parameters", {
  d <- ss_data()
  fa <- ss_try(frmtmb::frm(
    ss_f, family = bandit2arm_delta(subject = id, trial = trial,
                                    session = session), data = d))
  f1 <- frmtmb::frm(ss_f, family = bandit2arm_delta(subject = id,
                                                    trial = trial),
                    data = d[d$session == "a", ])
  f2 <- frmtmb::frm(ss_f, family = bandit2arm_delta(subject = id,
                                                    trial = trial),
                    data = d[d$session == "b", ])
  p <- if (is.null(fa)) NULL else fa$obj$env$last.par.best
  two <- if (is.null(p)) 0 else f1$obj$fn(p) + f2$obj$fn(p)
  one <- if (is.null(p)) NA_real_ else fa$obj$fn(p)
  # the objective at the SAME parameter vector: a sum of the two
  # sessions' objectives, which is what "the store restarts and the
  # parameters are shared" means, to the rounding of a reordered sum
  expect_equal(one, two, tolerance = 8 * .Machine$double.eps)
})

test_that("session = is the same likelihood as one learner per session", {
  # with no random effect, a subject's sessions share nothing but the
  # fixed parameters, so relabelling each session as its own learner
  # is the same model
  d <- ss_data()
  fa <- ss_try(frmtmb::frm(
    ss_f, family = bandit2arm_delta(subject = id, trial = trial,
                                    session = session), data = d))
  fb <- frmtmb::frm(ss_f, family = bandit2arm_delta(subject = idsess,
                                                    trial = trial), data = d)
  la <- if (is.null(fa)) NA_real_ else as.numeric(stats::logLik(fa))
  expect_equal(la, as.numeric(stats::logLik(fb)),
               tolerance = 8 * .Machine$double.eps)
})

test_that("the value store is back at its initial values at every boundary", {
  d <- ss_data()
  fa <- ss_try(frmtmb::frm(
    ss_f, family = bandit2arm_delta(subject = id, trial = trial,
                                    session = session), data = d))
  tr <- if (is.null(fa)) NULL else ss_try(frm_value_trace(fa))
  first <- if (is.null(tr)) NULL else tr[tr$trial == 1, ]
  # one first trial per subject per session, and every one of them
  # chose on the initial values, which are zero for this family
  expect_equal(NROW(first), 12L)
  expect_true(!is.null(first) && all(first$q1 == 0 & first$q2 == 0))
  expect_identical(names(tr)[1:3], c("subject", "session", "trial"))
})

test_that("the boundary is not a no-op: without session = the store carries", {
  d <- ss_data()
  fa <- ss_try(frmtmb::frm(
    ss_f, family = bandit2arm_delta(subject = id, trial = trial,
                                    session = session), data = d))
  dn <- d
  dn$trial <- dn$trial + 30L * (dn$session == "b")
  fc <- frmtmb::frm(ss_f, family = bandit2arm_delta(subject = id,
                                                    trial = trial), data = dn)
  la <- if (is.null(fa)) NA_real_ else as.numeric(stats::logLik(fa))
  lc <- as.numeric(stats::logLik(fc))
  expect_true(is.finite(la) && abs(la - lc) > 1)
})

test_that("the importance block stays the subject", {
  d <- ss_data()
  fa <- ss_try(frmtmb::frm(
    ss_f, family = bandit2arm_delta(subject = id, trial = trial,
                                    session = session), data = d))
  st <- if (is.null(fa)) NULL else
    frmtmb::single_response(fa, "a fit")$family$structure
  blk <- if (is.null(fa)) NULL else frmtmb::frame_block_of(fa$frame, "choice")
  args <- if (is.null(fa)) NULL else
    list(as.numeric(fa$frame$y$choice),
         lapply(frmtmb::eval_dpars(fa)$choice, as.numeric),
         fa$frame$aterm_values$choice, 1, blk, NULL)
  grp <- if (is.null(st)) NULL else as.numeric(do.call(st$loglik_group, args))
  row <- if (is.null(st)) NULL else as.numeric(do.call(st$loglik_row, args))
  tot <- if (is.null(st)) NA_real_ else as.numeric(do.call(st$loglik, args))
  # one value per SUBJECT, not per session, in the level order the core
  # aligns on, each the sum of that subject's rows over both sessions
  expect_equal(length(grp), nlevels(factor(d$id)))
  expect_identical(blk[["group"]], factor(d$id))
  byid <- if (is.null(row)) NA_real_ else
    unname(vapply(split(row, factor(d$id)), sum, 0))
  expect_equal(grp, byid, tolerance = 8 * .Machine$double.eps * 30)
  expect_equal(sum(grp), tot, tolerance = 8 * .Machine$double.eps * 30)
})

test_that("the stacked per-subject values keep replicate-major order", {
  d <- ss_data()
  fa <- ss_try(frmtmb::frm(
    ss_f, family = bandit2arm_delta(subject = id, trial = trial,
                                    session = session), data = d))
  st <- if (is.null(fa)) NULL else
    frmtmb::single_response(fa, "a fit")$family$structure
  blk <- if (is.null(fa)) NULL else frmtmb::frame_block_of(fa$frame, "choice")
  n <- nrow(d)
  ridx <- rep.int(seq_len(n), 3L)
  stk <- function(v) if (length(v) > 1L) v[ridx] else v
  g1 <- g3 <- NA_real_
  if (!is.null(st)) {
    y <- as.numeric(fa$frame$y$choice)
    dp <- lapply(frmtmb::eval_dpars(fa)$choice, as.numeric)
    av <- fa$frame$aterm_values$choice
    g1 <- as.numeric(st$loglik_group(y, dp, av, 1, blk, NULL))
    dp3 <- lapply(dp, stk)
    # the three replicates differ, so a mixed-up order would show
    dp3$alpha <- rep(dp$alpha, length.out = n) *
      rep(c(1, 1.5, 2), each = n)
    g3 <- as.numeric(st$loglik_group(y[ridx], dp3, lapply(av, stk), 1,
                                     blk, NULL))
  }
  ng <- nlevels(factor(d$id))
  expect_equal(length(g3), 3L * ng)
  expect_equal(g3[seq_len(ng)], g1, tolerance = 8 * .Machine$double.eps * 30)
  expect_false(isTRUE(all.equal(g3[seq_len(ng)], g3[ng + seq_len(ng)])))
})

test_that("the importance correction runs with sessions and (1 | id)", {
  skip_on_cran()
  d <- ss_data(ns = 8L, nt = 30L, seed = 72L)
  f <- frmtmb::bf(choice | reward(pay1, pay2) ~ 1 + (1 | id), tau ~ 1)
  fit <- ss_try(suppressWarnings(frmtmb::frm(
    f, family = bandit2arm_delta(subject = id, trial = trial,
                                 session = session),
    data = d, importance = 24)))
  # the correction verifies the family's per-group pieces against the
  # plain objective at its first freeze, so returning at all is that
  # check passing; one effective sample size per SUBJECT
  expect_equal(length(fit$importance$ess), 8L)
})

test_that("trial numbers need only be unique within a session", {
  d <- ss_data()
  expect_error(
    frmtmb::frm(ss_f, family = bandit2arm_delta(subject = id, trial = trial),
                data = d),
    "unique within a subject")
  dd <- d[c(1, seq_len(nrow(d))), ]
  # the message is read rather than matched by expect_error(), which in
  # edition 3 rethrows a non-matching error and would stop the block
  msg <- tryCatch({
    frmtmb::frm(ss_f, family = bandit2arm_delta(subject = id, trial = trial,
                                                session = session),
                data = dd)
    ""
  }, error = conditionMessage)
  expect_match(msg, "unique within a subject's session")
})

test_that("a missing session is refused by name", {
  d <- ss_data()
  d$session[3] <- NA
  # frm()'s own na.action drops the row before the family sees it, so
  # the refusal is reached through the generative route
  msg <- tryCatch({
    frm_task_simulate(bandit2arm_delta(subject = id, trial = trial,
                                       session = session), d,
                      pars = list(alpha = 0.4, tau = 3))
    ""
  }, error = conditionMessage)
  expect_match(msg, "session variable has 1 missing value")
})

test_that("frm_task_simulate() restarts the store at each session", {
  d <- ss_data()
  a <- ss_try(frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial, session = session), d,
    pars = list(alpha = 0.4, tau = 3), seed = 9)[[1L]]$choice)
  b <- frm_task_simulate(bandit2arm_delta(subject = idsess, trial = trial),
                         d, pars = list(alpha = 0.4, tau = 3),
                         seed = 9)[[1L]]$choice
  # the same draws, because the walk visits the same sequences in the
  # same order: sessions within a subject, subjects in level order
  expect_identical(a, b)
})

test_that("every family takes session =", {
  fams <- c("bandit2arm_delta", "bandit2arm_dual",
            "bandit4arm2_kalman_filter", "igt_orl", "igt_pvl_delta",
            "prl_fictitious", "rlddm", "ts_par7")
  has <- vapply(fams, function(f) {
    "session" %in% names(formals(getExportedValue("frmtmb.learn", f)))
  }, logical(1))
  expect_true(all(has), info = paste(fams[!has], collapse = ", "))
})

test_that("rlddm() restarts at each session too", {
  skip_if_not_installed("RWiener")
  d1 <- frm_task_design("bandit2arm", n_subject = 4L, n_trial = 30L,
                        seed = 73L)
  d2 <- frm_task_design("bandit2arm", n_subject = 4L, n_trial = 30L,
                        seed = 74L)
  d1$session <- 1L
  d2$session <- 2L
  d <- rbind(d1, d2)
  d$idsess <- interaction(d$id, d$session, lex.order = TRUE)
  d <- frm_task_simulate(
    rlddm(subject = idsess, trial = trial), d,
    pars = list(alpha = 0.4, drift = 3, bs = 1.6, ndt = 0.2, bias = 0.5),
    seed = 73)[[1L]]
  f <- frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1,
                  bs ~ 1, ndt ~ 1, bias ~ 1)
  fa <- ss_try(frmtmb::frm(f, family = rlddm(subject = id, trial = trial,
                                             session = session), data = d))
  fb <- frmtmb::frm(f, family = rlddm(subject = idsess, trial = trial),
                    data = d)
  la <- if (is.null(fa)) NA_real_ else as.numeric(stats::logLik(fa))
  expect_equal(la, as.numeric(stats::logLik(fb)), tolerance = 1e-12)
})

test_that("a session label reused in two separate runs is refused", {
  # Punch round 1, minor 1. Numbered 1..90 continuously, so the trial
  # column orders the sessions: a, b, a. Accepted silently, the two
  # runs of `a` joined into one sequence and the store carried across b.
  d <- frm_task_design("bandit2arm", n_subject = 3L, n_trial = 90L,
                       seed = 75L)
  d$session <- c("a", "b", "a")[(d$trial - 1L) %/% 30L + 1L]
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = 0.4, tau = 3), seed = 75)[[1L]]$choice
  msg <- tryCatch({
    frmtmb::frm(ss_f, family = bandit2arm_delta(subject = id, trial = trial,
                                                session = session),
                data = d)
    ""
  }, error = conditionMessage)
  expect_match(msg, "in two runs that are not adjacent")
  # the same labels in adjacent runs, a, a, b, fit
  d$session <- c("a", "a", "b")[(d$trial - 1L) %/% 30L + 1L]
  ok <- ss_try(frmtmb::frm(ss_f, family = bandit2arm_delta(
    subject = id, trial = trial, session = session), data = d))
  expect_false(is.null(ok))
})
