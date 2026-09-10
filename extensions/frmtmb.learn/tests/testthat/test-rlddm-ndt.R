## Item 1.0b of dev/extension-gaps-plan.md: the non-decision time's
## bound in rlddm().
##
## Through 0.3.0 the bound was ONE number, the fastest response in the
## whole data set, and rlddm() built the link for it itself. With 100
## learners that number is the fastest response of all of them, so a
## subject deviation on `ndt` is a deviation on a fraction of somebody
## else's floor and a learner whose true non-decision time is above the
## global minimum cannot be expressed at all. The scale tier measured
## the consequence: a maximum gradient of 6.36e+09, a Hessian that was
## not positive definite and four NaN standard errors.
##
## What changed: the bound comes from frmtmb.eam::ndt_bound() and
## ndt_group() makes it per group. A model without ndt_group() is the
## model 0.3.0 fitted, which dev/rlddm-findings.md records as identical
## in every digit across the two builds.

# 6 learners in two groups, whose true non-decision times are far
# enough apart that the SLOW group's truth sits ABOVE the fast group's
# fastest response. One bound over the whole data set cannot represent
# that, which is what makes this a before/after rather than a
# vocabulary change.
rlddm_two_group <- function(seed = 909L, ns = 6L, nt = 90L) {
  truth <- rep(c(0.18, 0.40), each = ns / 2L)
  d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                       seed = seed)
  i <- as.integer(d$id)
  s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
                         pars = list(alpha = 0.4, drift = 3, bs = 1.6,
                                     ndt = truth[i], bias = 0.5),
                         seed = seed)[[1L]]
  s$grp <- factor(ifelse(as.integer(s$id) <= ns / 2L, "fast", "slow"))
  attr(s, "ndt_truth") <- c(fast = 0.18, slow = 0.40)
  s
}

rlddm_one_row_per <- function(s) {
  s[match(levels(s$grp), as.character(s$grp)), , drop = FALSE]
}

test_that("ln_ndt_at reads a fraction only when the floor is there", {
  # The three branches, and the third is the one that matters: a
  # grouped model reached without its per-row bound would read the
  # fraction as seconds and report a converged fit several times too
  # fast. Constructed with the floor ABSENT rather than merely wrong.
  f <- frmtmb.learn:::ln_ndt_at
  expect_identical(f(list(ndt = 0.2)), 0.2)
  expect_identical(f(list(ndt = 0.5, ndt_floor = c(0.4, 0.6))),
                   c(0.2, 0.3))
  expect_error(f(list(ndt = 0.5, ndt_group = c(1, 2))),
               "did not, and reading the fraction as a time")
  # and the grouping present WITH the floor is not the refusing case
  expect_identical(f(list(ndt = 0.5, ndt_group = c(1, 2),
                          ndt_floor = c(0.4, 0.6))),
                   c(0.2, 0.3))
})

test_that("a family with no ndt_group carries one bound, in the link", {
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group()
  fit <- frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1,
                bs ~ 1, ndt ~ 1, bias = 0.5),
             family = rlddm(subject = id, trial = trial), data = s)
  fam <- frmtmb::single_response(fit)[["family"]]
  bd <- fam[["ndt_bound"]]
  expect_null(bd[["floors"]])
  expect_identical(bd[["ub"]], min(s$rt))
  expect_identical(fam[["links"]][["ndt"]][["name"]], "scaled_logit")
  # `ndt` is a TIME on the response scale, so the two agree
  p <- as.numeric(suppressWarnings(
    stats::predict(fit, dpar = "ndt", type = "response")))
  expect_equal(as.numeric(suppressWarnings(frmtmb.eam::ndt_time(fit))), p)
  expect_true(all(p < min(s$rt)))
  # no per-row bound is added to the addition-term values at all
  expect_false("ndt_floor" %in%
                 names(fit$frame$aterm_values[["rt"]]))
})

test_that("ndt_group makes the bound each group's own fastest response", {
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group()
  fit <- frm(bf(rt | dec(choice) + reward(pay1, pay2) +
                  ndt_group(grp) ~ 1, drift ~ 1, bs ~ 1, ndt ~ 0 + grp,
                bias = 0.5),
             family = rlddm(subject = id, trial = trial), data = s)
  fam <- frmtmb::single_response(fit)[["family"]]
  bd <- fam[["ndt_bound"]]
  own <- as.numeric(tapply(s$rt, s$grp, min))
  expect_equal(sort(unname(bd[["floors"]])), sort(own))
  expect_identical(unname(bd[["sizes"]]),
                   unname(as.integer(sort(table(s$grp)))))
  # frm() normalizes whatever family_finalize() returns through
  # get_link(), so the plain string arrives as the link object
  expect_identical(fam[["links"]][["ndt"]][["name"]], "logit")
  # the per-row bound rides with the addition-term values
  fl <- fit$frame$aterm_values[["rt"]][["ndt_floor"]]
  expect_length(fl, nrow(s))
  expect_equal(fl, own[as.integer(s$grp)])
  # predict() reports the FRACTION; ndt_time() reports the time
  one <- rlddm_one_row_per(s)
  fr <- as.numeric(suppressWarnings(
    stats::predict(fit, newdata = one, dpar = "ndt", type = "response")))
  tm <- as.numeric(suppressWarnings(frmtmb.eam::ndt_time(fit, newdata = one)))
  expect_true(all(fr > 0 & fr < 1))
  expect_equal(tm, fr * own)
  # and every group's fitted non-decision time is below its OWN floor
  expect_true(all(tm < own))
})

test_that("the per-group bound fits a group the global bound cannot", {
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group()
  truth <- attr(s, "ndt_truth")
  own <- as.numeric(tapply(s$rt, s$grp, min))
  one <- rlddm_one_row_per(s)
  # the design's own premise, asserted rather than assumed: the slow
  # group's truth is above the GLOBAL fastest response and below its own
  expect_gt(truth[["slow"]], min(s$rt))
  expect_true(all(truth < own))

  glob <- frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1,
                 bs ~ 1, ndt ~ 0 + grp, bias = 0.5),
              family = rlddm(subject = id, trial = trial), data = s)
  grp <- frm(bf(rt | dec(choice) + reward(pay1, pay2) +
                  ndt_group(grp) ~ 1, drift ~ 1, bs ~ 1, ndt ~ 0 + grp,
                bias = 0.5),
             family = rlddm(subject = id, trial = trial), data = s)
  hat_g <- as.numeric(suppressWarnings(
    stats::predict(glob, newdata = one, dpar = "ndt", type = "response")))
  hat_p <- as.numeric(suppressWarnings(
    frmtmb.eam::ndt_time(grp, newdata = one)))
  err_g <- abs(hat_g - truth) / truth
  err_p <- abs(hat_p - truth) / truth

  # Everything below is a comparison between two arms of the SAME run,
  # so nothing here is an absolute constant. Same data, same parameter
  # count, one bound against one per group.
  expect_identical(length(glob$opt$par), length(grp$opt$par))
  expect_gt(as.numeric(stats::logLik(grp)),
            as.numeric(stats::logLik(glob)))
  expect_lt(max(err_p), max(err_g))
  # the global arm is AT its wall: its slow group cannot go past the
  # global fastest response, and its fitted value is that number
  expect_lt(hat_g[[2L]], min(s$rt))
  expect_lt(abs(hat_g[[2L]] - min(s$rt)) / min(s$rt), max(err_p))
  # and the per-group arm goes past it, which is the whole item
  expect_gt(hat_p[[2L]], min(s$rt))
})

test_that("the grouped fit converges with usable standard errors", {
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group()
  fit <- frm(bf(rt | dec(choice) + reward(pay1, pay2) +
                  ndt_group(grp) ~ 1, drift ~ 1, bs ~ 1, ndt ~ 0 + grp,
                bias = 0.5),
             family = rlddm(subject = id, trial = trial), data = s,
             se = TRUE)
  dg <- frmtmb::diagnose(fit, quiet = TRUE)
  expect_identical(dg$convergence, 0L)
  expect_true(isTRUE(dg$pdHess))
  expect_length(dg$bad_se, 0L)
  expect_true(all(is.finite(fit$sdr$sd)))
  # the per-trial factorization is untouched by the rescaling
  tr <- frm_value_trace(fit)
  expect_equal(sum(log(tr$dens)), as.numeric(stats::logLik(fit)))
})

test_that("the per-row bound survives the stacked design", {
  # The importance correction evaluates the whole design once per draw,
  # STACKED, so `y`, every distributional parameter and every
  # addition-term value arrive nrep times as long. `ndt_floor` is one of
  # those values and is indexed by the same row numbers as the payoffs,
  # so it is stacked or not stacked with them; this constructs the case
  # rather than reasoning about it. imp_verify() compares the family's
  # per-group values against the plain objective at the first freeze, so
  # a fit returning at all is that check having passed.
  skip_on_cran()
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group(seed = 611L, ns = 6L, nt = 40L)
  fit <- suppressWarnings(frm(
    bf(rt | dec(choice) + reward(pay1, pay2) + ndt_group(id) ~ 1,
       drift ~ 1, bs ~ 1, ndt ~ 1 + (1 | id), bias = 0.5),
    family = rlddm(subject = id, trial = trial), data = s,
    importance = 16))
  expect_s3_class(fit, "frmtmb_fit")
  expect_equal(fit$importance$draws, 16)
  expect_length(fit$importance$ess, 6L)
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  own <- as.numeric(tapply(s$rt, s$id, min))
  one <- s[match(levels(s$id), as.character(s$id)), , drop = FALSE]
  expect_true(all(suppressWarnings(
    frmtmb.eam::ndt_time(fit, newdata = one)) < own))
})

test_that("the group's identity is its label, not its level index", {
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group()
  f <- bf(rt | dec(choice) + reward(pay1, pay2) + ndt_group(grp) ~ 1,
          drift ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5)
  fit <- frm(f, family = rlddm(subject = id, trial = trial), data = s)
  one <- rlddm_one_row_per(s)
  base <- as.numeric(suppressWarnings(
    frmtmb.eam::ndt_time(fit, newdata = one)))
  # droplevels() and relevel() both permute the level INDEX, which is
  # what frmtmb.eam's review found pairing every row with another
  # group's bound. The code is the LABEL's, so neither moves anything.
  slow <- droplevels(one[one$grp == "slow", , drop = FALSE])
  expect_equal(as.numeric(suppressWarnings(
    frmtmb.eam::ndt_time(fit, newdata = slow))), base[[2L]])
  rev <- one
  rev$grp <- factor(as.character(rev$grp), levels = c("slow", "fast"))
  expect_equal(as.numeric(suppressWarnings(
    frmtmb.eam::ndt_time(fit, newdata = rev))), base)
  # and a character column names the same groups, so it is the same fit
  s2 <- s
  s2$grp <- as.character(s2$grp)
  fit2 <- frm(f, family = rlddm(subject = id, trial = trial), data = s2)
  expect_equal(as.numeric(stats::logLik(fit2)),
               as.numeric(stats::logLik(fit)))
})

test_that("a bound the fit never saw is refused rather than substituted", {
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group()
  fit <- frm(bf(rt | dec(choice) + reward(pay1, pay2) +
                  ndt_group(grp) ~ 1, drift ~ 1, bs ~ 1, ndt ~ 1,
                bias = 0.5),
             family = rlddm(subject = id, trial = trial), data = s)
  one <- rlddm_one_row_per(s)
  nd <- one
  nd$grp <- factor(c("fast", "brandnew"))
  expect_error(frmtmb.eam::ndt_time(fit, newdata = nd), "was not fitted to")
  nd2 <- one
  nd2$grp <- NULL
  expect_error(frmtmb.eam::ndt_time(fit, newdata = nd2), "Supply that column")
  nd3 <- one
  nd3$grp <- factor(c("fast", NA), levels = levels(s$grp))
  expect_error(frmtmb.eam::ndt_time(fit, newdata = nd3), "every row needs one")
})

test_that("max_ndt and ndt_group together are refused", {
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group()
  expect_error(
    frm(bf(rt | dec(choice) + reward(pay1, pay2) + ndt_group(grp) ~ 1,
           drift ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
        family = rlddm(subject = id, trial = trial, max_ndt = 0.2),
        data = s, dry_run = "frame"),
    "both set the non-decision time's upper bound")
})

test_that("ndt_group on a family with no non-decision time is refused", {
  # THE CASE WHERE THE GUARDED THING IS ABSENT. Every family in this
  # package accepts every registered addition term, so a grouping
  # written on a softmax family would otherwise travel into the fit
  # unread and every row would be scored against nothing at all.
  # frmtmb.eam's frame check refuses it because no per-group table was
  # built, which is exactly the condition.
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group()
  s$pick <- s$choice + 1L
  expect_error(
    frm(bf(pick | reward(pay1, pay2) + ndt_group(grp) ~ 1, tau ~ 1),
        family = bandit2arm_delta(subject = id, trial = trial), data = s,
        dry_run = "frame"),
    "no family read it")
  expect_identical(frmtmb::frm_compat("bandit2arm_delta",
                                      "ndt_group()")$status, "refused")
  expect_identical(frmtmb::frm_compat("rlddm", "ndt_group()")$status,
                   "works")
})

test_that("a pinned ndt needs a bound, and is checked against it", {
  # WHAT THIS DOES AND DOES NOT PIN. An earlier version of this test
  # said an in-range `bf(ndt = 0.2)` was fitted at 0.0436 s through
  # 0.3.0. That was wrong and is retracted: frmtmb range-checks a
  # constant dpar at parse and TRANSFORMS it later against the settled
  # link, so 0.2 was fitted at exactly 0.2. What was wrong is the
  # check, which ran against the placeholder link and so let a constant
  # ABOVE the eventual bound through to become NaN in the objective.
  #
  # So there are two claims here and they pull opposite ways. The first
  # is a REMOVAL: `bf(ndt = )` on a bare family no longer works. The
  # second is the fix that pays for it.
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group()
  f <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift = 3,
          bs = 1.6, ndt = 0.2, bias = 0.5)
  # the removal, stated as a test so it cannot happen by accident
  expect_error(frm(f, family = rlddm(subject = id, trial = trial),
                   data = s, dry_run = "objective"),
               "bound is not set yet")
  # max_ndt settles the bound before the formula is parsed, and the
  # constant then decodes to itself where the fit reads it: out of the
  # parameter template, through the fitted family's own link
  o <- frm(f, family = rlddm(subject = id, trial = trial,
                             max_ndt = min(s$rt)),
           data = s, dry_run = "objective")
  rsp <- frmtmb::single_response(o)
  lk <- rsp[["family"]][["links"]][["ndt"]]
  hit <- Filter(function(lp) identical(lp[["dpar"]], "ndt"),
                o[["frame"]][["linpreds"]])
  eta <- as.numeric(
    o[["frame"]][["par_template"]][["betad"]][hit[[1L]][["idx"]]])
  expect_equal(lk$linkinv(eta), 0.2)

  # THE FIX. A constant above the settled bound was accepted through
  # 0.3.0, because the parse-time range check ran against the
  # placeholder `log` link where log(0.3) is finite, and reached the
  # objective as NaN. It is refused at parse now, and the message names
  # the constant and the link.
  f_bad <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift = 3,
              bs = 1.6, ndt = 2 * min(s$rt), bias = 0.5)
  expect_error(
    suppressWarnings(frm(f_bad,
                         family = rlddm(subject = id, trial = trial,
                                        max_ndt = min(s$rt)),
                         data = s, dry_run = "objective")),
    "is not in the range of")
})

test_that("a settled bound is kept when the family is finalized again", {
  # influence(), frm_simulate() and the prior-predictive path all run
  # family_finalize() over again on fewer rows. A bound re-derived there
  # would make a leave-one-out refit a refit of a DIFFERENT model.
  skip_if_not_installed("RWiener")
  s <- rlddm_two_group()
  fit <- frm(bf(rt | dec(choice) + reward(pay1, pay2) +
                  ndt_group(grp) ~ 1, drift ~ 1, bs ~ 1, ndt ~ 1,
                bias = 0.5),
             family = rlddm(subject = id, trial = trial), data = s)
  fam <- frmtmb::single_response(fit)[["family"]]
  drop <- s[-which.min(s$rt), ]
  # the dropped row IS one that SET a bound, and the bound the fewer
  # rows would give is a DIFFERENT number, so a re-derivation would be
  # visible rather than coincidentally equal
  key <- frmtmb.eam::ndt_bound_key(drop$grp)
  floors <- fam[["ndt_bound"]][["floors"]]
  expect_true(min(s$rt) %in% unname(floors))
  redrawn <- vapply(split(drop$rt, as.character(key)), min, numeric(1))
  expect_false(identical(sort(unname(redrawn)), sort(unname(floors))))
  again <- fam[["family_finalize"]](fam, drop$rt,
                                    list(ndt_group = key))
  expect_identical(again[["ndt_bound"]][["floors"]], floors)
  # and the re-finalized family adds ONE ndt_floor rather than a second
  expect_length(again[["aterm_data"]](drop$rt, list(ndt_group = key)),
                1L)
})
