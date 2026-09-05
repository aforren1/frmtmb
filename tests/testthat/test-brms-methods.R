# What the post-fit methods RETURN, against brms.
#
# The log-density tier in test-brms-likelihood.R proves that frmtmb's
# objective is the Stan program's log density at a point. It says
# nothing about fitted(), predict(), ranef(), conditional_effects() and
# the rest, and until this file nothing in the repository compared any
# of them against a brms fit: every other test that mentions brms beside
# a method is calling frmtmb's method next to a brms FORMULA object.
#
# The mechanism is a brms fit whose draws ARE frmtmb's estimates.
# Stan's Fixed_param algorithm returns its initial values as the draws,
# and the tier's translator already produces the constrained Stan-named
# parameter list that init takes, including the non-centered z blocks
# that carry frmtmb's conditional modes. Every draw is therefore
# identical and equal to frmtmb's estimate, so brms's posterior mean IS
# a point estimate and an exact comparison is legitimate: a difference
# in a returned value is a difference in the METHOD, never in the fit.
# brms_fixed_fit() in helper-brms-methods.R is that construction, and
# the first block below proves it before anything else uses it.
#
# Tolerances are relative and 1e-8. A comparison that needs anything
# looser is a finding, and a divergence is pinned by ASSERTING it, never
# by widening a tolerance or skipping a row. Every divergence this file
# asserts is written out in dev/brms-methods-tests.md with which package
# is right and what a user porting a brms script experiences.
#
# Stan compiles here, sharing the log-density tier's content-addressed
# cache, so a repository that runs both tiers compiles each program
# once. The whole file is opt-in:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")

# ---------------------------------------------------------------------
# The mechanism, proved before it is used
# ---------------------------------------------------------------------

test_that("the mechanism: brms draws ARE frmtmb's estimates", {
  skip_unless_brms_fit()

  s <- brms_shape("r1")

  # Nothing below means anything if the draws move. Fixed_param does not
  # move them, and this is the assertion that says so.
  expect_draws_degenerate(s$brmsfit)

  # The two claims the plan asks for before any comparison runs.
  pe <- brms::posterior_epred(s$brmsfit)
  expect_exact_num(pe[1, ], fitted(s$fit),
                   label = "posterior_epred vs fitted, row 1")

  ll <- brms::log_lik(s$brmsfit)
  expect_exact_num(ll[1, ], frm_row_loglik(s$fit),
                   label = "log_lik per row vs frmtmb's row density")

  # and the row sum is the log-density tier's own quantity, which ties
  # this tier to that one at a single number
  expect_lt(abs(sum(ll[1, ]) - as.numeric(logLik(s$fit))), 1e-8)
})

test_that("the mechanism holds for every shape in the matrix", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  for (nm in names(brms_methods_shapes)) {
    s <- brms_shape(nm)
    expect_draws_degenerate(s$brmsfit)
    ll <- brms::log_lik(s$brmsfit)
    expect_exact_num(ll[1, ], frm_row_loglik(s$fit),
                     label = paste("log_lik per row,", nm))
    # a row-aggregation error hides in the summed identity, so the sum
    # is checked only AFTER the per-row vector has been, and only where
    # the two quantities are the same one: see the block below
    if (!nm %in% brms_re_shapes()) {
      expect_lt(abs(sum(ll[1, ]) - as.numeric(logLik(s$fit))), 1e-6,
                label = paste("summed log_lik vs logLik,", nm))
    }
  }
})

test_that("summing log_lik is not logLik once there are random effects", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  # DIVERGENCE in the quantity, not in any arithmetic, and it is the
  # trap this file's per-row comparison exists to avoid. brms's
  # log_lik() is CONDITIONAL on the group-level values in the draw, so
  # its row sum is the conditional log-likelihood at frmtmb's own
  # conditional modes. frmtmb's logLik() is the MARGINAL likelihood,
  # the Laplace approximation with the modes integrated out. The gap is
  # the Laplace correction and it is large: 75.72 nats on the merged
  # (1 | q | g) shape, 55.19 on sleepstudy.
  #
  # The per-row comparison above is unaffected, because frm_row_loglik()
  # is conditional too. What would be wrong is porting
  # sum(log_lik(fit)) as a stand-in for logLik(fit).
  for (nm in brms_re_shapes()) {
    s <- brms_shape(nm)
    cond <- sum(brms::log_lik(s$brmsfit)[1, ])
    marg <- as.numeric(logLik(s$fit))
    # conditioning on the modes cannot lower the likelihood
    expect_gte(cond, marg - 1e-6)
    expect_exact_num(cond, sum(frm_row_loglik(s$fit)),
                     label = paste("brms log_lik sum is conditional,", nm))
  }
  # and on the two shapes where the correction is not near zero it is
  # far outside any tolerance a test could carry
  expect_gt(sum(brms::log_lik(brms_shape("r7")$brmsfit)[1, ]) -
              as.numeric(logLik(brms_shape("r7")$fit)), 50)
  expect_gt(sum(brms::log_lik(brms_shape("rC0")$brmsfit)[1, ]) -
              as.numeric(logLik(brms_shape("rC0")$fit)), 50)
})

# ---------------------------------------------------------------------
# The guard on the exclusion table
# ---------------------------------------------------------------------

test_that("no deferred exclusion is a repaired defect", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  # The agreement loops below run over lists that route the diverging
  # shapes AROUND them. That routing is the one place in this tier where
  # a fix to R/ could pass silently: repair the zero-inflated defect,
  # leave r16 out of brms_ce_shapes(), and every assertion still passes
  # while the tier has stopped covering the shape it was built for.
  #
  # So each exclusion whose class is "D" carries a probe, and this block
  # fails the moment one starts agreeing, naming the list to edit. The
  # "P" and "C" rows are permanent and carry none. Every D row has now
  # been repaired and dropped, so what this block asserts today is that
  # none has come back unprobed.
  ex <- brms_exclusions()

  # the table is the single source of truth for the lists, so a typo in
  # a key would silently exclude nothing at all
  registry <- names(brms_methods_shapes)
  for (i in seq_len(nrow(ex))) {
    shp <- sub(":.*$", "", ex$key[[i]])
    expect_true(shp %in% registry,
                label = paste("exclusion key names a registered shape:",
                              ex$key[[i]]))
  }
  expect_true(all(ex$class %in% c("D", "P", "C")))

  live <- ex[ex$class == "D", , drop = FALSE]
  # Every deferred DEFECT is repaired, so this is empty and the loop
  # below runs zero times. The lists themselves stay: their remaining
  # rows are paradigm differences and design choices, which do not
  # flip. Adding a D row back arms the probe again.
  expect_identical(nrow(live), 0L)
  expect_setequal(unique(ex$class), c("P", "C"))

  for (i in seq_len(nrow(live))) {
    agrees <- brms_exclusion_agrees(live$list[[i]], live$key[[i]])
    # a defect row in a list no probe covers would pass here unprobed,
    # so an undispatched list is a failure, not a pass
    if (is.na(agrees)) {
      fail(sprintf("finding %s: no live probe dispatches for list %s()",
                   live$finding[[i]], live$list[[i]]))
    } else if (isTRUE(agrees)) {
      fail(sprintf(paste0("finding %s looks fixed: %s now agrees with ",
                          "brms, so drop its row from brms_exclusions() ",
                          "and let %s() cover it again"),
                   live$finding[[i]], live$key[[i]], live$list[[i]]))
    } else {
      succeed()
    }
  }
})

# ---------------------------------------------------------------------
# The expectation surface: epred, linpred, and the dpars
# ---------------------------------------------------------------------

test_that("posterior_epred is frmtmb's fitted() on the response scale", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  for (nm in names(brms_methods_shapes)) {
    s <- brms_shape(nm)
    pe <- brms::posterior_epred(s$brmsfit)
    b <- if (length(dim(pe)) == 3L) pe[1, , ] else pe[1, ]
    expect_exact_num(b, as.matrix(fitted(s$fit)),
                     label = paste("posterior_epred vs fitted,", nm))
    # brms's own fitted() is the summary of the same quantity, so its
    # Estimate column has to land on the same numbers
    bf <- fitted(s$brmsfit)
    est <- if (length(dim(bf)) == 3L) {
      bf[, "Estimate", ]
    } else {
      bf[, "Estimate"]
    }
    expect_exact_num(est, as.matrix(fitted(s$fit)),
                     label = paste("brms fitted Estimate vs fitted,", nm))
  }
})

test_that("posterior_linpred is frmtmb's predict(type = 'link')", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  for (nm in brms_linpred_shapes()) {
    s <- brms_shape(nm)
    pl <- brms::posterior_linpred(s$brmsfit)
    b <- if (length(dim(pl)) == 3L) pl[1, , ] else pl[1, ]
    expect_exact_num(b, as.matrix(predict(s$fit, type = "link")),
                     label = paste("posterior_linpred vs link,", nm))
  }
})

test_that("a multi-column linear predictor is one vector to frmtmb", {
  skip_unless_brms_fit()

  # DIVERGENCE in shape, not in value. Where the mu predictor has more
  # than one column per observation, as a cs() term's thresholds and a
  # categorical family's K - 1 categories both do, brms's
  # posterior_linpred() returns draws x N x columns while frmtmb's
  # predict(type = "link") returns the N-vector of the mu predictor
  # alone. The extra columns are reachable through predict(dpar = ) for
  # the categorical case and are not reachable at all for cs().
  s <- brms_shape("r13")
  pl <- brms::posterior_linpred(s$brmsfit)
  expect_identical(dim(pl), c(10L, 300L, 2L))
  expect_length(as.numeric(predict(s$fit, type = "link")), 300L)
  # the columns ARE frmtmb's per-dpar link predictions
  for (k in seq_len(2L)) {
    dp <- paste0("mu", k + 1L)
    expect_exact_num(pl[1, , k], predict(s$fit, type = "link", dpar = dp),
                     label = paste("categorical linpred column", dp))
  }

  s2 <- brms_shape("r12e")
  expect_identical(dim(brms::posterior_linpred(s2$brmsfit)),
                   c(10L, 300L, 2L))
  expect_length(as.numeric(predict(s2$fit, type = "link")), 300L)
})

test_that("transform = TRUE is the inverse link, not always the mean", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  # brms's posterior_linpred(transform = TRUE) applies the mu link's
  # inverse and stops there. frmtmb's predict(type = "response") is the
  # MEAN. For most families those are the same number, and where they
  # are not, the difference is the family's own definition rather than
  # a disagreement: trials(n) multiplies by n, zero inflation
  # multiplies by (1 - zi), a multi-category family has no single mean.
  for (nm in brms_meanlink_shapes()) {
    s <- brms_shape(nm)
    plr <- brms::posterior_linpred(s$brmsfit, transform = TRUE)
    br <- if (length(dim(plr)) == 3L) plr[1, , ] else plr[1, ]
    expect_exact_num(br, as.matrix(predict(s$fit, type = "response")),
                     label = paste("linpred(transform) vs response,", nm))
  }

  # trials(n): brms's transformed predictor is the PROBABILITY and
  # frmtmb's response scale is the expected COUNT, so the ratio is n
  s <- brms_shape("r15")
  p <- brms::posterior_linpred(s$brmsfit, transform = TRUE)[1, ]
  expect_exact_num(p * s$data$n, predict(s$fit, type = "response"),
                   label = "binomial response scale is trials * p")
})

test_that("each dpar's epred is frmtmb's response-scale prediction", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  for (nm in names(brms_methods_shapes)) {
    s <- brms_shape(nm)
    for (dp in brms_dpars_of(s)) {
      pe <- brms::posterior_epred(s$brmsfit, dpar = dp)
      b <- if (length(dim(pe)) == 3L) pe[1, , ] else pe[1, ]
      expect_exact_num(b, as.matrix(predict(s$fit, type = "response",
                                            dpar = dp)),
                       label = paste0("epred dpar=", dp, ", ", nm))
    }
  }
})

test_that("a dpar's linpred agrees where the dpar has a predictor", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  for (nm in names(brms_methods_shapes)) {
    s <- brms_shape(nm)
    for (dp in brms_dpars_of(s)) {
      if (brms_dpar_is_scalar(s, dp)) next
      pl <- brms::posterior_linpred(s$brmsfit, dpar = dp)
      bl <- if (length(dim(pl)) == 3L) pl[1, , ] else pl[1, ]
      expect_exact_num(bl, as.matrix(predict(s$fit, type = "link",
                                             dpar = dp)),
                       label = paste0("linpred dpar=", dp, ", ", nm))
    }
  }
})

test_that("se() leaves a residual sigma both packages report as absent", {
  skip_unless_brms_fit()

  # WAS finding 15. With y | se(s) and no sigma = TRUE the residual
  # standard deviation beyond the known s is zero, and both packages
  # fit that model: the log-density tier's row 14c is exact. brms
  # declares sigma and holds it at 0. frmtmb used to report the log
  # link's inverse of its mapped-out coefficient, 1, which reads as an
  # estimate of a parameter the density does not have - and disagreed
  # with frmtmb's own sigma(), which has always answered 0.
  s <- brms_shape("r14c")
  expect_true(all(brms::posterior_epred(s$brmsfit, dpar = "sigma") == 0))
  # the coefficient is still the mapped-out link-scale zero: the repair
  # is at the reporting layer, not in the fit
  expect_identical(unname(fixef(s$fit)$sigma[["(Intercept)"]]), 0)
  expect_exact_num(predict(s$fit, type = "response", dpar = "sigma"),
                   rep(0, nrow(s$data)),
                   label = "frmtmb reports the unused sigma as 0")
  expect_identical(unname(sigma(s$fit)), 0)
  # and the density is the same one, which is what makes this a
  # reporting difference rather than a modeling one
  expect_exact_num(brms::log_lik(s$brmsfit)[1, ], frm_row_loglik(s$fit),
                   label = "se() density agrees per row")
})

test_that("a mixture's theta is the softmax on the response scale", {
  skip_unless_brms_fit()

  # WAS finding 1c. theta1's LINK is identity, because that is the
  # scale the multinomial logit inside the density works on; the
  # RESPONSE scale is now the softmax over the component predictors,
  # which for two components is plogis(). The likelihood is untouched,
  # which the per-row density comparison below asserts.
  s <- brms_shape("r17")
  expect_identical(family(s$fit)$links$theta1$name, "identity")

  eta <- as.numeric(predict(s$fit, type = "link", dpar = "theta1"))
  rv <- as.numeric(predict(s$fit, type = "response", dpar = "theta1"))
  be <- brms::posterior_epred(s$brmsfit, dpar = "theta1")[1, ]
  expect_exact_num(rv, be, label = "frmtmb theta1 IS brms's theta1")
  expect_exact_num(be, plogis(eta), label = "both are softmax(eta)")
  # it is a probability now: the old response scale reached 1.112 and
  # was outside [0, 1] on 1.25% of the rows
  expect_true(all(rv > 0 & rv < 1))
  # the link scale is still the predictor, and the two still differ,
  # which is what makes the reporting scale a choice and not a no-op
  expect_gt(max(abs(be - eta)), 0.1)
  # the density is the same one: the softmax is applied inside it, and
  # the reporting scale never reaches lpdf()
  expect_exact_num(brms::log_lik(s$brmsfit)[1, ], frm_row_loglik(s$fit),
                   label = "mixture density agrees per row")

  # and mixture_probs() is a different quantity again, the posterior
  # class responsibilities given y, so it is not the missing accessor
  mp <- as.matrix(mixture_probs(s$fit))
  expect_identical(colnames(mp), c("class1", "class2"))
  expect_gt(max(abs(mp[, 1] - be)), 0.1)

  # frmtmb refuses theta2 by name, listing what it has; brms answers
  # with the reference component's fixed zero predictor
  expect_error(predict(s$fit, type = "response", dpar = "theta2"),
               "Unknown dpar")
  expect_true(all(brms::posterior_epred(s$brmsfit, dpar = "theta2") == 0))
})

test_that("conditional_effects finds a mixture's covariate on theta", {
  skip_unless_brms_fit()

  # WAS finding 1d: frmtmb enumerated mu1's predictors, found y ~ 1 and
  # refused, naming the one linear predictor it had looked at. It now
  # falls back to every dpar of the response when the selected one has
  # nothing to plot, which is the only case that changes.
  s <- brms_shape("r17")
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))
  cf <- suppressWarnings(conditional_effects(s$fit))
  expect_identical(names(cf), names(cb))
  expect_identical(nrow(cf$x), nrow(cb$x))
  expect_exact_num(cb$x$estimate__, cf$x$estimate__,
                   label = "mixture ce default")

  # and dpar = "theta1" now agrees too, the softmax being the response
  # scale on both sides
  zb <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                   dpar = "theta1"))
  zf <- suppressWarnings(conditional_effects(s$fit, dpar = "theta1"))
  expect_identical(names(zb), "x")
  expect_identical(names(zf), "x")
  expect_exact_num(zb$x$estimate__, zf$x$estimate__,
                   label = "mixture ce dpar = theta1")
})

test_that("a dpar with no predictor puts linpred on different scales", {
  skip_unless_brms_fit()

  # DIVERGENCE, and a narrow one. brms declares an unmodeled dpar as a
  # scalar on its NATURAL scale, so posterior_linpred(dpar = ) hands
  # back that scalar: there is no linear predictor for it to be the link
  # scale of. frmtmb's predict(type = "link", dpar = ) returns the link
  # scale for every dpar alike, so it returns log(sigma).
  #
  # posterior_epred(dpar = ) agrees exactly on both kinds of dpar, so
  # only the link spelling is affected. See dev/brms-methods-tests.md.
  s <- brms_shape("rfac")
  expect_true(brms_dpar_is_scalar(s, "sigma"))

  nat <- brms::posterior_linpred(s$brmsfit, dpar = "sigma")[1, ]
  expect_exact_num(nat, predict(s$fit, type = "response", dpar = "sigma"),
                   label = "brms linpred of a scalar dpar is its value")
  expect_exact_num(log(nat), predict(s$fit, type = "link", dpar = "sigma"),
                   label = "frmtmb link of a scalar dpar is its log")
  expect_gt(max(abs(nat - as.numeric(predict(s$fit, type = "link",
                                             dpar = "sigma")))), 0.1)

  # and where the dpar IS modeled the two spellings coincide
  s1 <- brms_shape("r1")
  expect_false(brms_dpar_is_scalar(s1, "sigma"))
  expect_exact_num(brms::posterior_linpred(s1$brmsfit, dpar = "sigma")[1, ],
                   predict(s1$fit, type = "link", dpar = "sigma"),
                   label = "linpred of a modeled dpar")
})

# ---------------------------------------------------------------------
# The coefficient surface
# ---------------------------------------------------------------------

test_that("fixef point estimates agree, under brms's spelling", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  for (nm in names(brms_methods_shapes)) {
    s <- brms_shape(nm)
    bfe <- brms::fixef(s$brmsfit)
    flat <- brms_flatten_fixef(s$fit)
    common <- intersect(rownames(bfe), names(flat))
    expect_gt(length(common), 0)
    expect_exact_num(bfe[common, "Estimate"], flat[common],
                     label = paste("fixef Estimate,", nm))
    # the mechanism again, read off a summary column: identical draws
    # have no spread, so anything nonzero here means the fit moved
    expect_lt(max(abs(bfe[, "Est.Error"])), 1e-12)
  }
})

test_that("ranef and coef agree at the mapped conditional modes", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  for (nm in c("r7", "rC0", "rC16")) {
    s <- brms_shape(nm)
    bre <- brms::ranef(s$brmsfit)
    expect_gt(length(bre), 0)
    for (g in names(bre)) {
      F <- brms_ranef_block(s$fit, g)
      lv <- dimnames(bre[[g]])[[1]]
      expect_setequal(lv, rownames(F))
      # brms indexes coefficient THIRD and frmtmb second, so the
      # comparison is per coefficient rather than on a flattened array.
      # drop = FALSE throughout: a one-coefficient group collapses to a
      # vector otherwise and the level names go with it.
      for (cn in dimnames(bre[[g]])[[3]]) {
        b <- bre[[g]][lv, "Estimate", cn, drop = TRUE]
        expect_exact_num(b,
                         F[lv, brms_re_coef_to_frm(s$fit, colnames(F), cn)],
                         label = paste("ranef", nm, g, cn))
      }
    }
    # coef() is fixef + ranef in both packages, so it agrees wherever
    # ranef does; what differs is the container, asserted below
    bco <- coef(s$brmsfit)
    expect_setequal(names(bco), names(bre))
  }
})

test_that("ranef and coef key their lists the same way, as brms does", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")

  # WAS finding 11b: frmtmb keyed ranef() by the BLOCK and coef() by the
  # grouping factor, so ranef(fit)$Subject was NULL in a model where
  # coef(fit)$Subject was a data frame. Both are keyed by the factor
  # now, which is brms's and lme4's key; the block label rides along on
  # the matrix so two terms on one factor stay distinguishable.
  s <- brms_shape("rC0")
  expect_identical(names(brms::ranef(s$brmsfit)), "Subject")
  expect_identical(names(coef(s$brmsfit)), "Subject")
  expect_identical(names(ranef(s$fit)), "Subject")
  expect_identical(names(coef(s$fit)), "Subject")
  expect_false(is.null(ranef(s$fit)$Subject))
  expect_false(is.null(coef(s$fit)$Subject))
  expect_identical(attr(ranef(s$fit)$Subject, "term"), "Days | Subject")
  # VarCorr() still keys by the block, which is what tells two terms on
  # one factor apart
  expect_identical(names(VarCorr(s$fit)), "Days | Subject")

  # and brms broadcasts EVERY dpar's fixed effects over every grouping
  # factor, so its coef() carries two columns that do not vary across
  # levels; frmtmb carries the two that do
  expect_identical(dimnames(coef(s$brmsfit)$Subject)[[3]],
                   c("Intercept", "Days", "sigma_Intercept",
                     "sigma_Days"))
  expect_identical(colnames(coef(s$fit)$Subject),
                   c("(Intercept)", "Days"))
  const <- coef(s$brmsfit)$Subject[, "Estimate", "sigma_Intercept"]
  expect_lt(diff(range(const)), 1e-12)
})

test_that("coef() is the same generic for two different contracts", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")

  # DIVERGENCE. brms's coef() is group-level and nothing else: fixed
  # effects broadcast over each grouping factor's levels plus that
  # level's random effect. frmtmb's returns that when there are random
  # effects and falls back to stats::coef() when there are none
  # (R/methods-fit.R:475), so its return TYPE is a function of the
  # model rather than of the generic.
  three <- list(rfac = "numeric", r1 = "list", rC0 = "data.frame")
  for (nm in names(three)) {
    s <- brms_shape(nm)
    got <- coef(s$fit)
    if (identical(three[[nm]], "numeric")) {
      # a GLM-style fit: the mu vector, as stats::coef() gives it
      expect_type(got, "double")
      expect_null(dim(got))
      expect_exact_num(got, fixef(s$fit)$mu, label = "coef() is fixef mu")
    } else if (identical(three[[nm]], "list")) {
      # a second modeled dpar: the whole fixef() list
      expect_type(got, "list")
      expect_setequal(names(got), names(fixef(s$fit)))
    } else {
      # random effects present: brms's own quantity, per group
      expect_type(got, "list")
      expect_s3_class(got[[1]], "data.frame")
    }
  }

  # and brms answers only the third of the three: with no group-level
  # effects it refuses rather than returning an empty list
  expect_error(brms::ranef(brms_shape("rfac")$brmsfit),
               "does not contain group-level effects")
  expect_length(brms::ranef(brms_shape("rC0")$brmsfit), 1)
  expect_length(ranef(brms_shape("rfac")$fit), 0)
})

test_that("VarCorr standard deviations agree", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  for (nm in c("r7", "rC0", "rC16")) {
    s <- brms_shape(nm)
    bv <- brms::VarCorr(s$brmsfit)
    fv <- unclass(VarCorr(s$fit))
    for (g in names(bv)) {
      # rownames, not names(): a one-coefficient group drops to an
      # unnamed vector under [, "Estimate"] and takes its labels with it
      cn <- rownames(bv[[g]]$sd)
      bsd <- as.numeric(bv[[g]]$sd[, "Estimate"])
      keep <- vapply(names(fv), function(z) {
        identical(brms_block_group(z), g)
      }, logical(1))
      m <- as.matrix(Reduce(brms_blockdiag, fv[keep]))
      sd <- stats::setNames(sqrt(diag(m)), colnames(m))
      want <- vapply(cn, function(z) {
        brms_re_coef_to_frm(s$fit, colnames(m), z)
      }, character(1))
      expect_exact_num(bsd, sd[want],
                       label = paste("VarCorr sd,", nm, g))
    }
  }
})

# ---------------------------------------------------------------------
# residuals(): one type is exact, the other is brms's Monte Carlo of it
# ---------------------------------------------------------------------

test_that("brms's ordinary residual is frmtmb's response residual", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  for (nm in brms_resid_shapes()) {
    s <- brms_shape(nm)
    # brms warns that a residual is not meaningful under cens(); it is
    # right, and the two packages still form the same one, which is
    # what this block measures
    rb <- suppressWarnings(residuals(s$brmsfit, method = "posterior_epred",
                                     type = "ordinary", summary = FALSE))
    b <- if (length(dim(rb)) == 3L) rb[1, , ] else rb[1, ]
    expect_exact_num(b, as.matrix(residuals(s$fit, type = "response")),
                     label = paste("ordinary vs response residual,", nm))
  }
})

test_that("the pearson residuals divide by different quantities", {
  skip_unless_brms_fit()

  # DIVERGENCE, characterized rather than tolerated. frmtmb divides the
  # response residual by the MODEL's sigma, exactly. brms divides it by
  # the standard deviation of its posterior predictive DRAWS, which is a
  # Monte Carlo estimate of the same number, so the two agree only in
  # the limit and brms's answer moves with ndraws. brms also deprecates
  # the type. See dev/brms-methods-tests.md.
  s <- brms_shape("r1")
  ord <- residuals(s$brmsfit, method = "posterior_epred",
                   type = "ordinary", summary = FALSE)[1, ]
  rf <- as.numeric(residuals(s$fit, type = "pearson"))

  # frmtmb's denominator IS the model's sigma, to the last bit
  sig <- as.numeric(predict(s$fit, type = "response", dpar = "sigma"))
  expect_exact_num(rf, ord / sig, label = "frmtmb pearson denominator")

  # brms's is not, at any draw count a test can afford, and the gap
  # shrinks like a Monte Carlo error rather than staying put
  few <- suppressWarnings(residuals(s$brmsfit,
                                    method = "posterior_epred",
                                    type = "pearson",
                                    summary = FALSE))[1, ]
  gap_few <- max(abs(few - rf))
  expect_gt(gap_few, 1e-3)
  many <- brms_fixed_cached_n(s, 2000)
  gap_many <- max(abs(suppressWarnings(
    residuals(many, method = "posterior_epred", type = "pearson",
              summary = FALSE))[1, ] - rf))
  expect_lt(gap_many, gap_few)
})

# ---------------------------------------------------------------------
# predict(): the same word for two different quantities
# ---------------------------------------------------------------------

test_that("predict() is not the same estimand in the two packages", {
  skip_unless_brms_fit()

  # DIVERGENCE. brms's predict() summarizes posterior_predict(), so it
  # is a DRAW from the response distribution and its point column is a
  # Monte Carlo mean. frmtmb's predict(type = "response") is the
  # conditional mean itself, which is brms's fitted(). A ported script
  # that calls predict(fit) gets frmtmb's fitted() semantics and loses
  # the predictive spread. See dev/brms-methods-tests.md.
  s <- brms_shape("r1")

  # what frmtmb's predict() actually is
  expect_exact_num(predict(s$fit, type = "response"), fitted(s$fit),
                   label = "frmtmb predict(response) is fitted()")
  expect_exact_num(fitted(s$brmsfit)[, "Estimate"], fitted(s$fit),
                   label = "brms fitted() is frmtmb's predict(response)")

  # and what brms's predict() is: stochastic, so two calls differ
  set.seed(1)
  p1 <- predict(s$brmsfit)[, "Estimate"]
  p2 <- predict(s$brmsfit)[, "Estimate"]
  expect_false(isTRUE(all.equal(p1, p2)))
  # the same estimand underneath, reached only as the draws grow
  many <- brms_fixed_cached_n(s, 2000)
  gap_many <- max(abs(colMeans(brms::posterior_predict(many)) -
                        as.numeric(fitted(s$fit))))
  gap_few <- max(abs(colMeans(brms::posterior_predict(s$brmsfit)) -
                       as.numeric(fitted(s$fit))))
  expect_lt(gap_many, gap_few)

  # the containers differ too: brms returns a summary matrix, frmtmb a
  # bare vector, so nothing downstream can read a column by name
  expect_identical(colnames(predict(s$brmsfit)),
                   c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_null(dim(predict(s$fit, type = "response")))
})

# ---------------------------------------------------------------------
# conditional_effects()
# ---------------------------------------------------------------------

test_that("one-way conditional_effects agree on grid and estimate", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")

  for (nm in brms_ce_shapes()) {
    s <- brms_shape(nm)
    # the SAME call to both packages: brms_ce_args() adds
    # categorical = TRUE for a polytomous family, where brms refuses
    # its own default outright (nominal) or advises against it
    # (ordinal), and frmtmb's default is the layout it advises
    args <- brms_ce_args(s)
    cb <- suppressWarnings(suppressMessages(do.call(
      brms::conditional_effects, c(list(s$brmsfit), args))))
    cf <- suppressWarnings(suppressMessages(do.call(
      conditional_effects, c(list(s$fit), args))))
    expect_identical(names(cf), names(cb))
    for (e in names(cb)) {
      if (grepl(":", e, fixed = TRUE)) next
      expect_identical(nrow(cf[[e]]), nrow(cb[[e]]))
      expect_exact_num(cb[[e]][[e]], cf[[e]][[e]],
                       label = paste("ce grid,", nm, e))
      expect_exact_num(cb[[e]]$estimate__, cf[[e]]$estimate__,
                       label = paste("ce estimate,", nm, e))
    }
  }
})

test_that("conditional_effects draws the expected response, zero-inflated", {
  skip_unless_brms_fit()

  # WAS finding 1b, the defect this tier was built to catch. The
  # method = "epred", band = "wald" branch took its point estimate as
  # lp$link$linkinv(predict(type = "link")$fit), the inverse link of the
  # MU predictor, which is the expected response only for a family
  # whose mean is that. On a zero-inflated fit it plotted exp(eta)
  # where the mean is (1 - zi) * exp(eta): 15.2% high at the first grid
  # point and 268% at the last.
  for (nm in c("r16", "rC16")) {
    s <- brms_shape(nm)
    cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))$x
    cf <- suppressWarnings(conditional_effects(s$fit))$x
    nd <- data.frame(x = cf$x)
    if ("g" %in% names(s$data)) {
      nd$g <- s$data$g[1]
    }

    # brms's curve IS its posterior_epred, the expected response, and
    # frmtmb's is now the same number
    ep <- brms::posterior_epred(s$brmsfit, newdata = nd,
                                re_formula = NA)[1, ]
    expect_exact_num(cb$estimate__, ep,
                     label = paste("brms ce is epred,", nm))
    expect_exact_num(cf$estimate__, ep,
                     label = paste("frmtmb ce is epred,", nm))

    # and it is NOT the conditional mean of the count component, which
    # is what it used to be: the two differ by a factor of 3.7 at the
    # far end of the grid, so the fix is not a rounding
    mu <- brms::posterior_linpred(s$brmsfit, newdata = nd,
                                  transform = TRUE, re_formula = NA)[1, ]
    expect_exact_num(mu, predict(s$fit, newdata = nd,
                                 type = "conditional", re.form = ~ 0),
                     label = paste("exp(eta) is type=conditional,", nm))
    expect_gt(max(mu / cf$estimate__ - 1), 0.15)

    # the routes that always agreed still do
    expect_exact_num(predict(s$fit, newdata = nd, type = "response",
                             re.form = ~ 0), ep,
                     label = paste("predict(response) is epred,", nm))
    cfp <- suppressWarnings(conditional_effects(s$fit,
                                                method = "predict"))$x
    expect_exact_num(cfp$estimate__, ep,
                     label = paste("ce method=predict is epred,", nm))
    expect_exact_num(fitted(s$fit),
                     brms::posterior_epred(s$brmsfit)[1, ],
                     label = paste("frmtmb fitted() is epred,", nm))

    # the band is the delta method over EVERY dpar's coefficients
    # jointly, so it exists and is strictly positive on a log-scale
    # mean rather than being NA or crossing zero
    expect_true(all(is.finite(cf$se__)))
    expect_true(all(cf$lower__ > 0))
    expect_true(all(cf$lower__ < cf$estimate__ & cf$estimate__ < cf$upper__))
  }
})

test_that("a mo() predictor gets its levels, not a continuous grid", {
  skip_unless_brms_fit()

  # WAS finding 17. A monotonic effect is defined at the ordered LEVELS
  # of its variable and nowhere between them: the simplex assigns one
  # increment per step. frmtmb built the same 100-point numeric grid it
  # builds for any other numeric predictor and evaluated the monotonic
  # effect at 0.0303, 0.0606 and so on, where the model has no meaning.
  # Both packages now step by one over the observed range.
  s <- brms_shape("r2")
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))
  cf <- suppressWarnings(conditional_effects(s$fit))
  expect_setequal(names(cb), names(cf))

  expect_exact_num(cb$z$estimate__, cf$z$estimate__,
                   label = "the plain numeric effect agrees")

  expect_identical(nrow(cb$inc), 4L)
  expect_identical(nrow(cf$inc), 4L)
  expect_identical(cb$inc$inc, as.numeric(0:3))
  expect_exact_num(cf$inc$inc, as.numeric(0:3), label = "mo() grid")
  expect_exact_num(cb$inc$estimate__, cf$inc$estimate__,
                   label = "mo() estimate, elementwise")
  # the plain numeric predictor in the same model still gets the
  # 100-point grid, so the rule is about mo() and not about resolution
  expect_identical(nrow(cf$z), 100L)
})

test_that("a nonlinear predictor is refused a wald band", {
  skip_unless_brms_fit()

  # DIVERGENCE. brms plots the nonlinear model with no special
  # argument. frmtmb refuses, because its band is a delta-method
  # interval and predict() has no standard error for a nonlinear
  # predictor, and it names the three ways out. The refusal is
  # deliberate and its message is good; what it costs is that the
  # brms call does not port.
  s <- brms_shape("r5")
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))
  expect_identical(names(cb), "x")
  expect_identical(nrow(cb$x), 100L)
  expect_error(conditional_effects(s$fit),
               "cannot put a wald band on a nonlinear predictor")

  # one of the routes the message names does produce the curve, and it
  # is brms's curve
  cf <- suppressWarnings(conditional_effects(s$fit, method = "predict"))
  expect_identical(names(cf), "x")
  expect_exact_num(cb$x$estimate__, cf$x$estimate__,
                   label = "nonlinear ce under method = predict")
})

test_that("the hurdle families get the expected response too", {
  # NO Stan, no brms fit: the defect was a disagreement between two of
  # frmtmb's OWN methods, so the repair needs neither. Unlike the
  # zero-inflated case the error changed SIGN along the curve (14%
  # below the expected response at one end, 237% above at the other),
  # so no eyeball check of the plot would have read as "too high".
  # See dev/brms-methods-tests.md finding 1b.
  skip_on_cran()

  set.seed(13)
  n <- 300
  dh <- data.frame(x = rnorm(n))
  dh$y <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.3 * dh$x)), 0L,
                 1L + rpois(n, exp(0.6 + 0.4 * dh$x)))
  fh <- frm(frmtmb::bf(y ~ x, hu ~ x) + frmtmb::hurdle_poisson(),
            data = dh)

  ce <- suppressWarnings(conditional_effects(fh))$x
  nd <- data.frame(x = ce$x)
  epred <- as.numeric(predict(fh, newdata = nd, type = "response"))
  cond <- as.numeric(predict(fh, newdata = nd, type = "conditional"))

  # the curve IS the expected response, exactly
  expect_exact_num(ce$estimate__, epred, label = "hurdle ce is epred")
  # and it is not the conditional mean it used to be, on either side of
  # the grid: that ratio ran 0.862 at one end and 3.371 at the other
  ratio <- cond / epred
  expect_lt(min(ratio), 0.9)
  expect_gt(max(ratio), 3)

  expect_exact_num(fitted(fh), predict(fh, type = "response"),
                   label = "hurdle fitted() is predict(response)")

  # method = "predict" reaches the same mean up to Monte Carlo error,
  # which is the independent check on the analytic one
  set.seed(7)
  pm <- suppressWarnings(conditional_effects(fh, method = "predict"))$x
  expect_lt(max(abs(pm$estimate__ / epred - 1)), 0.15)
})

test_that("an unknown argument is named against conditional_effects()", {
  # Also frmtmb alone. The warning used to come from the
  # method = "epred", band = "wald" branch forwarding its dots to
  # predict(), so it named a function the user had not called - and the
  # ordinal and categorical branches, which do not forward, discarded
  # unknown arguments in complete silence. conditional_effects() now
  # checks its own dots, before any branch.
  skip_on_cran()

  set.seed(5)
  n <- 300
  do <- data.frame(x = rnorm(n), z = rnorm(n))
  do$y <- ordered(cut(0.9 * do$x + rlogis(n),
                      breaks = c(-Inf, -1, 0.5, Inf), labels = 1:3))
  fo <- frm(frmtmb::bf(y ~ x) + frmtmb::cumulative(), data = do)
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- 1 + 0.8 * dd$x - 0.4 * dd$z + rnorm(n)
  fg <- frm(frmtmb::bf(y ~ x + z) + gaussian(), data = dd)

  # one warning per call now, not one per internal predict() call
  w <- capture_warnings(conditional_effects(fg, nosucharg = 1))
  expect_length(w, 1L)
  expect_match(w, "conditional_effects\\(\\) is ignoring unknown")
  expect_match(w, "nosucharg")

  # the ordinal path warns for the same argument, where it used to be
  # silent for every one of these three
  wo <- capture_warnings(conditional_effects(fo, nosucharg = 1))
  expect_length(wo, 1L)
  expect_match(wo, "nosucharg")

  # and the two arguments that USED to be swallowed are arguments now,
  # so neither warns
  expect_silent(conditional_effects(fo, categorical = TRUE))
  expect_silent(conditional_effects(fo, int_conditions = list(x = c(-1, 1))))

  # categorical = is not a no-op any more: FALSE is brms's default
  # layout, one curve of expected category numbers
  expect_identical(names(conditional_effects(fo, categorical = TRUE)),
                   "x:cats__")
  expect_identical(names(conditional_effects(fo, categorical = FALSE)),
                   "x")
  expect_identical(nrow(conditional_effects(fo, categorical = FALSE)$x),
                   100L)

  # allow_new_levels is a real argument of the predict() underneath and
  # is passed through rather than reported
  expect_silent(conditional_effects(fg, allow_new_levels = TRUE))
})

test_that("conditional_effects returns brms's columns", {
  skip_unless_brms_fit()

  # WAS finding 4, and it was the data-level root of the faceting
  # defect in dev/brms-vignette-audit.md: brms's plot() facets on
  # cond__, and frmtmb's frame had no such column unless conditions =
  # was passed. The frame now carries what brms's carries, in brms's
  # order: the varied predictor, the other model variables at their
  # held values, cond__, effect1__, then the band.
  s <- brms_shape("r1")
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))$x
  cf <- suppressWarnings(conditional_effects(s$fit))$x

  expect_identical(names(cb),
                   c("x", "y", "z", "cond__", "effect1__",
                     "estimate__", "se__", "lower__", "upper__"))
  expect_identical(names(cf), names(cb))

  # the held covariate is at the same value on both sides, and it is
  # the mean, which is what the agreeing estimates already implied
  expect_exact_num(unique(cb$z), mean(s$data$z),
                   label = "brms holds a numeric covariate at its mean")
  expect_exact_num(unique(cf$z), mean(s$data$z),
                   label = "frmtmb holds it there too, and says so")
  expect_exact_num(cf$effect1__, cf$x, label = "effect1__ is the effect")
  # cond__ is a one-level factor when there is one condition set, as
  # brms's is
  expect_s3_class(cf$cond__, "factor")
  expect_identical(levels(cf$cond__), levels(cb$cond__))
})

test_that("the two-way grid agrees elementwise, order included", {
  skip_unless_brms_fit()

  # TWO findings, both repaired, and the strength of this block is that
  # it needs no alignment key at all now.
  #
  # 1. Row order (finding 5). brms varies the SECOND effect fastest;
  #    frmtmb varied the first, so no elementwise comparison of the two
  #    frames was meaningful and a script that indexed rows positionally
  #    read a different point.
  # 2. The held values (finding 2). brms evaluates at mean +- sd
  #    exactly and rounds only the LABEL in effect2__. frmtmb rounded
  #    the VALUE, signif(mean +- sd, 3), and evaluated there, so its
  #    curve was the model at a covariate value nobody chose and the
  #    error was set by the coefficient rather than by anything visible.
  s <- brms_shape("r1")
  cb <- suppressWarnings(
    brms::conditional_effects(s$brmsfit, effects = "x:z"))[["x:z"]]
  cf <- suppressWarnings(conditional_effects(s$fit,
                                             effects = "x:z"))[["x:z"]]
  expect_identical(nrow(cb), nrow(cf))

  # order: the first effect is the slowest on both sides
  expect_identical(cb$x[1], cb$x[2])
  expect_identical(cf$x[1], cf$x[2])
  expect_false(isTRUE(all.equal(cf$z[1], cf$z[2])))

  # held values: exact on both sides
  zf <- sort(unique(cf$z))
  expect_exact_num(zf, sort(mean(s$data$z) +
                              c(-1, 0, 1) * sd(s$data$z)),
                   label = "frmtmb holds the second effect at mean +- sd")
  expect_exact_num(sort(unique(cb$z)), zf,
                   label = "and at the same values brms uses")

  # so the whole grid compares elementwise, with no key
  expect_exact_num(cb$x, cf$x, label = "two-way x column")
  expect_exact_num(cb$z, cf$z, label = "two-way z column")
  expect_exact_num(cb$estimate__, cf$estimate__,
                   label = "two-way estimate, elementwise")

  # the rounding lives in the LABEL, which is the distinction frmtmb
  # used to collapse: brms orders the levels DESCENDING so a legend
  # reads from the top of the plot downward
  expect_true(is.factor(cf$effect2__))
  expect_identical(levels(cf$effect2__), levels(cb$effect2__))
  expect_identical(levels(cf$effect2__),
                   as.character(rev(sort(round(zf, 2)))))
})

test_that("the two-way factor moderator agrees elementwise too", {
  skip_unless_brms_fit()

  # A FACTOR moderator separates the two findings of the block above:
  # there is nothing to round, both packages take the observed levels,
  # and the row order alone used to put the frames 1.99 apart on the x
  # column. It agrees elementwise now.
  s <- brms_shape("rfac")
  cb <- suppressWarnings(
    brms::conditional_effects(s$brmsfit, effects = "x:f"))[["x:f"]]
  cf <- suppressWarnings(conditional_effects(s$fit,
                                             effects = "x:f"))[["x:f"]]
  expect_identical(nrow(cb), nrow(cf))
  expect_identical(as.character(cb$f), as.character(cf$f))
  expect_identical(cb$x[1], cb$x[2])
  expect_identical(cf$x[1], cf$x[2])

  expect_exact_num(cb$x, cf$x, label = "ce x:f grid, elementwise")
  expect_exact_num(cb$estimate__, cf$estimate__,
                   label = "ce x:f estimate, elementwise")
  # a factor moderator's effect2__ is the factor itself on both sides,
  # not a rounded label
  expect_identical(as.character(cf$effect2__), as.character(cb$effect2__))
})

test_that("conditional_effects(int_conditions =) conditions the effect", {
  skip_unless_brms_fit()

  # WAS finding 3. int_conditions was not an argument at all: it landed
  # in ... , reached predict() through the dots, and was reported there
  # as an unknown argument to a function the user had not called, while
  # the grid it was supposed to set came back unchanged. It is how
  # brms's own vignettes pick the levels of a moderator.
  s <- brms_shape("r1")
  ic <- list(z = c(-1, 0, 1))
  cb <- suppressWarnings(brms::conditional_effects(
    s$brmsfit, effects = "x:z", int_conditions = ic))[["x:z"]]
  cf <- suppressWarnings(conditional_effects(
    s$fit, effects = "x:z", int_conditions = ic))[["x:z"]]

  expect_identical(sort(unique(cb$z)), c(-1, 0, 1))
  expect_identical(sort(unique(cf$z)), c(-1, 0, 1))
  expect_identical(nrow(cf), nrow(cb))
  expect_exact_num(cb$estimate__, cf$estimate__,
                   label = "ce int_conditions estimate, elementwise")

  # it really moved the grid: the default holds z at mean +- sd
  base <- suppressWarnings(conditional_effects(s$fit,
                                               effects = "x:z"))[["x:z"]]
  expect_false(isTRUE(all.equal(cf$z, base$z)))
  expect_false(isTRUE(all.equal(cf$estimate__, base$estimate__)))

  # a function of the observed column works too, as in brms
  qf <- suppressWarnings(conditional_effects(
    s$fit, effects = "x:z",
    int_conditions = list(z = function(v) quantile(v, c(0.1, 0.9)))))[["x:z"]]
  expect_exact_num(sort(unique(qf$z)),
                   unname(quantile(s$data$z, c(0.1, 0.9))),
                   label = "int_conditions as a function")

  # and the VARIED variable takes one as well
  vf <- suppressWarnings(conditional_effects(
    s$fit, effects = "x", int_conditions = list(x = c(-2, 0, 2))))$x
  expect_exact_num(vf$x, c(-2, 0, 2), label = "int_conditions on effect 1")
})

test_that("conditional_effects(conditions =) agrees on values", {
  skip_unless_brms_fit()

  s <- brms_shape("r1")
  cnd <- data.frame(z = c(-1, 1))
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                   effects = "x",
                                                   conditions = cnd))$x
  cf <- suppressWarnings(conditional_effects(s$fit, effects = "x",
                                             conditions = cnd))$x
  expect_identical(nrow(cf), nrow(cb))
  expect_exact_num(cb$estimate__, cf$estimate__,
                   label = "ce conditions estimate")
  expect_true("cond__" %in% names(cf))
  expect_identical(as.character(unique(cf$cond__)),
                   as.character(unique(cb$cond__)))
  # and the conditioning variable is carried, at the value it was
  # conditioned on, so the frame says what the panel is
  expect_true("z" %in% names(cf))
  expect_true("z" %in% names(cb))
  expect_identical(cf$z, cb$z)
})

test_that("conditional_effects(dpar =) enumerates different effects", {
  skip_unless_brms_fit()

  # DIVERGENCE, in the effect LIST rather than in any value. Asked for
  # sigma, brms returns a panel for every population-level predictor in
  # the model, including the ones sigma does not depend on, where the
  # curve is flat. frmtmb returns a panel only for sigma's own
  # predictors. Values agree wherever both produce a panel.
  s <- brms_shape("r1")
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                   dpar = "sigma"))
  cf <- suppressWarnings(conditional_effects(s$fit, dpar = "sigma"))
  expect_identical(names(cb), c("x", "z"))
  expect_identical(names(cf), "x")
  expect_exact_num(cb$x$estimate__, cf$x$estimate__,
                   label = "ce dpar=sigma estimate on the shared effect")
  # brms's z panel is constant, which is what makes it uninformative
  # rather than wrong
  expect_lt(diff(range(cb$z$estimate__)), 1e-12)
  # the two record the dpar in different places
  expect_identical(attr(cb$x, "response"), "sigma")
  expect_identical(attr(cf$x, "response"), "y")
  expect_identical(attr(cf$x, "dpar"), "sigma")
})

test_that("an ordinal fit's conditional_effects reach both layouts", {
  skip_unless_brms_fit()

  # WAS finding 6b, three differences at once. frmtmb's DEFAULT is
  # brms's categorical = TRUE layout; frmtmb accepted and IGNORED
  # categorical = , so the other layout could not be asked for at all;
  # and the effect was keyed "x" where brms keys the categorical layout
  # "x:cats__", so neither result could be indexed with the other's
  # name. The argument is honored now and the keys match; what stays is
  # the DEFAULT, and it is deliberate - brms's own default warns that
  # it is treating an ordered factor as continuous and asks the user to
  # set categorical = TRUE, which is what frmtmb does without being
  # asked.
  s <- brms_shape("r12a")
  bdef <- suppressWarnings(brms::conditional_effects(s$brmsfit))
  bcat <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                     categorical = TRUE))
  fdef <- suppressWarnings(conditional_effects(s$fit))
  fcat <- suppressWarnings(conditional_effects(s$fit, categorical = TRUE))
  fnum <- suppressWarnings(conditional_effects(s$fit, categorical = FALSE))

  expect_identical(names(bdef), "x")
  expect_identical(names(bcat), "x:cats__")
  expect_identical(names(fdef), "x:cats__")
  expect_identical(names(fcat), "x:cats__")
  expect_identical(names(fnum), "x")

  # the default IS the categorical layout, and it is now brms's frame
  # elementwise: same key, same rows, same order, same probabilities
  expect_identical(fcat[["x:cats__"]], fdef[["x:cats__"]])
  expect_identical(nrow(fdef[["x:cats__"]]), nrow(bcat[["x:cats__"]]))
  expect_identical(as.character(fdef[["x:cats__"]]$cats__),
                   as.character(bcat[["x:cats__"]]$cats__))
  expect_exact_num(bcat[["x:cats__"]]$estimate__,
                   fdef[["x:cats__"]]$estimate__,
                   label = "ordinal ce probabilities, elementwise")

  # categorical = FALSE is brms's default summary, the expected
  # CATEGORY NUMBER, and it agrees with brms's default exactly
  expect_identical(nrow(fnum$x), nrow(bdef$x))
  expect_exact_num(bdef$x$estimate__, fnum$x$estimate__,
                   label = "ordinal ce expected category")
  expect_gt(min(fnum$x$estimate__), 1)
  expect_lt(max(fnum$x$estimate__), 3)
  # and it is the probability-weighted category, which is the identity
  # that says what the number means
  P <- fdef[["x:cats__"]]
  k <- as.integer(P$cats__)
  expect_exact_num(as.numeric(tapply(k * P$estimate__, P$x, sum)),
                   fnum$x$estimate__,
                   label = "expected category is sum(k * p_k)")
  # its band is a delta-method band of its own, not a copy of anything
  expect_true(all(fnum$x$se__ > 0))
  expect_true(all(fnum$x$lower__ < fnum$x$estimate__))
})

test_that("a nominal per-category display needs a bootstrap band", {
  skip_unless_brms_fit()

  # PARADIGM DIFFERENCE, and the reason r13 stays out of the one-way
  # loop. The per-category display of a NOMINAL family has no
  # thresholds, so the ordinal delta method (ord_prob_se) does not
  # apply and there is no analytic standard error to draw a wald band
  # from; frmtmb refuses by name and points at band = "boot". brms
  # summarizes posterior draws and needs no Jacobian.
  s <- brms_shape("r13")
  expect_error(conditional_effects(s$fit, categorical = TRUE),
               "no analytic standard error for the category")

  # under that band the ESTIMATE is still the fit's own, not a draw
  # mean, so it is brms's curve exactly - two refits are enough to
  # show it, because only lower__/upper__/se__ come from them
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                   categorical = TRUE))
  cf <- suppressWarnings(conditional_effects(s$fit, categorical = TRUE,
                                             band = "boot", boot = 2,
                                             seed = 1))
  expect_identical(names(cf), names(cb))
  e <- names(cb)[[1]]
  expect_identical(nrow(cf[[e]]), nrow(cb[[e]]))
  expect_identical(as.character(cf[[e]]$cats__),
                   as.character(cb[[e]]$cats__))
  expect_exact_num(cb[[e]]$estimate__, cf[[e]]$estimate__,
                   label = "nominal ce probabilities, elementwise")
})

test_that("conditional_effects(method =) takes brms's vocabulary too", {
  skip_unless_brms_fit()

  # WAS a design choice with a one-line cost: a ported call spelling
  # brms's method = "posterior_epred" failed at match.arg(), which
  # named the choices but not the rename. Both spellings resolve now,
  # and to the same numbers.
  s <- brms_shape("rfac")
  b <- suppressWarnings(brms::conditional_effects(
    s$brmsfit, effects = "x", method = "posterior_epred"))$x
  f <- suppressWarnings(conditional_effects(s$fit, effects = "x",
                                            method = "epred"))$x
  fa <- suppressWarnings(conditional_effects(s$fit, effects = "x",
                                             method = "posterior_epred"))$x
  expect_exact_num(b$estimate__, f$estimate__,
                   label = "ce epred under both spellings")
  expect_identical(fa, f)
  # posterior_predict too, and the one brms name with no frmtmb display
  # is refused by name rather than by a list of two choices
  expect_identical(
    suppressWarnings(conditional_effects(s$fit, effects = "x",
                                         method = "posterior_predict"))$x$x,
    suppressWarnings(conditional_effects(s$fit, effects = "x",
                                         method = "predict"))$x$x)
  expect_error(conditional_effects(s$fit, effects = "x",
                                   method = "posterior_linpred"),
               "no frmtmb spelling")
  expect_error(conditional_effects(s$fit, effects = "x", method = "nope"),
               "should be one of")
})

# ---------------------------------------------------------------------
# hypothesis()
# ---------------------------------------------------------------------

test_that("hypothesis point estimates agree on every expression", {
  skip_unless_brms_fit()

  s <- brms_shape("r1")
  for (h in c("x = 0", "z = 0", "Intercept = 0", "x - z = 0",
              "sigma_x = 0", "sigma_Intercept = 0", "2 * x + z = 1")) {
    hb <- brms::hypothesis(s$brmsfit, h)$hypothesis$Estimate
    hf <- as.data.frame(hypothesis(s$fit, h))$estimate
    expect_exact_num(hb, hf, label = paste("hypothesis", h))
  }
})

test_that("hypothesis reaches sd and cor by brms's names", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")

  # The natural-scale group-level quantities, both spellings: the
  # class/group pair and the fully qualified brms variable name. All
  # four agree exactly, which is the best result in this file: frmtmb
  # accepts brms's naming for a quantity it stores completely
  # differently.
  s <- brms_shape("rC0")
  for (v in c("sd_Subject__Intercept", "sd_Subject__Days",
              "cor_Subject__Intercept__Days")) {
    cls <- sub("^([a-z]+)_.*$", "\\1", v)
    grp <- sub("^[a-z]+_([^_]+)__.*$", "\\1", v)
    co <- sub("^[a-z]+_[^_]+__", "", v)
    hb <- brms::hypothesis(s$brmsfit, paste0(co, " = 0"),
                           class = cls, group = grp)$hypothesis$Estimate
    hf <- as.data.frame(hypothesis(s$fit, paste0(co, " = 0"),
                                   class = cls, group = grp))$estimate
    expect_exact_num(hb, hf, label = paste("hypothesis", cls, co))
    # and the same number through the fully qualified name
    # class = NULL, or brms prefixes the name with "b_" and cannot
    # find it
    hb2 <- brms::hypothesis(s$brmsfit, paste0(v, " = 0"),
                            class = NULL)$hypothesis$Estimate
    hf2 <- as.data.frame(hypothesis(s$fit,
                                    paste0(v, " = 0")))$estimate
    expect_exact_num(hb2, hf2, label = paste("hypothesis full name", v))
    expect_exact_num(hb2, hb, label = paste("both spellings", v))
  }
})

test_that("re_formula: prediction agrees, and so does the population curve", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")

  s <- brms_shape("rC0")

  # the prediction surface agrees on BOTH settings
  expect_exact_num(brms::posterior_epred(s$brmsfit, re_formula = NA)[1, ],
                   predict(s$fit, type = "response", re.form = ~ 0),
                   label = "epred re_formula = NA")
  expect_exact_num(brms::posterior_epred(s$brmsfit,
                                         re_formula = NULL)[1, ],
                   predict(s$fit, type = "response", re.form = NULL),
                   label = "epred re_formula = NULL")

  # and so does conditional_effects at the default re_formula = NA
  cbn <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                    re_formula = NA))$Days
  cfn <- suppressWarnings(conditional_effects(s$fit,
                                              re_formula = NA))$Days
  expect_exact_num(cbn$estimate__, cfn$estimate__,
                   label = "ce re_formula = NA")

  # WAS finding 6d. At re_formula = NULL frmtmb conditioned on the
  # FIRST OBSERVED level, Subject 308, and the frame carried no
  # grouping column to say so; the curve was 86.5 away from the
  # population one at its furthest point, and which subject it belonged
  # to depended on factor level order. It now conditions on a NEW
  # group, as brms does, and the frame says NA.
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                   re_formula = NULL))$Days
  cf <- suppressWarnings(conditional_effects(s$fit,
                                             re_formula = NULL))$Days
  expect_true("Subject" %in% names(cb))
  expect_true("Subject" %in% names(cf))
  expect_true(all(is.na(cb$Subject)))
  expect_true(all(is.na(cf$Subject)))

  # PARADIGM DIFFERENCE in what a new group's curve IS. brms draws that
  # group's random effects from the fitted covariance in every
  # posterior draw, so its curve is stochastic around the population
  # one (253.28 at Days = 0 with a slope of 7.40 on one run, against
  # the population 252.86 and 10.09). A maximum-likelihood fit has the
  # MODE, which is zero, so the curve is the population curve exactly
  # and the random-effect variance goes into the band instead.
  expect_exact_num(cf$estimate__, cfn$estimate__,
                   label = "a new group's curve is the population curve")
  expect_true(all(cf$se__ > cfn$se__))
  # brms's stochastic curve lands near it but not on it
  expect_lt(max(abs(cb$estimate__ / cf$estimate__ - 1)), 0.25)
  expect_gt(max(abs(cb$estimate__ - cf$estimate__)), 1e-6)

  # an OBSERVED group is conditions = , which says which one in the
  # frame and reproduces that level exactly
  fe <- fixef(s$fit)$mu
  re <- ranef(s$fit)[["Subject"]]
  lvl1 <- rownames(re)[1]
  cg <- suppressWarnings(conditional_effects(
    s$fit, re_formula = NULL,
    conditions = list(Subject = lvl1)))$Days
  expect_identical(unique(as.character(cg$Subject)), lvl1)
  expect_exact_num(cg$estimate__,
                   fe[["(Intercept)"]] + re[lvl1, "(Intercept)"] +
                     (fe[["Days"]] + re[lvl1, "Days"]) * cg$Days,
                   label = "ce conditions = names the group exactly")
  expect_gt(max(abs(cg$estimate__ - cfn$estimate__)), 1)
})

test_that("an unknown argument is reported against the function called", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")

  # frmtmb spells brms's re_formula as re.form on predict(), and brms's
  # spelling reaches ... and is dropped with a warning that names it.
  # That is the right behavior for a DIRECT call.
  s <- brms_shape("rC0")
  expect_warning(predict(s$fit, type = "response", re_formula = NULL),
                 "ignoring unknown arguments to predict\\(\\): re_formula")

  # conditional_effects() used to forward its dots to that same
  # predict(), so an argument IT did not know was reported against a
  # function the user had not called - and only on the branches that
  # forward. It checks its own dots now, and int_conditions is one of
  # its arguments rather than an unknown one.
  expect_silent(conditional_effects(s$fit, effects = "Days",
                                    int_conditions = list(Days = c(1, 2))))
  expect_warning(conditional_effects(s$fit, effects = "Days",
                                     nosucharg = 1),
                 "conditional_effects\\(\\) is ignoring unknown")
})

test_that("hypothesis returns a different object in each package", {
  skip_unless_brms_fit()

  # DIVERGENCE, structural. brms returns a brmshypothesis LIST whose
  # $hypothesis is the table; frmtmb returns the table itself, with
  # lower-case columns. The documented brms idiom,
  # hypothesis(fit)$hypothesis$Estimate, therefore reaches a character
  # vector on a frmtmb fit and errors on the second $.
  s <- brms_shape("r1")
  hb <- brms::hypothesis(s$brmsfit, "x = 0")
  hf <- hypothesis(s$fit, "x = 0")

  expect_s3_class(hb, "brmshypothesis")
  expect_s3_class(hf, "frmtmb_hypothesis")
  expect_s3_class(hf, "data.frame")
  expect_identical(names(hb$hypothesis),
                   c("Hypothesis", "Estimate", "Est.Error", "CI.Lower",
                     "CI.Upper", "Evid.Ratio", "Post.Prob", "Star"))
  expect_identical(names(as.data.frame(hf)),
                   c("hypothesis", "estimate", "se", "lwr", "upr",
                     "z", "p"))
  expect_type(hf$hypothesis, "character")
  expect_error(hf$hypothesis$Estimate, "atomic")
  # brms rewrites the expression, frmtmb keeps it verbatim
  expect_identical(hb$hypothesis$Hypothesis, "(x) = 0")
  expect_identical(as.data.frame(hf)$hypothesis, "x = 0")
})

# ---------------------------------------------------------------------
# Structure only: the columns that depend on posterior spread
# ---------------------------------------------------------------------

test_that("the interval columns are a documented difference in kind", {
  skip_unless_brms_fit()

  # NOT a divergence to fix. brms's Est.Error and quantile columns
  # summarize DRAWS, and under Fixed_param there is nothing to
  # summarize, so they are exactly zero and the band collapses onto the
  # estimate. frmtmb's se__ is a Wald standard error from the observed
  # information, which is a frequentist quantity with no draws behind
  # it and is nonzero. The shapes agree; the meanings do not, and this
  # records that rather than forcing agreement.
  s <- brms_shape("r1")
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))$x
  cf <- suppressWarnings(conditional_effects(s$fit))$x

  expect_identical(cb$se__, rep(0, nrow(cb)))
  expect_identical(cb$lower__, cb$estimate__)
  expect_identical(cb$upper__, cb$estimate__)
  expect_true(all(cf$se__ > 0))
  expect_true(all(cf$lower__ < cf$estimate__))
  expect_true(all(cf$upper__ > cf$estimate__))

  # the same in the coefficient table
  expect_lt(max(abs(brms::fixef(s$brmsfit)[, "Est.Error"])), 1e-12)
  expect_true(all(as.data.frame(hypothesis(s$fit, "x = 0"))$se > 0))
})

test_that("loo and bayes_R2 refuse on a maximum-likelihood fit", {
  skip_unless_brms_fit()

  # NOT a divergence to fix either: frmtmb states the reason in R/loo.R
  # and names the route. What this pins is that the refusal is a refusal
  # and not a wrong number, and that brms answers where frmtmb declines.
  s <- brms_shape("r1")
  expect_error(loo(s$fit), "posterior quantity")
  expect_error(bayes_R2(s$fit), "per posterior draw")

  lb <- suppressWarnings(brms::loo(s$brmsfit))
  expect_true("elpd_loo" %in% rownames(lb$estimates))
  expect_true("SE" %in% colnames(lb$estimates))
  br2 <- suppressWarnings(brms::bayes_R2(s$brmsfit))
  expect_identical(colnames(br2),
                   c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
})

test_that("pp_check returns a plot on both sides", {
  skip_unless_brms_fit()
  skip_if_not_installed("ggplot2")

  s <- brms_shape("r1")
  pb <- suppressWarnings(brms::pp_check(s$brmsfit, ndraws = 5))
  pf <- suppressWarnings(pp_check(s$fit, ndraws = 5))
  expect_s3_class(pb, "ggplot")
  expect_s3_class(pf, "ggplot")
})

# ---------------------------------------------------------------------
# newdata
# ---------------------------------------------------------------------

test_that("newdata: both agree on the values and on what is refused", {
  skip_unless_brms_fit()

  s <- brms_shape("r1")
  dd <- s$data
  cases <- list(
    plain = dd[1:5, ],
    permuted = dd[1:5, c("z", "y", "x")],
    extra_column = cbind(dd[1:5, ], junk = 1:5)
  )
  for (cs in names(cases)) {
    pe <- brms::posterior_epred(s$brmsfit, newdata = cases[[cs]])
    fp <- predict(s$fit, newdata = cases[[cs]], type = "response")
    expect_exact_num(pe[1, ], fp, label = paste("newdata", cs))
  }
  # a column the model needs is refused by both, in different words
  drop <- dd[1:5, setdiff(names(dd), "z"), drop = FALSE]
  expect_error(brms::posterior_epred(s$brmsfit, newdata = drop))
  expect_error(predict(s$fit, newdata = drop, type = "response"))
})

test_that("newdata: dropping a factor level, and adding one", {
  skip_unless_brms_fit()

  s <- brms_shape("rfac")
  # a newdata that uses two of the three levels, with the third dropped
  # from the factor entirely, is accepted by both and gives the same
  # numbers: the contrast coding comes from the FIT, not from newdata
  nd <- droplevels(subset(s$data, f != "c"))[1:6, ]
  expect_identical(levels(nd$f), c("a", "b"))
  expect_exact_num(brms::posterior_epred(s$brmsfit, newdata = nd)[1, ],
                   predict(s$fit, newdata = nd, type = "response"),
                   label = "newdata with a level dropped")

  # a character column where the fit saw a factor is accepted by both
  ndc <- nd
  ndc$f <- as.character(ndc$f)
  expect_exact_num(brms::posterior_epred(s$brmsfit, newdata = ndc)[1, ],
                   predict(s$fit, newdata = ndc, type = "response"),
                   label = "newdata with a character column")

  # a level the fit never saw is refused by both, in different words
  ndn <- nd
  levels(ndn$f) <- c("a", "zz")
  expect_error(brms::posterior_epred(s$brmsfit, newdata = ndn),
               "New factor levels are not allowed")
  expect_error(predict(s$fit, newdata = ndn, type = "response"),
               "new levels")
})
