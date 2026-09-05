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

test_that("every deferred exclusion still has a live defect", {
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
  # "P" and "C" rows are permanent and carry none.
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
  expect_setequal(unique(ex$class), c("D", "P", "C"))

  live <- ex[ex$class == "D", , drop = FALSE]
  # if this ever reaches zero, every deferred defect is fixed and the
  # lists themselves should be gone, not merely empty
  expect_gt(nrow(live), 0)

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

test_that("se() leaves a residual sigma each package reports differently", {
  skip_unless_brms_fit()

  # DIVERGENCE in a number neither density uses. With y | se(s) and no
  # sigma = TRUE, the residual standard deviation beyond the known s is
  # zero, and both packages fit that model: the log-density tier's row
  # 14c is exact. What each REPORTS for the unused parameter differs.
  # brms declares sigma and holds it at 0, so it enters in quadrature
  # and switches off. frmtmb leaves its sigma dpar at the link-scale
  # zero, which the log link turns into 1 on the response scale.
  s <- brms_shape("r14c")
  expect_true(all(brms::posterior_epred(s$brmsfit, dpar = "sigma") == 0))
  expect_identical(unname(fixef(s$fit)$sigma[["(Intercept)"]]), 0)
  expect_exact_num(predict(s$fit, type = "response", dpar = "sigma"),
                   rep(1, nrow(s$data)),
                   label = "frmtmb reports the unused sigma as 1")
  # and the density is the same one, which is what makes this a
  # reporting difference rather than a modeling one
  expect_exact_num(brms::log_lik(s$brmsfit)[1, ], frm_row_loglik(s$fit),
                   label = "se() density agrees per row")
})

test_that("a mixture's theta is not a probability on frmtmb's scale", {
  skip_unless_brms_fit()

  # DIVERGENCE, and the same class as the zero-inflated one: a response
  # scale that is not the quantity it names. frmtmb declares theta1's
  # link as IDENTITY, so predict(type = "response", dpar = "theta1")
  # hands back the linear predictor. brms applies the softmax over the
  # component predictors, which for two components is plogis(). The
  # likelihood is unaffected, because the log-density tier proves the
  # objective is identical: the softmax IS applied inside the density.
  s <- brms_shape("r17")
  expect_identical(family(s$fit)$links$theta1$name, "identity")

  eta <- as.numeric(predict(s$fit, type = "link", dpar = "theta1"))
  expect_exact_num(predict(s$fit, type = "response", dpar = "theta1"), eta,
                   label = "frmtmb theta1 response scale is the predictor")
  be <- brms::posterior_epred(s$brmsfit, dpar = "theta1")[1, ]
  expect_exact_num(be, plogis(eta), label = "brms theta1 is softmax(eta)")
  expect_gt(max(abs(be - eta)), 0.1)

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

test_that("conditional_effects refuses a mixture with theta ~ x", {
  skip_unless_brms_fit()

  # DIVERGENCE. brms finds the covariate wherever it sits and returns a
  # panel for it; frmtmb looks only at mu1, finds y ~ 1, and stops.
  s <- brms_shape("r17")
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))
  expect_identical(names(cb), "x")
  expect_identical(nrow(cb$x), 100L)
  expect_error(conditional_effects(s$fit), "No plottable predictors")

  # named explicitly, both produce the panel and then disagree by the
  # identity-versus-softmax difference of the block above
  zb <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                   dpar = "theta1"))
  zf <- suppressWarnings(conditional_effects(s$fit, dpar = "theta1"))
  expect_identical(names(zb), "x")
  expect_identical(names(zf), "x")
  expect_gt(max(abs(zb$x$estimate__ - zf$x$estimate__)), 0.01)
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

test_that("ranef and coef key their lists differently", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")

  # DIVERGENCE, structural, and frmtmb is inconsistent with itself.
  # brms keys both lists by the GROUPING FACTOR. frmtmb keys ranef() by
  # the BLOCK and coef() by the grouping factor, so in one model
  # ranef(fit)$Subject is NULL while coef(fit)$Subject is not.
  s <- brms_shape("rC0")
  expect_identical(names(brms::ranef(s$brmsfit)), "Subject")
  expect_identical(names(coef(s$brmsfit)), "Subject")
  expect_identical(names(ranef(s$fit)), "Days | Subject")
  expect_identical(names(coef(s$fit)), "Subject")
  expect_null(ranef(s$fit)$Subject)
  expect_false(is.null(coef(s$fit)$Subject))

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
    cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))
    cf <- suppressWarnings(conditional_effects(s$fit))
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

test_that("conditional_effects plots the wrong mean when zero-inflated", {
  skip_unless_brms_fit()

  # DIVERGENCE, and the one this tier was built to catch. frmtmb's
  # DEFAULT conditional_effects() curve on a zero-inflated fit is the
  # CONDITIONAL mean exp(eta_mu), not the expected response
  # (1 - zi) * exp(eta_mu). The cause is in R/conditional-effects.R,
  # where the method = "epred", band = "wald" branch takes its point
  # estimate as lp$link$linkinv(predict(type = "link")$fit), the same
  # call that supplies the delta-method standard error. For a family
  # whose mean IS the inverse link of its mu predictor that is the
  # expected response; for one whose mean is not, it is a different
  # quantity. See dev/brms-methods-tests.md.
  for (nm in c("r16", "rC16")) {
    s <- brms_shape(nm)
    cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))$x
    cf <- suppressWarnings(conditional_effects(s$fit))$x
    nd <- data.frame(x = cf$x)
    if ("g" %in% names(s$data)) {
      nd$g <- s$data$g[1]
    }

    # brms's curve IS its posterior_epred, the expected response
    ep <- brms::posterior_epred(s$brmsfit, newdata = nd,
                                re_formula = NA)[1, ]
    expect_exact_num(cb$estimate__, ep,
                     label = paste("brms ce is epred,", nm))

    # frmtmb's is the conditional mean, exp(eta), which is brms's
    # transformed linear predictor and frmtmb's own type = "conditional"
    mu <- brms::posterior_linpred(s$brmsfit, newdata = nd,
                                  transform = TRUE, re_formula = NA)[1, ]
    expect_exact_num(cf$estimate__, mu,
                     label = paste("frmtmb ce is exp(eta),", nm))
    expect_exact_num(cf$estimate__,
                     predict(s$fit, newdata = nd, type = "conditional",
                             re.form = ~ 0),
                     label = paste("frmtmb ce is type=conditional,", nm))

    # so the default curve sits materially above the mean of Y
    expect_gt(max(cf$estimate__ / cb$estimate__ - 1), 0.15)

    # while frmtmb's OWN fitted(), predict(response) and the
    # non-default method = "predict" all give brms's answer: the
    # package disagrees with itself on one path only
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
  }
})

test_that("a mo() predictor gets a continuous grid, not its levels", {
  skip_unless_brms_fit()

  # DIVERGENCE. A monotonic effect is defined at the ordered LEVELS of
  # its variable and nowhere between them: the simplex assigns one
  # increment per step. brms plots the four levels. frmtmb builds the
  # same 100-point numeric grid it builds for any other numeric
  # predictor and evaluates the monotonic effect at 0.0303, 0.0606 and
  # so on, which the model does not define.
  #
  # The other effect in the same model, the plain numeric z, agrees
  # exactly, so this is about mo() and not about the grid machinery.
  s <- brms_shape("r2")
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))
  cf <- suppressWarnings(conditional_effects(s$fit))
  expect_setequal(names(cb), names(cf))

  expect_exact_num(cb$z$estimate__, cf$z$estimate__,
                   label = "the plain numeric effect agrees")

  expect_identical(nrow(cb$inc), 4L)
  expect_identical(cb$inc$inc, as.numeric(0:3))
  expect_identical(nrow(cf$inc), 100L)
  expect_gt(length(setdiff(cf$inc$inc, 0:3)), 90)
  # the four points both packages define agree, so the extra 96 are the
  # whole of the difference
  m <- match(cb$inc$inc, cf$inc$inc)
  expect_false(anyNA(m))
  expect_exact_num(cb$inc$estimate__, cf$inc$estimate__[m],
                   label = "mo() agrees at the levels it is defined on")
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

test_that("the conditional_effects defect is live on hurdle too", {
  # NO Stan, no brms fit. The defect is a disagreement between two of
  # frmtmb's OWN methods, so it needs neither, and it is the reason this
  # block is gated on nothing but CRAN.
  #
  # conditional_effects() takes its point estimate as the inverse link
  # of the mu predictor; a hurdle family's mean is not that, so the
  # curve is the conditional mean of the truncated count component
  # rather than the expected response. Unlike the zero-inflated case the
  # error changes SIGN along the curve, and unlike the zero-inflated
  # case method = "predict" is not a way out: hurdle_poisson has no
  # simulator. See dev/brms-methods-tests.md finding 1b.
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

  # the curve IS the conditional mean, exactly, which is what says the
  # cause is the same one
  expect_exact_num(ce$estimate__,
                   predict(fh, newdata = nd, type = "conditional"),
                   label = "hurdle ce is type=conditional")

  # and it is not the expected response, on either side of the grid
  ratio <- ce$estimate__ / epred
  expect_lt(min(ratio), 0.9)
  expect_gt(max(ratio), 3)

  # frmtmb's own fitted() and predict(response) agree with each other,
  # so the package disagrees with itself on the plotting path alone
  expect_exact_num(fitted(fh), predict(fh, type = "response"),
                   label = "hurdle fitted() is predict(response)")

  # the zero-inflated workaround is not available here
  expect_error(conditional_effects(fh, method = "predict"),
               "needs a family with a simulator")
})

test_that("the ordinal branch forwards no dots, so nothing warns", {
  # Also frmtmb alone. Finding 3 says int_conditions is ignored but at
  # least warns; that is true of the gaussian path only. The warning
  # comes from the method = "epred", band = "wald" branch forwarding its
  # dots to predict(). The ordinal and categorical branches do not
  # forward, so an unknown argument there is discarded in silence, which
  # makes finding 6b's ignored `categorical =` genuinely silent.
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

  # The gaussian path warns and names predict(). It warns ONCE PER
  # predict() call, and a wald band makes two, so the warnings are
  # captured rather than matched with expect_warning(), which would
  # consume the first and let the second escape the test.
  w <- capture_warnings(conditional_effects(fg, nosucharg = 1))
  expect_gte(length(w), 1L)
  expect_true(all(grepl("unknown arguments to predict\\(\\): nosucharg",
                        w)))

  # the ordinal path does not warn at all, for any of the three
  expect_silent(conditional_effects(fo, nosucharg = 1))
  expect_silent(conditional_effects(fo, categorical = TRUE))
  expect_silent(conditional_effects(fo,
                                    int_conditions = list(x = c(-1, 1))))

  # and categorical = TRUE really is a no-op, not merely a quiet one
  expect_identical(conditional_effects(fo, categorical = TRUE)$x,
                   conditional_effects(fo)$x)
})

test_that("conditional_effects returns different columns", {
  skip_unless_brms_fit()

  # DIVERGENCE, structural, and it is the data-level root of the
  # faceting defect in dev/brms-vignette-audit.md. brms's data frame
  # carries the held-constant covariates, cond__ and effect1__;
  # frmtmb's carries the varying predictor and the band only. Any brms
  # code that facets on cond__ or reads effect1__ has nothing to read.
  s <- brms_shape("r1")
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit))$x
  cf <- suppressWarnings(conditional_effects(s$fit))$x

  expect_identical(names(cb),
                   c("x", "y", "z", "cond__", "effect1__",
                     "estimate__", "se__", "lower__", "upper__"))
  expect_identical(names(cf),
                   c("x", "estimate__", "se__", "lower__", "upper__"))
  expect_false("cond__" %in% names(cf))
  expect_false("effect1__" %in% names(cf))
  # the held value brms records IS the mean of the other covariate, and
  # frmtmb agrees on the estimate, so it holds it there too without
  # saying so
  expect_exact_num(unique(cb$z), mean(s$data$z),
                   label = "brms holds a numeric covariate at its mean")
})

test_that("the two-way grid differs in order and in held value", {
  skip_unless_brms_fit()

  # TWO DIVERGENCES on one call.
  #
  # 1. Row order. brms varies the SECOND effect fastest; frmtmb varies
  #    the first. The two frames hold the same 300 points and no
  #    elementwise comparison of them is meaningful.
  # 2. The held values themselves. brms evaluates at mean +- sd exactly
  #    and rounds only the LABEL it puts in effect2__. frmtmb rounds the
  #    VALUE, signif(mean +- sd, 3) at R/conditional-effects.R
  #    ce_second_values(), and evaluates there, so its curve is the
  #    model at a slightly different covariate value.
  s <- brms_shape("r1")
  cb <- suppressWarnings(
    brms::conditional_effects(s$brmsfit, effects = "x:z"))[["x:z"]]
  cf <- suppressWarnings(conditional_effects(s$fit,
                                             effects = "x:z"))[["x:z"]]
  expect_identical(nrow(cb), nrow(cf))

  # order: brms repeats x while z moves, frmtmb repeats z while x moves
  expect_identical(cb$x[1], cb$x[2])
  expect_false(isTRUE(all.equal(cf$x[1], cf$x[2])))
  expect_identical(cf$z[1], cf$z[2])

  # held values: exact on one side, signif(, 3) on the other
  zb <- sort(unique(cb$z))
  zf <- sort(unique(cf$z))
  expect_exact_num(zb, sort(mean(s$data$z) +
                              c(-1, 0, 1) * sd(s$data$z)),
                   label = "brms holds the second effect at mean +- sd")
  expect_identical(zf, sort(signif(mean(s$data$z) +
                                     c(-1, 0, 1) * sd(s$data$z), 3)))
  expect_false(isTRUE(all.equal(zb, zf)))

  # and the estimates then differ, at the size of the rounding, once
  # the rows are aligned so that order is not what is being measured
  # WHEN THE ROUNDING IS FIXED, expect_false(anyNA(m)) below fails too,
  # and for the right reason: the key rounds brms's z to match frmtmb's
  # rounded one, so an unrounded cf$z stops matching. That failure reads
  # like an unrelated alignment bug, so it is named here. Four
  # expectations flip together on this fix, not three.
  kb <- paste(format(cb$x, digits = 12), format(signif(cb$z, 3),
                                                digits = 12))
  kf <- paste(format(cf$x, digits = 12), format(cf$z, digits = 12))
  m <- match(kf, kb)
  expect_false(anyNA(m))
  gap <- max(abs(cb$estimate__[m] - cf$estimate__))
  expect_gt(gap, 1e-6)
  expect_lt(gap, 1e-2)

  # brms's effect2__ is the ROUNDED LABEL of an unrounded value, which
  # is the distinction frmtmb collapses
  expect_true(is.factor(cb$effect2__))
  # brms orders the levels DESCENDING, so that a legend reads from
  # the top of the plot downward
  expect_identical(levels(cb$effect2__),
                   as.character(rev(sort(round(zb, 2)))))
})

test_that("the two-way order differs for a factor moderator too", {
  skip_unless_brms_fit()

  # The order divergence and the rounding one are independent, and a
  # FACTOR moderator separates them: there is nothing to round, both
  # packages take sort(unique(f)), and the row order still differs.
  s <- brms_shape("rfac")
  cb <- suppressWarnings(
    brms::conditional_effects(s$brmsfit, effects = "x:f"))[["x:f"]]
  cf <- suppressWarnings(conditional_effects(s$fit,
                                             effects = "x:f"))[["x:f"]]
  expect_identical(nrow(cb), nrow(cf))
  expect_setequal(as.character(cb$f), as.character(cf$f))
  expect_identical(cb$x[1], cb$x[2])
  expect_false(isTRUE(all.equal(cf$x[1], cf$x[2])))

  # aligned on (x, f) the estimates are exact, which is what says the
  # only difference here is order
  kb <- paste(format(cb$x, digits = 12), as.character(cb$f))
  kf <- paste(format(cf$x, digits = 12), as.character(cf$f))
  m <- match(kf, kb)
  expect_false(anyNA(m))
  expect_exact_num(cb$estimate__[m], cf$estimate__,
                   label = "ce x:f estimate once aligned")
})

test_that("conditional_effects(int_conditions =) is accepted and ignored", {
  skip_unless_brms_fit()

  # DIVERGENCE. int_conditions is not an argument of
  # conditional_effects.frmtmb_fit and the string does not occur
  # anywhere under R/, so it lands in ... and has no effect on the grid.
  # It is not swallowed in silence: conditional_effects() forwards its
  # dots to predict(), and predict() warns. What the warning says is
  # "ignoring unknown arguments to predict(): int_conditions", naming a
  # function the user did not call, so the message does not connect the
  # ignored argument to the plot that came back wrong.
  s <- brms_shape("r1")
  ic <- list(z = c(-1, 0, 1))
  cb <- suppressWarnings(brms::conditional_effects(
    s$brmsfit, effects = "x:z", int_conditions = ic))[["x:z"]]
  expect_identical(sort(unique(cb$z)), c(-1, 0, 1))

  expect_warning(conditional_effects(s$fit, effects = "x:z",
                                     int_conditions = ic),
                 "unknown arguments to predict\\(\\): int_conditions")

  base <- suppressWarnings(conditional_effects(s$fit,
                                               effects = "x:z"))[["x:z"]]
  cf <- suppressWarnings(conditional_effects(s$fit, effects = "x:z",
                                             int_conditions = ic))[["x:z"]]
  expect_identical(cf$z, base$z)
  expect_identical(cf$estimate__, base$estimate__)
  expect_false(any(unique(cf$z) %in% c(-1, 1)))
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
  # frmtmb DOES emit cond__ here, and only here: the column exists when
  # conditions are given and not otherwise, while brms always has it
  expect_true("cond__" %in% names(cf))
  expect_identical(as.character(unique(cf$cond__)),
                   as.character(unique(cb$cond__)))
  # but the conditioning variable itself is still not carried
  expect_false("z" %in% names(cf))
  expect_true("z" %in% names(cb))
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

test_that("an ordinal fit's conditional_effects differ three ways", {
  skip_unless_brms_fit()

  # DIVERGENCE, and it is three at once. frmtmb's DEFAULT for an
  # ordinal fit is brms's categorical = TRUE layout; frmtmb accepts and
  # ignores categorical =, so the other layout cannot be asked for; and
  # the effect is keyed "x" by frmtmb and "x:cats__" by brms, so
  # neither result can be indexed with the other's name. See
  # dev/brms-methods-tests.md.
  s <- brms_shape("r12a")
  bdef <- suppressWarnings(brms::conditional_effects(s$brmsfit))
  bcat <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                     categorical = TRUE))
  fdef <- suppressWarnings(conditional_effects(s$fit))
  fcat <- suppressWarnings(conditional_effects(s$fit, categorical = TRUE))

  expect_identical(names(bdef), "x")
  expect_identical(names(bcat), "x:cats__")
  expect_identical(names(fdef), "x")
  expect_identical(names(fcat), "x")

  # categorical = is accepted and does nothing
  expect_identical(fcat$x, fdef$x)

  # frmtmb's default IS brms's categorical layout: per-category
  # probabilities on a three-times-longer grid, with a cats__ column
  expect_identical(nrow(fdef$x), nrow(bcat[["x:cats__"]]))
  expect_true("cats__" %in% names(fdef$x))
  expect_false("cats__" %in% names(bdef$x))
  expect_setequal(as.character(fdef$x$cats__),
                  as.character(bcat[["x:cats__"]]$cats__))

  # and the probabilities agree once the rows are aligned, which is
  # what makes the difference a layout and not an arithmetic one
  kb <- paste(format(bcat[["x:cats__"]]$x, digits = 12),
              as.character(bcat[["x:cats__"]]$cats__))
  kf <- paste(format(fdef$x$x, digits = 12), as.character(fdef$x$cats__))
  m <- match(kf, kb)
  expect_false(anyNA(m))
  expect_exact_num(bcat[["x:cats__"]]$estimate__[m], fdef$x$estimate__,
                   label = "ordinal ce probabilities once aligned")

  # brms's default summary, the expected category number, has no
  # frmtmb spelling: its estimates are on the category scale
  expect_gt(min(bdef$x$estimate__), 1)
  expect_lt(max(bdef$x$estimate__), 3)
})

test_that("conditional_effects(method =) uses a different vocabulary", {
  skip_unless_brms_fit()

  # DIVERGENCE in the argument, not in the answer. Where both
  # spellings resolve the values are identical.
  s <- brms_shape("rfac")
  expect_error(conditional_effects(s$fit, effects = "x",
                                   method = "posterior_epred"),
               "should be one of")
  b <- suppressWarnings(brms::conditional_effects(
    s$brmsfit, effects = "x", method = "posterior_epred"))$x
  f <- suppressWarnings(conditional_effects(s$fit, effects = "x",
                                            method = "epred"))$x
  expect_exact_num(b$estimate__, f$estimate__,
                   label = "ce epred under both spellings")
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

test_that("re_formula: prediction agrees, conditional_effects does not", {
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

  # DIVERGENCE in conditional_effects, and only there. At the default
  # re_formula = NA the two are identical.
  cbn <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                    re_formula = NA))$Days
  cfn <- suppressWarnings(conditional_effects(s$fit,
                                              re_formula = NA))$Days
  expect_exact_num(cbn$estimate__, cfn$estimate__,
                   label = "ce re_formula = NA")

  # At re_formula = NULL brms conditions on a NEW group, which is why
  # its grouping column is NA and its curve moves from call to call.
  # frmtmb conditions on the FIRST observed level and does not say so.
  cb <- suppressWarnings(brms::conditional_effects(s$brmsfit,
                                                   re_formula = NULL))$Days
  cf <- suppressWarnings(conditional_effects(s$fit,
                                             re_formula = NULL))$Days
  expect_true("Subject" %in% names(cb))
  expect_true(all(is.na(cb$Subject)))
  expect_false("Subject" %in% names(cf))

  # frmtmb's curve IS level one's, exactly
  fe <- fixef(s$fit)$mu
  re <- ranef(s$fit)[[1]]
  lvl1 <- rownames(re)[1]
  expect_exact_num(cf$estimate__,
                   fe[["(Intercept)"]] + re[lvl1, "(Intercept)"] +
                     (fe[["Days"]] + re[lvl1, "Days"]) * cf$Days,
                   label = "frmtmb ce re_formula = NULL is level one")
  # and it is not the population curve, so the choice is visible
  expect_gt(max(abs(cf$estimate__ - cfn$estimate__)), 1)
})

test_that("an unknown argument is reported against predict()", {
  skip_unless_brms_fit()
  skip_if_not_installed("lme4")

  # frmtmb spells brms's re_formula as re.form on predict(), and brms's
  # spelling reaches ... and is dropped with a warning that names it.
  # That is the right behavior for a direct call.
  s <- brms_shape("rC0")
  expect_warning(predict(s$fit, type = "response", re_formula = NULL),
                 "ignoring unknown arguments to predict\\(\\): re_formula")

  # conditional_effects() forwards its dots to the same predict(), so an
  # argument it does not know is reported the same way, against
  # predict(), a function the user did not call. The warning is there;
  # what it does not say is that the plot ignored what was asked for.
  expect_warning(conditional_effects(s$fit, effects = "Days",
                                     int_conditions = list(Days = c(1, 2))),
                 "ignoring unknown arguments to predict\\(\\)")
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
