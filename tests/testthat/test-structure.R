# The structured-family protocol's constructor and the core accessors
# that read it. The families that RIDE the protocol are covered by the
# mixture tests here and, for the two that now live outside this
# package, by the extension's own test-structure-latent.R; what is
# asserted here is the contract a family written outside this package
# must meet.

ll_ok <- function(y, dpars, aterms, weights, block, extra) sum(y)

test_that("every slot defaults to the rowwise behavior", {
  expect_error(frmtmb_structure(loglik = "not a function"), "loglik")
  st <- frmtmb_structure(loglik = ll_ok)
  expect_s3_class(st, "frmtmb_structure")
  expect_null(st[["frame_vars"]])
  expect_null(st[["fitted_mean"]])
  expect_false(st[["keep_na"]])
})

test_that("a structure with no loglik declares capabilities only", {
  # the lca() shape: a likelihood that DOES factorize per row, carried
  # by the family's own lpdf, plus refusals that belong to the family
  st <- frmtmb_structure(supports = list(osa = FALSE))
  expect_s3_class(st, "frmtmb_structure")
  expect_null(st[["loglik"]])
  expect_false(structure_allows(st, "osa"))
  # with no loglik the defaults OPEN: everything already works for a
  # rowwise family, so it names only what it refuses
  expect_true(structure_allows(st, "reml"))
  expect_true(structure_allows(st, "conditional_effects"))
  expect_identical(sum(!st[["supports"]]), 1L)
})

test_that("every capability starts refused", {
  st <- frmtmb_structure(loglik = ll_ok)
  expect_false(any(st[["supports"]]))
  expect_length(st[["supports"]], length(frmtmb_structure_flags))
  # a family with no structure at all is rowwise and refuses nothing
  expect_true(structure_allows(NULL, "osa"))
  expect_false(structure_allows(st, "osa"))
})

test_that("unit is one noun phrase, or the generic one", {
  expect_error(frmtmb_structure(loglik = ll_ok, unit = c("a", "b")),
               "noun phrase")
  expect_error(frmtmb_structure(loglik = ll_ok, unit = ""), "noun phrase")
  expect_identical(
    structure_unit(frmtmb_structure(loglik = ll_ok, unit = "a sequence")),
    "a sequence")
  expect_true(nzchar(structure_unit(frmtmb_structure(loglik = ll_ok))))
})

test_that("the constructor validates slot types", {
  for (arg in c("frame_vars", "check_spec", "frame_block", "check_frame",
                "check_fit", "fitted_mean", "fitted_var", "latent_probs",
                "sim_ctx")) {
    a <- list(loglik = ll_ok)
    a[[arg]] <- 1
    expect_error(do.call(frmtmb_structure, a), "must be a function")
  }
  expect_error(frmtmb_structure(loglik = ll_ok, keep_na = NA), "TRUE or FALSE")
  expect_error(frmtmb_structure(loglik = ll_ok, keep_na = "yes"),
               "TRUE or FALSE")
})

test_that("supports names come from a closed vocabulary", {
  expect_error(frmtmb_structure(loglik = ll_ok, supports = list(nope = TRUE)),
               "unknown capability flag")
  expect_error(frmtmb_structure(loglik = ll_ok, supports = list(TRUE)),
               "must name every flag")
  expect_error(frmtmb_structure(loglik = ll_ok, supports = list(osa = "yes")),
               "must be TRUE or FALSE")
  expect_error(
    frmtmb_structure(loglik = ll_ok,
                     supports = list(osa = TRUE, osa = FALSE)),
    "names a flag twice")
  st <- frmtmb_structure(loglik = ll_ok, supports = c(osa = TRUE))
  expect_true(st[["supports"]][["osa"]])
  expect_false(st[["supports"]][["reml"]])
})

test_that("a refusal must explain a flag that is actually refused", {
  expect_error(frmtmb_structure(loglik = ll_ok, refusals = list(zz = "x")),
               "explains no known capability flag")
  expect_error(
    frmtmb_structure(loglik = ll_ok, supports = list(osa = TRUE),
                     refusals = list(osa = "x")),
    "could never be shown")
  expect_error(frmtmb_structure(loglik = ll_ok, refusals = list(osa = 1)),
               "one non-empty string")
  expect_error(frmtmb_structure(loglik = ll_ok, refusals = list("x")),
               "must name the flag")
})

test_that("a loglik that drops the weights the core passes is named", {
  expect_warning(
    frmtmb_structure(loglik = function(y, dpars, aterms, block, extra) 1),
    "weights")
  expect_silent(frmtmb_structure(loglik = ll_ok))
  # a dots-taking loglik has not dropped anything
  expect_silent(frmtmb_structure(loglik = function(y, ...) 1))
})

test_that("a refusal is raised in the family's own words", {
  st <- frmtmb_structure(
    loglik = ll_ok,
    refusals = list(osa = "no one-step residual here",
                    re_form = "no re.form here",
                    re_form.simulate = "not when simulating"))
  expect_error(structure_gate(st, "osa", "GENERIC"),
               "no one-step residual here", fixed = TRUE)
  # a context-specific message wins, and falls back to the bare flag
  expect_error(structure_gate(st, "re_form", "GENERIC",
                              context = "simulate"),
               "not when simulating", fixed = TRUE)
  expect_error(structure_gate(st, "re_form", "GENERIC"),
               "no re.form here", fixed = TRUE)
  # a flag refused without a sentence of its own gets the generic one
  expect_error(structure_gate(st, "deviance", "GENERIC"), "GENERIC",
               fixed = TRUE)
  # nothing is raised for a supported flag or for a rowwise family
  expect_silent(structure_gate(NULL, "osa", "GENERIC"))
})

test_that("frmtmb_family() takes a structure and nothing else", {
  st <- frmtmb_structure(loglik = ll_ok)
  fam <- frmtmb_family("demo", "mu", list(mu = "identity"),
                       lpdf = function(y, dpars, aterms) 0, structure = st)
  expect_identical(fam_structure(fam), st)
  expect_null(fam_structure(as_frmtmb_family(stats::gaussian())))
  expect_error(
    frmtmb_family("demo", "mu", list(mu = "identity"),
                  lpdf = function(y, dpars, aterms) 0,
                  structure = list(loglik = ll_ok)),
    "frmtmb_structure")
})

test_that("a group-level mixture carries the protocol", {
  fg <- mixture(gaussian(), gaussian(), groups = ~g)
  st <- fam_structure(fg)
  expect_s3_class(st, "frmtmb_structure")
  # the four flags that follow from the per-row mean being rowwise
  for (flag in c("conditional_effects", "newdata_response",
                 "se_fit_response", "re_form")) {
    expect_true(st[["supports"]][[flag]], label = flag)
  }
  # the group branch registers no observation vector, so there is
  # nothing for one-step-ahead residuals to step through
  expect_false(st[["supports"]][["osa"]])
  expect_true(is.function(st[["loglik"]]))
})

test_that("a rowwise mixture-type family declares capabilities only", {
  # mixture() and mixture_mvn() have a likelihood that DOES factorize
  # per row, so their structures carry no loglik; what they carry is
  # the multimodality refusal that used to be one gate in fit.R naming
  # every mixture-type family. The third such family lives in
  # frmtmb.latent now and asserts this same row of the table there.
  for (fam in list(mixture(gaussian(), gaussian()),
                   mixture_mvn(K = 2, D = 2))) {
    st <- fam_structure(fam)
    expect_s3_class(st, "frmtmb_structure")
    expect_null(st[["loglik"]])
    expect_false(st[["supports"]][["reml"]])
    expect_false(st[["supports"]][["profile"]])
    # everything else stays exactly as available as it was
    expect_true(st[["supports"]][["quadrature"]])
    expect_true(st[["supports"]][["cluster_robust"]])
    expect_true(is.function(st[["latent_probs"]]))
  }
  expect_true(
    fam_structure(mixture(gaussian(), gaussian()))[["supports"]][["osa"]])
})

test_that("structure_supports_all() is what such a family starts from", {
  # the exported starting point for a capability-only structure, which
  # is what an out-of-tree rowwise family builds on. Every flag TRUE
  # except the named ones, and the flag set is frmtmb's own, so a
  # capability added later arrives switched on rather than silently
  # taken away.
  s <- structure_supports_all(reml = FALSE, osa = FALSE)
  expect_false(s[["reml"]])
  expect_false(s[["osa"]])
  expect_true(s[["profile"]])
  expect_true(s[["quadrature"]])
  expect_true(s[["cluster_robust"]])
  expect_setequal(
    names(s),
    names(fam_structure(mixture(gaussian(), gaussian()))[["supports"]]))
})

test_that("the frame carries one block per structured response", {
  set.seed(4)
  ng <- 12L
  m <- 5L
  g <- rep(seq_len(ng), each = m)
  cls <- stats::rbinom(ng, 1L, 0.5)
  dd <- data.frame(y = stats::rnorm(ng * m, c(-2, 2)[cls + 1L][g], 0.5),
                   g = factor(g))
  fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
             data = dd)
  blk <- fit$frame[["blocks"]][["y"]]
  expect_type(blk, "list")
  expect_identical(length(blk[["first"]]), ng)
  expect_identical(blk[["levels"]], levels(dd$g))
  # the block is data, so a refit rebuilds an identical one
  expect_null(fit$frame[["mix_g"]])

  # a rowwise fit carries no block at all
  f2 <- frm(bf(y ~ 1) + gaussian(), data = dd)
  expect_length(f2$frame[["blocks"]], 0L)
})

# ---------------------------------------------------------------------
# HOW FINELY THE LIKELIHOOD FACTORIZES.
#
# `loglik` is one number for the whole response, and three consumers
# need the pieces: the importance correction resamples a grouping level,
# loo() leaves a unit out, and a deviance residual compares one row with
# its saturated fit. `loglik_row` and `loglik_group` carry the pieces a
# family HAS; `unit` still says what may be left out, which is a
# different question and keeps its own answer.
# ---------------------------------------------------------------------

test_that("a factorization slot needs a likelihood to factorize", {
  # with loglik = NULL the family keeps its own rowwise lpdf, which the
  # core already evaluates a row at a time: a slot here would be a
  # second definition of the same numbers
  expect_error(frmtmb_structure(loglik_row = ll_ok), "loglik_row")
  expect_error(frmtmb_structure(loglik_group = ll_ok),
               "with loglik = NULL", fixed = TRUE)
  expect_error(frmtmb_structure(loglik = ll_ok, loglik_row = "no"),
               "loglik_row")
})

test_that("both factorization slots warn when they drop the weights", {
  no_w <- function(y, dpars, aterms) sum(y)
  expect_warning(frmtmb_structure(loglik = ll_ok, loglik_row = no_w),
                 "loglik_row =) has no `weights`", fixed = TRUE)
  expect_warning(frmtmb_structure(loglik = ll_ok, loglik_group = no_w),
                 "loglik_group =) has no `weights`", fixed = TRUE)
  # and a slot that takes them is silent
  expect_silent(frmtmb_structure(loglik = ll_ok, loglik_row = ll_ok,
                                 loglik_group = ll_ok))
})

test_that("the declared slots are what print() reports", {
  st <- frmtmb_structure(loglik = ll_ok, loglik_group = ll_ok,
                         unit = "one subject's sequence")
  expect_null(st[["loglik_row"]])
  expect_type(st[["loglik_group"]], "closure")
  expect_identical(st[["unit"]], "one subject's sequence")
  expect_output(print(st), "loglik_group")
})

test_that("a block's grouping is read in one fixed order", {
  # a factor is read on its LEVELS and anything else on its sorted
  # values, so a family that numbers its units 1, 2, 10 does not get 10
  # before 2 through a character sort
  expect_null(structure_group_codes(list(n = 3)))
  expect_identical(structure_group_codes(list(group = c(1L, 10L, 2L))),
                   c(1L, 3L, 2L))
  f <- factor(c("b", "a"), levels = c("b", "a"))
  expect_identical(structure_group_codes(list(group = f)), c(1L, 2L))
})

test_that("the block a factorization rests on is checked at assembly", {
  fam <- frmtmb_family("toy", dpars = "mu", links = list(mu = "identity"),
                       lpdf = function(y, dpars, aterms) 0 * y)
  grp <- frmtmb_structure(loglik = ll_ok, loglik_group = ll_ok)
  row <- frmtmb_structure(loglik = ll_ok, loglik_row = ll_ok)
  expect_error(check_structure_block(grp, list(n = 4), fam, "y", 4L),
               "carries no `group`")
  # a per-ROW slot may leave it out; only the correction then refuses
  expect_null(check_structure_block(row, list(n = 4), fam, "y", 4L))
  expect_error(
    check_structure_block(grp, list(group = 1:3), fam, "y", 4L),
    "for 3 rows, and the frame has 4", fixed = TRUE)
  expect_error(
    check_structure_block(grp, list(group = c(1L, NA, 2L, 2L)), fam, "y", 4L),
    "leaves `group` missing on 1 row")
  expect_error(
    check_structure_block(grp, list(group = matrix(1:4, 2)), fam, "y", 4L),
    "one entry per row")
  expect_null(check_structure_block(grp, list(group = factor(letters[c(1, 1, 2, 2)])),
                                    fam, "y", 4L))
  # a family with no structure at all is not asked for any of this
  expect_null(check_structure_block(NULL, NULL, fam, "y", 4L))
})

test_that("a deviance residual needs the saturated half, and says so", {
  # the fitted half alone cannot make a unit deviance: the saturated
  # log-density is zero for a Bernoulli trial and is not for a Poisson
  # count, and nothing the core can see tells the two apart
  st <- frmtmb_structure(
    loglik = ll_ok,
    loglik_row = function(y, dpars, aterms, weights, block, extra) {
      rep(-0.5, length(y))
    },
    supports = list(deviance = TRUE))
  fit <- list(frame = list(y = list(y = c(1, 0, 1)),
                           aterm_values = list(y = list())))
  rspec <- list(resp_name = "y",
                family = list(family = "toy"))
  local_mocked_bindings(eval_dpars = function(fit, ...) list(y = list()),
                        fit_extras = function(fit) NULL)
  expect_error(structure_unit_deviance(fit, rspec, st, list(n = 3)),
               "attr(x, \"saturated\")", fixed = TRUE)
  st2 <- frmtmb_structure(loglik = ll_ok, supports = list(deviance = TRUE))
  expect_error(structure_unit_deviance(fit, rspec, st2, list(n = 3)),
               "frmtmb_structure(loglik_row = )", fixed = TRUE)
})

test_that("a saturated value below the fitted one is refused, not clamped", {
  # clamping the unit deviance at zero would report a residual of zero,
  # which reads as a PERFECT fit at exactly the rows where the family
  # has its two log-densities the wrong way round
  st <- frmtmb_structure(
    loglik = ll_ok,
    loglik_row = function(y, dpars, aterms, weights, block, extra) {
      structure(rep(-0.5, length(y)), saturated = c(0, -2, 0))
    },
    supports = list(deviance = TRUE))
  fit <- list(frame = list(y = list(y = c(1, 0, 1)),
                           aterm_values = list(y = list())))
  rspec <- list(resp_name = "y", family = list(family = "toy"))
  local_mocked_bindings(eval_dpars = function(fit, ...) list(y = list()),
                        fit_extras = function(fit) NULL)
  expect_error(structure_unit_deviance(fit, rspec, st, list(n = 3)),
               "saturated log-density BELOW the fitted one at 1 row")
  # and the honest case still returns the unit deviances
  st2 <- frmtmb_structure(
    loglik = ll_ok,
    loglik_row = function(y, dpars, aterms, weights, block, extra) {
      structure(rep(-0.5, length(y)), saturated = rep(0, 3))
    },
    supports = list(deviance = TRUE))
  expect_equal(structure_unit_deviance(fit, rspec, st2, list(n = 3)),
               rep(1, 3))
})

test_that("the allow-list names every term it refuses", {
  fam <- frmtmb_family("toy", dpars = "mu", links = list(mu = "identity"),
                       accepts_aterms = "weights",
                       lpdf = function(y, dpars, aterms) 0 * y)
  one <- expect_error(check_accepted_aterms(
    list(family = fam, aterms = list(trials = quote(n))), list()))
  expect_match(conditionMessage(one), "term `trials()` is not one",
               fixed = TRUE)
  two <- expect_error(check_accepted_aterms(
    list(family = fam,
         aterms = list(trials = quote(n), vint1 = quote(z))), list()))
  expect_match(conditionMessage(two), "terms `trials()`, `vint()` are not",
               fixed = TRUE)
  expect_match(conditionMessage(two), "This family takes `weights()`.",
               fixed = TRUE)
})

# ---------------------------------------------------------------------
# WHO APPLIES THE ROW WEIGHTS.
#
# One answer, for every slot and every consumer: the FAMILY does. The
# core passes them in and never multiplies the result again. The
# deviance path used to do both - pass them into loglik_row AND multiply
# the returned unit deviance by them - so a family following the
# documented instruction reported residuals exactly sqrt(w) too large.
# ---------------------------------------------------------------------

# gaussian, written as a structure, with the weights applied inside the
# slot as the documentation instructs. `saturated` carries the same
# weight, so the unit deviance is weighted once.
struct_wgauss <- function() {
  ld <- function(y, dpars) {
    stats::dnorm(y, dpars[["mu"]], dpars[["sigma"]], log = TRUE)
  }
  frmtmb_family(
    "wgauss", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    accepts_aterms = "weights",
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars[["mu"]], dpars[["sigma"]], log = TRUE)
    },
    init_dpars = list(mu = function(y, aterms) mean(y),
                      sigma = function(y, aterms) stats::sd(y)),
    structure = frmtmb_structure(
      loglik = function(y, dpars, aterms, weights, block, extra) {
        sum(weights * RTMB::dnorm(y, dpars[["mu"]], dpars[["sigma"]],
                                  log = TRUE))
      },
      loglik_row = function(y, dpars, aterms, weights, block, extra) {
        out <- weights * ld(y, dpars)
        # the saturated fit puts mu at y; it carries the weight too
        attr(out, "saturated") <- weights *
          ld(y, list(mu = y, sigma = dpars[["sigma"]]))
        out
      },
      # the structured deviance path runs inside the fitted_mean
      # branch of residuals(), so a family wanting it declares one
      fitted_mean = function(fit, block) {
        as.numeric(eval_dpars(fit)[["y"]][["mu"]])
      },
      supports = list(deviance = TRUE, conditional_effects = TRUE)))
}

test_that("the family applies the row weights, and the core does not", {
  skip_on_cran()
  set.seed(808)
  dd <- data.frame(x = rnorm(60))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
  dd$w <- rep(c(1, 2, 4, 3), length.out = 60)
  wf <- bf(y | weights(w) ~ x)
  fit <- frm(wf, data = dd, family = struct_wgauss())
  ref <- frm(wf + gaussian(), data = dd)
  # the same weighted likelihood, so the same fit
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-8)
  expect_equal(unlist(fixef(fit)), unlist(fixef(ref)), tolerance = 1e-6)

  d_str <- as.numeric(residuals(fit, type = "deviance"))
  d_ref <- as.numeric(residuals(ref, type = "deviance"))
  # The toy's unit deviance is the SCALED one (its log-densities carry
  # sigma), so it differs from gaussian()'s by a constant 1 / sigma. A
  # constant is all it may differ by: the ratio must not track the
  # weights, which is what the double application made it do.
  ratio <- d_str / d_ref
  expect_equal(max(ratio) - min(ratio), 0, tolerance = 1e-6)
  sigma <- as.numeric(eval_dpars(fit)[["y"]][["sigma"]])[1L]
  expect_equal(mean(ratio), 1 / sigma, tolerance = 1e-6)
  # That constancy IS the regression pin. Under the double
  # application the ratio was sqrt(w) / sigma, which over weights of
  # 1, 2, 3 and 4 spreads by a factor of two.
  spread <- max(sqrt(dd$w)) / min(sqrt(dd$w))
  expect_gt(spread, 1.9)
})

test_that("an unused level in the block's grouping is refused by name", {
  fam <- frmtmb_family("toy", dpars = "mu", links = list(mu = "identity"),
                       lpdf = function(y, dpars, aterms) 0 * y)
  grp <- frmtmb_structure(loglik = ll_ok, loglik_group = ll_ok)
  g <- factor(c("a", "a", "b", "b"), levels = c("a", "MID", "b"))
  err <- expect_error(check_structure_block(grp, list(group = g), fam,
                                            "y", 4L))
  expect_match(conditionMessage(err), "1 unused level(s) ('MID')",
               fixed = TRUE)
  expect_match(conditionMessage(err), "toy", fixed = TRUE)
  # and the codes the core works in are contiguous whatever arrives, so
  # the permutation the correction builds cannot carry an NA
  expect_identical(structure_group_codes(list(group = g)),
                   c(1L, 1L, 2L, 2L))
})

test_that("a per-row deviance without a fitted mean names the missing half", {
  # the magnitude is in loglik_row(); the sign needs fitted_mean(). The
  # rowwise path used to answer that the family has no unit deviance
  ld <- function(y, dpars) {
    stats::dnorm(y, dpars[["mu"]], dpars[["sigma"]], log = TRUE)
  }
  fam <- frmtmb_family(
    "nomeangauss", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars[["mu"]], dpars[["sigma"]], log = TRUE)
    },
    init_dpars = list(mu = function(y, aterms) mean(y),
                      sigma = function(y, aterms) stats::sd(y)),
    structure = frmtmb_structure(
      loglik = function(y, dpars, aterms, weights, block, extra) {
        sum(weights * RTMB::dnorm(y, dpars[["mu"]], dpars[["sigma"]],
                                  log = TRUE))
      },
      loglik_row = function(y, dpars, aterms, weights, block, extra) {
        out <- weights * ld(y, dpars)
        attr(out, "saturated") <- weights *
          ld(y, list(mu = y, sigma = dpars[["sigma"]]))
        out
      },
      supports = list(deviance = TRUE)))
  set.seed(3)
  d <- data.frame(y = rnorm(40, 2, 1.5))
  fit <- frm(bf(y ~ 1), family = fam, data = d)
  expect_error(residuals(fit, type = "deviance"),
               "declares no fitted_mean")
})
