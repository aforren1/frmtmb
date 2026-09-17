# The one place two spellings survive, and the boundary of it. brms is
# the tiebreaker on every argument name, so lme4's re.form was dropped
# from the fit surface. It stays on the FIVE draws methods where brms
# ACCEPTS it, which is not the same as the four that DECLARE it:
# predictive_interval.brmsfit() declares neither spelling and its whole
# body is posterior_predict(object, ...), so the alias reaches a formal
# one frame down. What brms accepts is the test; what it declares is
# not. dev/argspell-brms-accepts.R measures both.
#
# pp_check() is the one method here that does NOT take the alias, and
# that is now settled rather than open. brms HONORS re.form on
# pp_check(), through the same dots forwarding, while warning
# "unrecognized and ignored" as it does so. Warn-then-honor is the
# failure this item exists to stop, so the spelling is refused here
# rather than copied: matching brms means matching what brms decided,
# not reproducing a leak. predictive_interval() is the contrast, where
# brms honors the alias silently and deliberately and this package
# follows it. dev/argspell-brms-accepts.R records the mechanism.
#
# frmtmb keeps the fit-side half of this suite (pp_check() on a fit, and
# the rule that the fit surface speaks brms alone).

# ---- from tests/testthat/test-api-spellings.R ----

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

sp_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(11)
      dd <- data.frame(x = stats::rnorm(60),
                       g = factor(rep(1:6, 10)))
      dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x +
                             stats::rnorm(6, 0, 0.5)[dd$g], 1)
      fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 3)))
      cache <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache
  }
})

test_that("re_formula and re.form give identical draws-method output", {
  cs <- sp_case()
  ds <- cs$ds
  nd <- data.frame(x = c(-1, 0, 1),
                   g = factor(1, levels = levels(cs$dd$g)))

  for (rf in list(NA, NULL, ~0)) {
    expect_equal(posterior_epred(ds, newdata = nd, re_formula = rf,
                                 ndraws = 12),
                 posterior_epred(ds, newdata = nd, re.form = rf,
                                 ndraws = 12))
    expect_equal(posterior_linpred(ds, newdata = nd, re_formula = rf,
                                   ndraws = 12),
                 posterior_linpred(ds, newdata = nd, re.form = rf,
                                   ndraws = 12))
  }

  # the simulating methods need the same seed to be comparable at all
  set.seed(4)
  a <- posterior_predict(ds, newdata = nd, re_formula = NA, ndraws = 12)
  set.seed(4)
  b <- posterior_predict(ds, newdata = nd, re.form = NA, ndraws = 12)
  expect_equal(a, b)

  set.seed(5)
  a <- predictive_interval(ds, re_formula = NA, ndraws = 12)
  set.seed(5)
  b <- predictive_interval(ds, re.form = NA, ndraws = 12)
  expect_equal(a, b)

  set.seed(6)
  a <- predictive_error(ds, re_formula = NA, ndraws = 12)
  set.seed(6)
  b <- predictive_error(ds, re.form = NA, ndraws = 12)
  expect_equal(a, b)
})

test_that("the draws methods still default to NULL", {
  cs <- sp_case()
  ds <- cs$ds

  # the pinned pre-change default of every draws method: condition on
  # each draw's own random effects
  expect_equal(posterior_epred(ds, ndraws = 12),
               posterior_epred(ds, re.form = NULL, ndraws = 12))
  expect_equal(posterior_epred(ds, ndraws = 12),
               posterior_epred(ds, re_formula = NULL, ndraws = 12))
  expect_equal(posterior_linpred(ds, ndraws = 12),
               posterior_linpred(ds, re.form = NULL, ndraws = 12))

  # and NULL is a different quantity from NA, so the assertion above is
  # not vacuous: a default silently flipped to NA would show up here
  expect_false(isTRUE(all.equal(
    posterior_epred(ds, ndraws = 12),
    posterior_epred(ds, re_formula = NA, ndraws = 12))))

  set.seed(9)
  a <- predictive_error(ds, ndraws = 12)
  set.seed(9)
  b <- predictive_error(ds, re.form = NULL, ndraws = 12)
  expect_equal(a, b)
})

# Does brms ACCEPT the alias on this generic? Declared as a formal, or
# reached one frame down by a body that forwards its dots to a brmsfit
# method that declares it. Comparing DECLARED formals alone got
# predictive_interval() wrong, because brms declares neither spelling
# there and honors both.
ar_brms_accepts <- function(gen) {
  bm <- tryCatch(getFromNamespace(paste0(gen, ".brmsfit"), "brms"),
                 error = function(e) NULL)
  if (is.null(bm)) return(NA)
  if ("re.form" %in% names(formals(bm))) return(TRUE)
  b <- deparse(body(bm))
  if (!any(grepl("...", b, fixed = TRUE))) return(FALSE)
  targets <- c("posterior_predict", "posterior_epred",
               "posterior_linpred")
  hit <- targets[vapply(targets,
                        function(t) any(grepl(t, b, fixed = TRUE)), TRUE)]
  any(vapply(hit, function(t) {
    m <- tryCatch(getFromNamespace(paste0(t, ".brmsfit"), "brms"),
                  error = function(e) NULL)
    !is.null(m) && "re.form" %in% names(formals(m))
  }, TRUE))
}

test_that("the alias detector reads brms's behavior, not its formals", {
  skip_if_not_installed("brms")
  # constructed in both cases before it is used on anything: a method
  # that declares it, one that only forwards to it, and one that does
  # neither
  expect_true(ar_brms_accepts("posterior_predict"))
  expect_true(ar_brms_accepts("predictive_interval"))
  expect_false("re.form" %in%
                 names(formals(getFromNamespace("predictive_interval.brmsfit",
                                                "brms"))))
  expect_true(is.na(ar_brms_accepts("no_such_generic")))
})

test_that("the dual spelling covers the methods brms accepts it on", {
  # the contract, in the form that can go wrong both ways: a method that
  # loses an alias brms accepts, and a method that keeps one brms does
  # not.
  dual <- c("posterior_epred.frmtmb_draws", "posterior_linpred.frmtmb_draws",
            "posterior_predict.frmtmb_draws",
            "predictive_error.frmtmb_draws",
            "predictive_interval.frmtmb_draws")
  ns <- asNamespace("frmtmb.sample")
  for (nm in dual) {
    fo <- formals(getFromNamespace(nm, "frmtmb.sample"))
    expect_true(all(c("re_formula", "re.form") %in% names(fo)),
                label = paste0(nm, " has both spellings"))
    # both default to the "not supplied" marker, so either alone is a
    # setting and neither NULL nor NA is mistaken for one
    expect_true(frmtmb:::is_arg_unset(eval(fo[["re_formula"]], ns)))
    expect_true(frmtmb:::is_arg_unset(eval(fo[["re.form"]], ns)))
  }
  # pp_check() does not take it, deliberately: brms honors it there and
  # warns that it did not, and this package refuses rather than copy a
  # warn-then-honor path. `ar_brms_accepts()` reports TRUE for pp_check,
  # which is why this is asserted separately from the loop below rather
  # than derived from it.
  fo <- formals(getFromNamespace("pp_check.frmtmb_draws", "frmtmb.sample"))
  expect_false("re.form" %in% names(fo))

  skip_if_not_installed("brms")
  for (nm in dual) {
    gen <- sub("[.]frmtmb_draws$", "", nm)
    expect_true(isTRUE(ar_brms_accepts(gen)),
                label = paste("brms accepts re.form on", gen))
  }
})

test_that("the re_formula switch takes effect in sample", {
  skip_if_not_installed("bayesplot")
  cs <- sp_case()

  # non-vacuity: the switch takes effect IN-SAMPLE (review finding: it
  # used to be consulted only under newdata), so NA must differ from
  # the default NULL wherever the model has random effects, under the
  # same RNG seed
  set.seed(9)
  dNA <- posterior_predict(cs$ds, ndraws = 5, re_formula = NA)
  set.seed(9)
  dNU <- posterior_predict(cs$ds, ndraws = 5)
  expect_false(isTRUE(all.equal(dNA, dNU)))
  set.seed(9)
  dNA2 <- posterior_predict(cs$ds, ndraws = 5, re.form = NA)
  expect_equal(dNA, dNA2)
  set.seed(9)
  eNA <- predictive_error(cs$ds, ndraws = 5, re_formula = NA)
  set.seed(9)
  eNU <- predictive_error(cs$ds, ndraws = 5)
  expect_false(isTRUE(all.equal(eNA, eNU)))
  set.seed(9)
  pNA <- pp_check(cs$ds, ndraws = 5, re_formula = NA)
  set.seed(9)
  pNU <- pp_check(cs$ds, ndraws = 5)
  expect_false(isTRUE(all.equal(pNA$data, pNU$data)))
})

test_that("giving both spellings is refused, not resolved", {
  cs <- sp_case()
  ds <- cs$ds

  expect_error(posterior_epred(ds, re_formula = NA, re.form = NA),
               "two spellings of ONE setting")
  expect_error(posterior_linpred(ds, re_formula = NA, re.form = NA),
               "posterior_linpred\\(\\)")
  expect_error(posterior_predict(ds, re_formula = NA, re.form = NA),
               "posterior_predict\\(\\)")
  expect_error(predictive_error(ds, re_formula = NA, re.form = NA),
               "predictive_error\\(\\)")
  expect_error(predictive_interval(ds, re_formula = NA, re.form = NA),
               "predictive_interval\\(\\)")
  # pp_check() no longer takes the alias, and refuses it by name rather
  # than reporting a clash of spellings
  expect_error(pp_check(ds, re.form = NA), "re.form", fixed = TRUE)

  # the refusal names both spellings and says which one the function is
  # named after, so it can be acted on without reading the manual
  msg <- tryCatch(posterior_epred(ds, re_formula = NA, re.form = NA),
                  error = conditionMessage)
  expect_match(msg, "`re_formula`")
  expect_match(msg, "`re.form`")

  # agreeing values are refused too: the point is that the call did not
  # say which name it meant, not that the two disagreed
  expect_error(posterior_epred(ds, re_formula = NULL, re.form = NULL),
               "two spellings of ONE setting")
})

## ---- brms's positional slots ----------------------------------------
#
# The other half of the seam. A brms-named method speaks brms's argument
# NAMES, and it must also put them in brms's ORDER, or a positional call
# ported from brms answers a different question with nothing said.
# dev/samplegen-findings.md item 4 found ten methods diverging at a
# positional slot, two of them silently.

brms_ten <- c("as.mcmc", "log_lik", "mcmc_plot", "posterior_epred",
              "posterior_interval", "posterior_linpred",
              "posterior_predict", "pp_mixture", "predictive_error",
              "psis")

# the arguments a caller can reach positionally: everything before the
# method's own `...`
positional_args <- function(f) {
  a <- names(formals(f))
  a[seq_len(match("...", a, nomatch = length(a) + 1L) - 1L)]
}

# the review's own criterion: the first position at which the two
# names differ, over the positions both sides have
first_divergence <- function(b, o) {
  k <- min(length(b), length(o))
  if (!k) return(NA_integer_)
  d <- which(b[seq_len(k)] != o[seq_len(k)])
  if (length(d)) d[[1L]] else NA_integer_
}

test_that("ten brms-facing methods take brms's arguments in brms's order", {
  # generated from the installed brms rather than typed, so it cannot
  # drift from the package it is matching
  skip_if_not_installed("brms")
  tb <- get(".__S3MethodsTable__.", envir = asNamespace("frmtmb.sample"),
            inherits = FALSE)
  bad <- character()
  for (nm in brms_ten) {
    b <- positional_args(get(paste0(nm, ".brmsfit"),
                             envir = asNamespace("brms")))
    o <- positional_args(get(paste0(nm, ".frmtmb_draws"), envir = tb))
    p <- first_divergence(b, o)
    if (!is.na(p)) {
      bad <- c(bad, sprintf("%s: position %d is brms's `%s` and this ",
                            nm, p, b[[p]]))
    }
  }
  expect_equal(bad, character())

  # the guard is only a guard if it read ten methods and can name a
  # wrong one. The vector below is the signature log_lik() shipped with,
  # and its second position was ndraws where brms's is newdata
  expect_length(brms_ten, 10L)
  expect_equal(
    first_divergence(
      positional_args(get("log_lik.brmsfit", envir = asNamespace("brms"))),
      c("object", "ndraws", "resp")),
    2L)
})

test_that("brms's positional calls mean here what they mean in brms", {
  cs <- sp_case()
  ds <- cs$ds
  nd <- data.frame(x = c(-1, 0, 1),
                   g = factor(1, levels = levels(cs$dd$g)))

  # position 3 of the predictive methods is re_formula, so brms's
  # `re_formula = NA` idiom drops the random effects rather than naming
  # a response
  expect_equal(posterior_epred(ds, nd, NA),
               posterior_epred(ds, newdata = nd, re_formula = NA))
  expect_false(isTRUE(all.equal(posterior_epred(ds, nd, NA),
                                posterior_epred(ds, newdata = nd))))
  expect_equal(posterior_linpred(ds, FALSE, nd, NA),
               posterior_linpred(ds, newdata = nd, re_formula = NA))
  set.seed(4)
  a <- posterior_predict(ds, nd, NA)
  set.seed(4)
  b <- posterior_predict(ds, newdata = nd, re_formula = NA)
  expect_equal(a, b)

  # the two that used to ANSWER a different question. brms's answer to
  # both is a refusal, because that slot is `pars` and brms takes only
  # NA or a character vector there
  expect_error(as.mcmc(ds, TRUE), "must be NA or a character vector")
  expect_error(posterior_interval(ds, 0.9),
               "must be NA or a character vector")
  # and the call brms does answer in that slot works
  # posterior_interval() reaches brms's as.matrix(pars =), which warns
  # that `pars` is deprecated, as brms's does
  expect_warning(pi <- posterior_interval(ds, "^b_x$"), "deprecated")
  expect_equal(rownames(pi), "b_x")
  expect_equal(colnames(as.mcmc(ds, "^b_x$")[[1L]]), "b_x")
  # the six whose slots used to diverge: bayes_R2's fourth is robust,
  # posterior_summary's second is pars, hypothesis's third is class
  expect_equal(bayes_R2(ds, NULL, TRUE, TRUE),
               bayes_R2(ds, robust = TRUE))
  expect_equal(colnames(bayes_R2(ds, NULL, TRUE, TRUE)),
               c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_equal(rownames(suppressWarnings(posterior_summary(ds, "^b_x$"))),
               "b_x")
  expect_equal(hypothesis(ds, "Intercept > 0", "sd", "g")$class, "sd_g")
  # and the summary = FALSE family returns the draws, as in brms
  expect_equal(dim(fixef(ds, FALSE)), c(ndraws(ds), 2L))
  expect_equal(length(dim(ranef(ds, FALSE)$g)), 3L)
  expect_true(is.matrix(VarCorr(ds, NULL, FALSE)$g$sd))
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("bayesplot")
  expect_s3_class(mcmc_plot(ds, "^b_x$"), "ggplot")
})

test_that("the slots this package cannot answer refuse by name", {
  cs <- sp_case()
  ds <- cs$ds
  nd <- data.frame(x = c(-1, 0, 1),
                   g = factor(1, levels = levels(cs$dd$g)))
  expect_error(log_lik(ds, nd), "does not take newdata")
  expect_error(log_lik(ds, re_formula = NA), "does not take re_formula")
  expect_error(psis(ds, nd), "does not take newdata")
  # brms's own default, spelled out positionally, is the no-op these
  # methods already do and is NOT refused
  expect_silent(ll <- log_lik(ds, NULL, NULL, NULL, 5L))
  expect_equal(nrow(ll), 5L)
  expect_error(as.mcmc(ds, NA, FALSE, FALSE, TRUE),
               "nothing to include")
  expect_error(posterior_predict(ds, negative_rt = TRUE), "wiener")
  expect_error(posterior_linpred(ds, incl_thres = TRUE), "thresholds")
})

test_that("predictive_error takes brms's newdata and method slots", {
  cs <- sp_case()
  ds <- cs$ds
  ndy <- data.frame(x = c(-1, 0, 1),
                    g = factor(1, levels = levels(cs$dd$g)),
                    y = c(0.2, 0.4, 0.6))
  pe <- predictive_error(ds, ndy, NA, method = "posterior_epred")
  expect_equal(pe, sweep(-posterior_epred(ds, newdata = ndy,
                                          re_formula = NA),
                         2L, ndy$y, "+"))
  # the response must be there to subtract, and the refusal says so
  expect_error(predictive_error(ds, ndy[, c("x", "g")]),
               "needs the observed response")
})

test_that("draw_ids names the draws, in brms's own position", {
  cs <- sp_case()
  ds <- cs$ds
  ids <- c(2L, 5L, 11L)
  ep <- posterior_epred(ds, draw_ids = ids)
  expect_equal(nrow(ep), 3L)
  expect_equal(ep, posterior_epred(ds)[ids, , drop = FALSE])
  expect_equal(log_lik(ds, draw_ids = ids), log_lik(ds)[ids, ,
                                                        drop = FALSE])
  expect_error(posterior_epred(ds, ndraws = 3, draw_ids = ids),
               "only one of them")
  expect_error(posterior_epred(ds, draw_ids = c(0L, 1L)),
               "whole numbers between 1 and")
})
