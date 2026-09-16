# The two-dialect argument seam on the DRAWS side: a brms-named method
# takes re_formula and accepts lme4's re.form as an alias of it, and
# giving both at once is refused rather than resolved.
#
# frmtmb keeps the fit-side half of this suite (pp_check() on a fit,
# and the rule that frmtmb's own fit surface speaks lme4's spelling
# alone). The blocks here are the draws-side half.

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

test_that("every brms-named draws method carries both spellings", {
  # a structural guard, so a method added later to this family cannot
  # quietly ship one spelling: the list is the contract
  dual <- c("posterior_epred.frmtmb_draws", "posterior_linpred.frmtmb_draws",
            "posterior_predict.frmtmb_draws",
            "predictive_interval.frmtmb_draws",
            "predictive_error.frmtmb_draws",
            "pp_check.frmtmb_draws")
  for (nm in dual) {
    fo <- formals(getFromNamespace(nm, "frmtmb.sample"))
    expect_true(all(c("re_formula", "re.form") %in% names(fo)),
                label = paste0(nm, " has both spellings"))
    # both default to the "not supplied" marker, so either alone is a
    # setting and neither NULL nor NA is mistaken for one
    ns <- asNamespace("frmtmb.sample")
    expect_true(frmtmb:::is_arg_unset(eval(fo[["re_formula"]], ns)))
    expect_true(frmtmb:::is_arg_unset(eval(fo[["re.form"]], ns)))
  }
})

test_that("re_formula and re.form give identical pp_check() output", {
  skip_if_not_installed("bayesplot")
  cs <- sp_case()

  set.seed(7)
  p1 <- pp_check(cs$fit, ndraws = 5, re_formula = NA)
  set.seed(7)
  p2 <- pp_check(cs$fit, ndraws = 5, re.form = NA)
  expect_equal(p1$data, p2$data)

  set.seed(8)
  d1 <- pp_check(cs$ds, ndraws = 5, re_formula = NA)
  set.seed(8)
  d2 <- pp_check(cs$ds, ndraws = 5, re.form = NA)
  expect_equal(d1$data, d2$data)

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
  expect_error(predictive_interval(ds, re_formula = NA, re.form = NA),
               "predictive_interval\\(\\)")
  expect_error(predictive_error(ds, re_formula = NA, re.form = NA),
               "predictive_error\\(\\)")
  expect_error(pp_check(ds, re_formula = NA, re.form = NA),
               "pp_check\\(\\)")
  expect_error(pp_check(cs$fit, re_formula = NA, re.form = NA),
               "pp_check\\(\\)")

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
  expect_equal(rownames(posterior_interval(ds, "^x$")), "x")
  expect_equal(colnames(as.mcmc(ds, "^x$")[[1L]]), "x")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("bayesplot")
  expect_s3_class(mcmc_plot(ds, "^x$"), "ggplot")
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
