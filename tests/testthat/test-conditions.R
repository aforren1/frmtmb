# Classed conditions, brms contract 1 (?frmtmb-conditions).
#
# brms raises every refusal as a `brms_error`, so `tryCatch(brms_error
# = )` catches any of them. frmtmb keeps the same contract with
# `frmtmb_error`, `frmtmb_warning` and `frmtmb_message`. The census in
# test-conditions-census.R keeps every call site on the helpers; this
# file asserts what the helpers deliver, and that the switch from the
# base functions changed the class of a condition and nothing else.

test_that("a refusal from frm() is caught by class", {
  d <- data.frame(y = c(1, 2, 3))
  e <- tryCatch(frm(y ~ 1, data = d, family = "not_a_family"),
                frmtmb_error = function(e) e)
  expect_s3_class(e, "frmtmb_error")
  expect_identical(class(e), c("frmtmb_error", "error", "condition"))
  expect_match(conditionMessage(e), "is not a supported family")
})

test_that("the class vectors match the brms shape", {
  e <- tryCatch(frm_stop("x"), error = identity)
  w <- tryCatch(frm_warning("x"), warning = identity)
  m <- tryCatch(frm_message("x"), message = identity)
  expect_identical(class(e), c("frmtmb_error", "error", "condition"))
  expect_identical(class(w), c("frmtmb_warning", "warning", "condition"))
  expect_identical(class(m), c("frmtmb_message", "message", "condition"))
  # a caller-supplied class goes in front, as brms's .subclass does
  e <- tryCatch(frm_stop("x", class = "frmtmb_fit_error"), error = identity)
  expect_identical(class(e), c("frmtmb_fit_error", "frmtmb_error", "error",
                               "condition"))
})

test_that("a condition raised in another namespace gains its subclass", {
  # an extension is found from where the raising function was defined;
  # stats stands in for an extension namespace, and a function defined
  # outside any namespace gets frmtmb's classes only
  f <- function() frm_stop("x")
  environment(f) <- asNamespace("stats")
  e <- tryCatch(f(), error = identity)
  expect_identical(class(e), c("stats_error", "frmtmb_error", "error",
                               "condition"))
  g <- function() frm_warning("x")
  environment(g) <- asNamespace("stats")
  w <- tryCatch(g(), warning = identity)
  expect_identical(class(w), c("stats_warning", "frmtmb_warning",
                               "warning", "condition"))
  h <- function() frm_stop("x")
  environment(h) <- globalenv()
  expect_identical(class(tryCatch(h(), error = identity)),
                   c("frmtmb_error", "error", "condition"))
})

# The base function and the helper, raised from the same shape of
# caller, must agree on message, call and printed text.

test_that("frm_stop() records the message and call stop() records", {
  # several pieces, a number, a newline and a vector, which is what
  # .makeMessage() has to paste the same way
  parts_a <- function() stop("a ", 1L, " b\n", c("c", "d"), call. = FALSE)
  parts_b <- function() frm_stop("a ", 1L, " b\n", c("c", "d"),
                                 call. = FALSE)
  expect_identical(conditionMessage(tryCatch(parts_b(), error = identity)),
                   conditionMessage(tryCatch(parts_a(), error = identity)))
  # the ordinary shape of a call site: a named function calling it
  # directly, with and without call.
  site_a <- function(x) stop("bad x: ", x)
  site_b <- function(x) frm_stop("bad x: ", x)
  ea <- tryCatch(site_a(3), error = identity)
  eb <- tryCatch(site_b(3), error = identity)
  expect_identical(conditionMessage(eb), conditionMessage(ea))
  expect_identical(conditionCall(eb), quote(site_b(3)))
  expect_identical(conditionCall(ea), quote(site_a(3)))
  # the text a user sees printed, which try() renders from both
  ta <- try(site_a(3), silent = TRUE)
  tb <- try(site_b(3), silent = TRUE)
  expect_identical(sub("site_b", "site_a", as.character(tb)),
                   as.character(ta))
  quiet_a <- function() stop("quiet", call. = FALSE)
  quiet_b <- function() frm_stop("quiet", call. = FALSE)
  expect_identical(as.character(try(quiet_b(), silent = TRUE)),
                   as.character(try(quiet_a(), silent = TRUE)))
  # the absent case: an unclassed stop() is not caught by class
  expect_error(tryCatch(site_a(3), frmtmb_error = function(e) NULL),
               "bad x: 3")
})

test_that("frm_stop() records the call when forced inside another frame", {
  # a promise forced deeper in the stack is where sys.call() and stop()
  # could part ways, so it is checked rather than assumed
  force_it <- function(x) x
  pa <- function() force_it(stop("p"))
  pb <- function() force_it(frm_stop("p"))
  ea <- tryCatch(pa(), error = identity)
  eb <- tryCatch(pb(), error = identity)
  # the same frame, whose call text differs only by the function name
  expect_identical(conditionCall(ea), quote(force_it(stop("p"))))
  expect_identical(conditionCall(eb), quote(force_it(frm_stop("p"))))
})

test_that("frm_warning() records the message and call warning() does", {
  wa <- function(x) warning("careful: ", x)
  wb <- function(x) frm_warning("careful: ", x)
  ca <- tryCatch(wa(2), warning = identity)
  cb <- tryCatch(wb(2), warning = identity)
  expect_identical(conditionMessage(cb), conditionMessage(ca))
  expect_identical(conditionCall(ca), quote(wa(2)))
  expect_identical(conditionCall(cb), quote(wb(2)))
  qa <- function() warning("q", call. = FALSE)
  qb <- function() frm_warning("q", call. = FALSE)
  expect_null(conditionCall(tryCatch(qb(), warning = identity)))
  expect_null(conditionCall(tryCatch(qa(), warning = identity)))
  # muffling and the return value behave as warning()'s do
  expect_identical(withCallingHandlers(wb(1), warning = function(w) {
    invokeRestart("muffleWarning")
  }), "careful: 1")
  expect_warning(wb(1), "careful: 1", class = "frmtmb_warning")
  expect_no_warning(suppressWarnings(wb(1), classes = "frmtmb_warning"))
})

test_that("frm_message() prints what message() prints", {
  ma <- function() message("note ", 1, appendLF = TRUE)
  mb <- function() frm_message("note ", 1, appendLF = TRUE)
  expect_identical(capture.output(mb(), type = "message"),
                   capture.output(ma(), type = "message"))
  na <- function() message("no newline", appendLF = FALSE)
  nb <- function() frm_message("no newline", appendLF = FALSE)
  expect_identical(conditionMessage(tryCatch(nb(), message = identity)),
                   conditionMessage(tryCatch(na(), message = identity)))
  expect_message(mb(), "note 1", class = "frmtmb_message")
  expect_silent(suppressMessages(mb()))
})

# Punch round 1 of item 2.6e: the cases no call site reached, where the
# helpers had parted from the base functions.

test_that("a condition object is signalled as the base functions do", {
  e0 <- simpleError("boom", call = quote(f(1)))
  sa <- function() stop(e0)
  sb <- function() frm_stop(e0)
  ea <- tryCatch(sa(), error = identity)
  eb <- tryCatch(sb(), error = identity)
  expect_identical(conditionMessage(eb), conditionMessage(ea))
  expect_identical(conditionCall(eb), conditionCall(ea))
  expect_identical(class(eb), c("frmtmb_error", class(e0)))
  expect_identical(as.character(try(sb(), silent = TRUE)),
                   as.character(try(sa(), silent = TRUE)))
  w0 <- simpleWarning("careful")
  wa <- tryCatch(warning(w0), warning = identity)
  wb <- tryCatch(frm_warning(w0), warning = identity)
  expect_identical(conditionMessage(wb), conditionMessage(wa))
  expect_null(conditionCall(wb))
  expect_identical(class(wb), c("frmtmb_warning", class(w0)))
  m0 <- simpleMessage("note\n")
  mb <- tryCatch(frm_message(m0), message = identity)
  expect_identical(conditionMessage(mb), "note\n")
  expect_identical(class(mb), c("frmtmb_message", class(m0)))
  # a class the object already has is not added twice
  e1 <- errorCondition("x", class = "frmtmb_error")
  expect_identical(class(tryCatch(frm_stop(e1), error = identity)),
                   c("frmtmb_error", "error", "condition"))
})

test_that("a helper passed as a value gets its caller's subclass", {
  # lapply() calls FUN from its own frame, which belongs to base; the
  # refusal belongs to the function that handed the helper over
  f <- function() lapply("x", frm_stop)
  environment(f) <- asNamespace("stats")
  expect_identical(class(tryCatch(f(), error = identity)),
                   c("stats_error", "frmtmb_error", "error", "condition"))
  g <- function() Map(frm_warning, "x")
  environment(g) <- asNamespace("stats")
  expect_identical(class(tryCatch(g(), warning = identity)),
                   c("stats_warning", "frmtmb_warning", "warning",
                     "condition"))
  h <- function() lapply("x", frm_stop)
  environment(h) <- globalenv()
  expect_identical(class(tryCatch(h(), error = identity)),
                   c("frmtmb_error", "error", "condition"))
})

test_that("package = names the subclass in place of the caller's", {
  f <- function() frm_stop("x", call. = FALSE, package = "frmtmb.eam")
  expect_identical(class(tryCatch(f(), error = identity)),
                   c("frmtmb_eam_error", "frmtmb_error", "error",
                     "condition"))
  g <- function() frm_warning("x", package = "frmtmb")
  environment(g) <- asNamespace("stats")
  expect_identical(class(tryCatch(g(), warning = identity)),
                   c("frmtmb_warning", "warning", "condition"))
})

test_that("immediate. = TRUE prints the warning at once", {
  skip_on_cran()
  # testthat muffles every warning a test raises, so the default handler
  # that does the printing only runs in a separate R process
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "suppressMessages(library(frmtmb))",
    "g <- function() warning(\"imm\", immediate. = TRUE)",
    "f <- function() frm_warning(\"imm\", immediate. = TRUE)",
    "g(); message(\"after base\")",
    "f(); message(\"after frm\")"), script)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                  shQuote(script), stdout = TRUE,
                                  stderr = TRUE))
  expect_identical(out, c("Warning in g() : imm", "after base",
                          "Warning in f() : imm", "after frm"))
})

test_that("frm_match_arg() matches as match.arg() does", {
  f <- function(type = c("response", "link", "linear")) frm_match_arg(type)
  expect_identical(f(), "response")
  expect_identical(f("link"), match.arg("link", c("response", "link",
                                                  "linear")))
  # partial matching, which brms relies on through match.arg()
  expect_identical(f("resp"), "response")
  expect_identical(f(NULL), "response")
  g <- function(x) frm_match_arg(x, c("a", "b"), several.ok = TRUE)
  expect_identical(g(c("b", "a")), c("b", "a"))
})

test_that("frm_match_arg() refuses by name, with the value and choices", {
  f <- function(type = c("response", "link")) frm_match_arg(type)
  e <- tryCatch(f("bogus"), error = identity)
  expect_identical(class(e), c("frmtmb_error", "error", "condition"))
  expect_identical(conditionMessage(e),
                   paste0("`type` must be one of \"response\", \"link\", ",
                          "not character \"bogus\""))
  expect_identical(conditionCall(e), quote(f("bogus")))
  expect_error(f(1), "`type` must be one of .*not numeric \"1\"",
               class = "frmtmb_error")
  expect_error(f(c("response", "x")), "not a character of length 2",
               class = "frmtmb_error")
  # the refusal is the calling function's, so an extension's carries
  # the extension's subclass
  environment(f) <- asNamespace("stats")
  expect_s3_class(tryCatch(f("bogus"), error = identity), "stats_error")
})

test_that("user-reachable refusals once left to base R are frmtmb_error", {
  # review of item 2.6e: every one of these was a simpleError from
  # match.arg(), stopifnot() or stats::model.frame(). Each regex is the
  # part of the new message that names what was wrong.
  set.seed(20260917)
  d <- data.frame(y = rnorm(40), x = rnorm(40), g = gl(8, 5),
                  f = factor(rep(c("a", "b"), 20)))
  fit <- frm(y ~ x, data = d, family = gaussian())
  fitf <- frm(y ~ x + f, data = d, family = gaussian())
  multi <- structure(list(), class = "frmtmb_multiple")
  short <- rnorm(3)
  dl <- d
  dl$lc <- replicate(40, list(1:2), simplify = FALSE)
  cases <- list(
    predict_type = list(quote(predict(fit, type = "bogus")),
                        "`type` must be one of \"link\""),
    fitted_scale = list(quote(fitted(fit, scale = "bogus")), "`scale`"),
    residuals_type = list(quote(residuals(fit, type = "bogus")), "`type`"),
    confint_method = list(quote(confint(fit, method = "bogus")),
                          "`method`"),
    drop1_test = list(quote(drop1(fit, test = "bogus")), "`test`"),
    hypothesis_method = list(quote(hypothesis(fit, "x = 0",
                                              method = "bogus")),
                             "`method`"),
    hypothesis_scope = list(quote(hypothesis(fit, "x = 0",
                                             scope = "bogus")), "`scope`"),
    ce_band = list(quote(conditional_effects(fit, band = "bogus")),
                   "`band`"),
    ce_method = list(quote(conditional_effects(fit, method = "bogus")),
                     "`method`"),
    control_nlev = list(quote(frmtmb_control(check_nlev_1 = "bogus")),
                        "`check_nlev_1`"),
    control_olre = list(quote(frmtmb_control(check_olre = "bogus")),
                        "`check_olre`"),
    default_prior_route = list(quote(default_prior(y ~ x, data = d,
                                                   route = "bogus")),
                               "`route`"),
    validate_prior_route = list(quote(validate_prior(
      set_prior("normal(0, 1)"), y ~ x, data = d, route = "bogus")),
      "`route`"),
    vcov_cluster_type = list(quote(vcov_cluster(fit, ~ g, type = "bogus")),
                             "`type`"),
    periodogram_taper = list(quote(frm_periodogram(rnorm(64),
                                                   taper = "bogus")),
                             "`taper`"),
    periodogram_detrend = list(quote(frm_periodogram(rnorm(64),
                                                     detrend = "bogus")),
                               "`detrend`"),
    anova_multiple_method = list(quote(anova(multi, method = "bogus")),
                                 "`method`"),
    anova_multiple_use = list(quote(anova(multi, use = "bogus")), "`use`"),
    set_prior_numeric = list(quote(set_prior(1)), "one string"),
    set_prior_normal = list(quote(set_prior("normal(0)")),
                            "gives 1 argument; normal[(][)] takes 2"),
    set_prior_student_t = list(quote(set_prior("student_t(3, 0)")),
                               "student_t[(]nu, mu, sigma[)]"),
    set_prior_cauchy = list(quote(set_prior("cauchy(0)")), "takes 2"),
    set_prior_exponential = list(quote(set_prior("exponential(1, 2)")),
                                 "gives 2 arguments; exponential"),
    set_prior_exponential_rate = list(quote(set_prior("exponential(-1)")),
                                      "must be positive, not -1"),
    set_prior_lkj = list(quote(set_prior("lkj()")), "gives 0 arguments"),
    set_prior_logistic = list(quote(set_prior("logistic(0)")), "takes 2"),
    set_prior_gamma = list(quote(set_prior("gamma(1)")), "takes 2"),
    set_prior_inv_gamma = list(quote(set_prior("inv_gamma(1)")), "takes 2"),
    set_prior_beta = list(quote(set_prior("beta(1)")), "takes 2"),
    prior_plus_string = list(quote(set_prior("normal(0, 1)") + "x"),
                             "one side is character \"x\""),
    frm_prior_list = list(quote(frm(y ~ x, data = d, prior = list(1))),
                          "`prior` must be"),
    diagnose_nonfit = list(quote(diagnose(lm(y ~ x, d))),
                           "fitted by frm[(][)], not a lm"),
    frmtmb_family_bad = list(quote(frmtmb_family(
      family = 1, dpars = "mu", links = list(mu = "identity"),
      lpdf = function(y, dpars) y)), "family =[)] must be one string"),
    check_custom_family = list(quote(check_custom_family(1)),
                               "not numeric \"1\""),
    cluster_scores_nonfit = list(quote(cluster_scores(lm(y ~ x, d), ~ g)),
                                 "cluster_scores[(][)] needs"),
    vcov_cluster_nonfit = list(quote(vcov_cluster(lm(y ~ x, d), ~ g)),
                               "vcov_cluster[(][)] needs"),
    register_prior_defaults = list(
      quote(frmtmb:::frmtmb_register_prior_defaults(1)), "needs a function"),
    predict_missing_column = list(quote(predict(fit,
                                                newdata = d[, c("y", "g")])),
                                  "Variable 'x' missing from newdata"),
    predict_new_level = list(quote(predict(fitf, newdata = transform(
      d[1:4, ], f = factor(c("a", "zz", "a", "b"))))),
      "a level of `f` that the fit did not see: 'zz'"),
    frm_unknown_variable = list(quote(frm(y ~ nope, data = d)),
                                "The model uses `nope`"),
    frm_data_not_df = list(quote(frm(y ~ x, data = "notdf")),
                           "`data` must be a data frame"),
    frm_env_length = list(quote(frm(bf(y ~ b0 * short, b0 ~ 1, nl = TRUE),
                                    data = d, dry_run = "frame",
                                    start = list(beta = 1))),
                          "has 3 rows and `data` has 40"),
    frm_list_column = list(quote(frm(y ~ lc, data = dl)),
                           "`lc` is a list")
  )
  n <- 0L
  for (nm in names(cases)) {
    e <- tryCatch(eval(cases[[nm]][[1L]]), error = identity)
    expect_identical(class(e), c("frmtmb_error", "error", "condition"),
                     label = nm)
    expect_match(conditionMessage(e), cases[[nm]][[2L]], label = nm)
    n <- n + 1L
  }
  expect_identical(n, length(cases))
  # a newdata that lacks nothing and has no new level still predicts,
  # which is the case a refusal must not fire on
  expect_no_error(predict(fitf, newdata = d[1:4, ]))
})

test_that("the frame check searches an environment data as model.frame()", {
  # review recheck R1: model.frame() evaluates in an environment `data`
  # and its parents, never in the formula environment
  set.seed(20260917)
  d <- data.frame(y = rnorm(60), x = rnorm(60), xq = rnorm(60))
  pe <- new.env()
  assign("xq", d$xq, envir = pe)
  de <- new.env(parent = pe)
  assign("y", d$y, envir = de)
  assign("x", d$x, envir = de)
  fit <- frm(y ~ x + xq, data = de)
  expect_identical(as.numeric(logLik(fit)),
                   as.numeric(logLik(frm(y ~ x + xq, data = d))))
})

test_that("an unsupported term and `.` are refused as what they are", {
  # review recheck R2: the frame check reported both as an unknown
  # variable
  set.seed(20260917)
  d <- data.frame(y = rnorm(60), x = rnorm(60))
  Wm <- diag(60)
  expect_error(frm(y ~ x + sar(Wm), data = d, data2 = list(Wm = Wm)),
               "sar[(][)] is a brms autocorrelation term that frmtmb",
               class = "frmtmb_error")
  expect_error(frm(y ~ x + fcor(Vm), data = d, data2 = list(Vm = Wm)),
               "fcor[(][)] is a brms autocorrelation term",
               class = "frmtmb_error")
  expect_error(frm(y ~ ., data = d), "A formula with `.` is not supported",
               class = "frmtmb_error")
})
