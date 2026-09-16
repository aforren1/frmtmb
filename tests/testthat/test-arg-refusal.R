# Plan item 2.5e. Two things are asserted here and they are different
# claims.
#
# The INSTANCE: `fitted()` took only `...`, so `fitted(fit, re_formula =
# NA)` returned the conditional fit and said nothing. On the fixture
# below that is 1.291651 at the first row, bit-identical to
# `fitted(fit)`, where the population answer the caller asked for is
# 1.954262 and the largest gap over the 120 rows is 3.065438. A
# misspelled argument was swallowed the same way. Measured against
# frmtmb 0.57.0 by dev/argspell-defect.R, seeds 2505 and 25051.
#
# The CLASS: a method whose `...` is never read accepts every name at
# all, so it will take every argument brms adds in future in the same
# silence. The last block is the structural guard against that, and it
# is the reason this file exists rather than three assertions bolted
# onto test-methods.R.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

# One random-intercept gaussian fit, built once. The group effects are
# large relative to sigma on purpose: the defect is invisible when the
# conditional and population answers are close.
ar_case <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(2505)
      dd <- data.frame(x = stats::rnorm(120),
                       g = factor(rep(1:8, each = 15)), y = 0)
      dd$y <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                           newparams = list(Intercept = 1, x = 0.5,
                                            sigma = 0.4,
                                            sd_g__Intercept = 1.2),
                           nsim = 1, seed = 25051)[[1]]
      cache <<- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
    }
    cache
  }
})

# ---- the instance ----------------------------------------------------

test_that("fitted() honors re_formula instead of swallowing it", {
  fit <- ar_case()
  cond <- fitted(fit)
  pop <- fitted(fit, re_formula = NA)

  # the two answers are different numbers, which is the whole defect
  expect_false(isTRUE(all.equal(unname(cond), unname(pop))))

  # and the population one is predict()'s, exactly
  expect_equal(unname(pop),
               unname(predict(fit, re_formula = NA, type = "response")))
  expect_equal(unname(cond),
               unname(predict(fit, type = "response")))

  # the gap at every row is that row's fitted group mode, so the size of
  # the difference is accounted for rather than merely observed
  b <- ranef(fit)[["g"]][, 1L]
  expect_equal(unname(cond - pop),
               unname(b[as.integer(fit[["frame"]][["data_frame"]][["g"]])]),
               tolerance = 1e-8)
})

test_that("fitted() reaches newdata, scale, dpar and resp", {
  fit <- ar_case()
  nd <- fit[["frame"]][["data_frame"]][1:5, ]

  expect_equal(unname(fitted(fit, newdata = nd)),
               unname(fitted(fit))[1:5])
  expect_equal(unname(fitted(fit, scale = "linear")),
               unname(predict(fit, type = "link")))
  expect_equal(unname(fitted(fit, dpar = "sigma")),
               unname(predict(fit, dpar = "sigma", type = "response")))
  expect_equal(unname(fitted(fit, resp = "y")), unname(fitted(fit)))
  # match.arg() names the permitted values rather than the argument
  expect_error(fitted(fit, scale = "latent"), "should be one of")
})

# ---- the refusal -----------------------------------------------------

test_that("a misspelled argument is refused and named", {
  fit <- ar_case()

  expect_error(fitted(fit, re_frmula = NA), "re_frmula")
  expect_error(fitted(fit, re_frmula = NA), "Did you mean")
  expect_error(predict(fit, re_frmula = NA), "re_frmula")
  expect_error(coef(fit, robst = TRUE), "robst")
  expect_error(ranef(fit, condVarr = TRUE), "condVarr")
  expect_error(logLik(fit, REML = TRUE), "REML")

  # an argument with no name at all has nowhere to go either. nobs()
  # is the case to use: fitted() has positional formals, so a bare 7
  # matches `newdata` and is refused for being the wrong type instead.
  expect_error(nobs(fit, 7), "no name")
})

test_that("an argument brms has and a point fit cannot honor says why", {
  fit <- ar_case()
  for (a in c("ndraws", "draw_ids", "summary", "robust", "probs", "sort",
              "nlpar")) {
    cl <- list(quote(stats::fitted), quote(fit))
    cl[[a]] <- 1
    msg <- tryCatch(eval(as.call(cl)), error = conditionMessage)
    expect_true(is.character(msg), info = a)
    expect_match(msg, a, fixed = TRUE, info = a)
  }
  # the draws-specific ones say where the argument does work
  expect_match(tryCatch(fitted(fit, ndraws = 10), error = conditionMessage),
               "frmtmb.sample")
})

test_that("a base-R path that forwards its own argument still works", {
  fit <- ar_case()
  # stats::step(), add1.default(), drop1.default() and sigma.default()
  # all call nobs(object, use.fallback = TRUE). `use.fallback` is
  # stats::nobs.default's own formal, so refusing it broke step(), which
  # worked before the dots were refused. Nothing else in the suite
  # reaches step().
  expect_identical(nobs(fit, use.fallback = TRUE), nobs(fit))
  expect_no_error(suppressWarnings(stats::step(fit, trace = 0)))
  # and the refusal is still on: a name stats does NOT pass is refused
  expect_error(nobs(fit, use.fallbck = TRUE), "use.fallbck")
})

test_that("a brms argument is refused with its reason, not as unknown", {
  fit <- ar_case()
  # the principle fitted_no_draws states, applied to the rest of the
  # surface: each of these is a REAL argument of the brmsfit method of
  # the same generic, so "no such argument" would send the caller after
  # a typo that is not there. All of them were accepted in silence
  # before the dots were refused.
  cases <- list(
    list(quote(print(fit, digits = 3)), "digits", "not implemented"),
    list(quote(print(summary(fit), digits = 3)), "digits",
         "not implemented"),
    list(quote(summary(fit, prob = 0.9)), "prob", "confint"),
    list(quote(summary(fit, priors = TRUE)), "priors", "prior_summary"),
    list(quote(fixef(fit, pars = "x")), "pars", "regular expression"),
    list(quote(coef(fit, robust = TRUE)), "robust", "draws"),
    list(quote(ranef(fit, groups = "g")), "groups", "named list"),
    list(quote(VarCorr(fit, probs = 0.5)), "probs", "quantiles"),
    list(quote(nobs(fit, resp = "y")), "resp", "same nobs"),
    list(quote(family(fit, resp = "y")), "resp", "NAMED LIST"),
    list(quote(vcov(fit, correlation = TRUE)), "correlation", "cov2cor")
  )
  for (cs in cases) {
    msg <- tryCatch(eval(cs[[1L]]), error = conditionMessage)
    expect_true(is.character(msg), info = deparse(cs[[1L]]))
    expect_match(msg, cs[[2L]], fixed = TRUE, info = deparse(cs[[1L]]))
    expect_match(msg, cs[[3L]], fixed = TRUE, info = deparse(cs[[1L]]))
  }
  # and the method name in the message is one the caller could type:
  # print.summary.frmtmb_fit used to report `print.summary()`
  msg <- tryCatch(print(summary(fit), digits = 3), error = conditionMessage)
  expect_match(msg, "print()", fixed = TRUE)
  expect_false(grepl("print.summary(", msg, fixed = TRUE))
})

test_that("predict() and simulate() take brms's re_formula, not lme4's", {
  for (nm in c("predict.frmtmb_fit", "simulate.frmtmb_fit")) {
    fo <- names(formals(getFromNamespace(nm, "frmtmb")))
    expect_true("re_formula" %in% fo, info = nm)
    expect_false("re.form" %in% fo, info = nm)
  }
  fo <- names(formals(frm_bootstrap))
  expect_true("re_formula" %in% fo)
  expect_false("re.form" %in% fo)

  fit <- ar_case()
  # the dropped lme4 spellings are refused, not ignored
  expect_error(predict(fit, re.form = NA), "re.form", fixed = TRUE)
  expect_error(predict(fit, allow.new.levels = TRUE),
               "allow.new.levels", fixed = TRUE)
  expect_error(simulate(fit, nsim = 1, re.form = NA), "re.form",
               fixed = TRUE)
})

# ---- the class -------------------------------------------------------

# A method that takes `...` and never mentions it again swallows
# everything. `frm_check_dots(...)` is one such mention, so the test is
# "the body refers to its dots somehow", which a pass-through method
# satisfies by forwarding them.
#
# The test reads the PARSE TREE, not the deparsed text. A text search
# for the three characters matches them inside a STRING, and one method
# was scored as a dots user on the strength of
# `cat("... ", length(v) - n, " more\n")`: `print.frmtmb_par_template()`
# swallowed everything while the detector called it clean, so the
# generated "0 methods still swallowing" was false.
ar_dots_names <- c("...", "...length", "...names", "...elt")

# `..1` through `..N`, by pattern: enumerating three of them scored
# `..4` and `..11` as swallowers. Nothing in either package uses one
# today, which is why it would not have been caught by a run.
ar_is_dots_name <- function(x) {
  x %in% ar_dots_names | grepl("^[.][.][0-9]+$", x)
}

# One shape is exempt and is recognised rather than listed: a method
# whose whole body is an unconditional refusal cannot swallow anything,
# because EVERY call to it errors. Guarding those made the message
# worse, and on one it made it wrong:
# `conditional_effects(frm_multiple_result, "x")` reported "1 argument
# with no name" where "no pooled version" is the answer.
ar_refusers <- c("stop", "fit_no_draws", "multiple_no_draws")

# Under covr every statement arrives as
# `if (TRUE) { covr:::count(key); <statement> }`, so a one-call body no
# longer LOOKS like one call and all 20 refusers read as swallowers in
# the coverage job only. This removes exactly that wrapper, at the top
# level, which is the only level the shape checks read; anything else
# is returned untouched, so a real branch still fails them.
ar_uncovr <- function(e) {
  if (is.call(e) && identical(e[[1L]], as.name("if")) &&
        length(e) == 3L && isTRUE(e[[2L]])) {
    i <- e[[3L]]
    if (is.call(i) && identical(i[[1L]], as.name("{")) &&
          length(i) == 3L && is.call(i[[2L]]) &&
          identical(i[[2L]][[1L]], quote(covr:::count))) {
      return(i[[3L]])
    }
  }
  e
}

ar_body <- function(fn) {
  b <- ar_uncovr(body(fn))
  if (is.call(b) && identical(b[[1L]], as.name("{")) && length(b) > 1L) {
    for (k in 2:length(b)) b[[k]] <- ar_uncovr(b[[k]])
  }
  b
}

ar_refuses_always <- function(fn) {
  b <- ar_body(fn)
  head_is_refusal <- function(e) {
    is.call(e) && as.character(e[[1L]])[1L] %in% ar_refusers
  }
  if (!is.call(b) || !identical(as.character(b[[1L]]), "{")) {
    return(head_is_refusal(b))
  }
  length(b) == 2L && head_is_refusal(b[[2L]])
}

ar_swallows <- function(fn) {
  fo <- names(formals(fn))
  if (!"..." %in% fo) return(FALSE)
  if (any(ar_is_dots_name(all.names(body(fn))))) return(FALSE)
  !ar_refuses_always(fn)
}

# Generics another package owns and calls itself. emmeans, insight and
# marginaleffects pass arguments of their own choosing into these
# methods, so accepting unknown names IS the contract and refusing them
# would break the caller, not the caller's typo.
#
# An exemption is a hole in the guard, so each name here has evidence
# behind it. `dev/argspell-exempt2.R` counts, per generic, the WRITTEN
# call sites in those three packages that pass a name our method has no
# formal for. The list was thirteen and is nine:
# `recover_data.frmtmb_fit` left because it READS its dots, so exempting
# it did nothing, and three left because they had no need by either
# measurement AND survive the runtime harness.
#
# `get_coef` and `set_coef` are here on RUNTIME evidence and not on a
# count. Both static passes scored them 0, and guarding them broke four
# marginaleffects entry points, because the call is assembled rather
# than written: `get_coef() has no argument 'variables' (and 2 more:
# numderiv, internal_call)`. `dev/argspell-interop.R` is the harness
# and `dev/argspell-interop-log.txt` the run. A count of zero written
# call sites is not evidence that none exists.
ar_dots_exempt <- c(
  "emm_basis.frmtmb_fit", "get_predict.frmtmb_fit",
  "get_vcov.frmtmb_fit", "get_varcov.frmtmb_fit",
  "get_parameters.frmtmb_fit", "find_formula.frmtmb_fit",
  "find_random.frmtmb_fit",
  # not found by either static count; found by RUNNING marginaleffects
  "get_coef.frmtmb_fit", "set_coef.frmtmb_fit"
)

test_that("the swallow detector fires on a swallower and not on a user", {
  # constructed in its ABSENT case first: a guard that cannot report a
  # positive has never been shown to work
  expect_true(ar_swallows(function(object, ...) object))
  expect_false(ar_swallows(function(object, ...) list(...)))
  expect_false(ar_swallows(function(object) object))

  # the case that got past the text detector: three dots inside a
  # STRING are not a use of the dots. This is the absent case for the
  # parse-tree rewrite, and it is the one that was missing.
  expect_true(ar_swallows(function(object, ...) {
    message("a ... b")
    object
  }))
  expect_true(ar_swallows(function(object, ...) {
    cat("... ", 1, " more\n", sep = "")
    object
  }))
  # and the two real uses that a crude regex would have to get right
  # for the rewrite to be safe both ways
  expect_false(ar_swallows(function(x, i, ...) unclass(x)[[i, ...]]))
  expect_false(ar_swallows(function(object, ...) {
    match.call(expand.dots = FALSE)$...
  }))
  expect_false(ar_swallows(function(object, ...) ...length()))

  # the unconditional-refusal shape, in both its spellings, and the
  # near miss that must NOT be exempted: a refusal reached only on a
  # branch leaves the other branches swallowing
  expect_false(ar_swallows(function(object, ...) stop("no")))
  expect_false(ar_swallows(function(object, ...) {
    stop("no")
  }))
  expect_true(ar_swallows(function(object, ...) {
    if (is.null(object)) stop("no")
    object
  }))
})

test_that("R's own S3 contract is honored, and the scan finds the next hop", {
  # The rule "refuse what the METHOD does not have" collides with R's
  # "a method tolerates what its GENERIC carries". Three rounds of
  # one-off .allow patches each uncovered the next hop: nobs(use.fallback)
  # from step(), drop1(scale =, trace =) from the same call,
  # formula(env =) from as.formula(), terms(data =) from
  # model.matrix.default(). This walks R's own packages for the pattern
  # instead of waiting for the next report.
  scan <- testthat::test_path("..", "..", "dev",
                              "argspell-contract-scan.R")
  skip_if(!file.exists(scan), "dev/ not reachable from the test dir")
  sys.source(scan, envir = environment())

  scan_pkgs <- c("stats", "base", "utils", "graphics", "grDevices",
                 "methods")
  scan_pkgs <- Filter(function(p) requireNamespace(p, quietly = TRUE),
                      scan_pkgs)
  pkgs <- Filter(function(p) nzchar(system.file(package = p)),
                 c("frmtmb", "frmtmb.sample"))
  acc <- frm_contract_scan(pkgs, scan_pkgs,
                           contract = frmtmb:::s3_contract_args)

  # non-vacuous in two ways: the scan really reached the generics, and
  # with the table EMPTIED it really reports the hops the table exists
  # for. A scan that found nothing either way would pass silently.
  expect_gt(acc$n_generics, 40L)
  bare <- frm_contract_scan(pkgs, scan_pkgs, contract = NULL)
  expect_true(any(bare$hits$generic == "formula" &
                    bare$hits$arg == "env"))
  expect_true(any(bare$hits$generic == "nobs" &
                    bare$hits$arg == "use.fallback"))

  miss <- if (nrow(acc$hits)) {
    paste(acc$hits$generic, acc$hits$arg, acc$hits$where)
  } else {
    character()
  }
  expect_identical(miss, character())
})

test_that("the contract tolerates R's names and refuses everything else", {
  fit <- ar_case()
  # what the table buys, end to end
  expect_no_error(stats::as.formula(fit))
  expect_no_error(formula(fit, env = baseenv()))
  expect_no_error(terms(fit, data = fit[["frame"]][["data_frame"]]))
  expect_no_error(model.frame(fit,
                              data = fit[["frame"]][["data_frame"]]))
  expect_no_error(confint(fit, trace = FALSE, parm = "x"))
  expect_no_error(residuals(fit, na.rm = TRUE))
  expect_no_error(coef(fit, complete = TRUE))
  expect_no_error(nobs(fit, use.fallback = TRUE))
  expect_identical(nobs(fit, use.fallback = TRUE), nobs(fit))

  # and what it does NOT buy: a name R never passes is still refused,
  # on the very same methods
  expect_error(formula(fit, envv = baseenv()), "envv")
  expect_error(nobs(fit, nosucharg = 1), "nosucharg")
  expect_error(terms(fit, re.form = NA), "re.form", fixed = TRUE)
  expect_error(coef(fit, complet = TRUE), "complet")
  # an `.unsupported` reason still wins over the contract
  expect_error(print(fit, digits = 3), "not implemented")
})

test_that("frm_check_dots() is only ever called from a dots-taker", {
  # It reads the CALLING FRAME for the name and the formals to report.
  # One frame down, from a helper the method delegates to, it names the
  # helper and prints an empty formal list. Today no call site does
  # that; this fails if one ever appears, instead of degrading quietly.
  #
  # Both roots hang off the package root, and that root must PROVE it is
  # frmtmb's source. Under R CMD check `../..` is the .Rcheck directory,
  # and the old `../../../..` for frmtmb.sample reached a real
  # extensions/ only when the check directory happened to sit two levels
  # inside the repository, as it does on CI. That run scanned the
  # extension alone, found 28 call sites and no definition, and failed;
  # a source run scanned core alone, because the same path pointed above
  # the repository.
  pkg_root <- testthat::test_path("..", "..")
  desc <- file.path(pkg_root, "DESCRIPTION")
  is_src <- file.exists(desc) && dir.exists(file.path(pkg_root, "R")) &&
    identical(unname(read.dcf(desc, "Package")[1L, 1L]), "frmtmb")
  skip_if(!is_src, "frmtmb's own sources are not beside the tests")
  roots <- c(file.path(pkg_root, "R"),
             file.path(pkg_root, "extensions", "frmtmb.sample", "R"))
  files <- unlist(lapply(roots, function(d) {
    if (dir.exists(d)) list.files(d, pattern = "[.]R$", full.names = TRUE)
  }))

  # Counted two ways rather than walked: `all.names()` finds EVERY
  # mention of the helper anywhere in a file, and the structured pass
  # finds only those that are a top-level statement of a top-level
  # function that takes `...`. Equal counts mean there is no third
  # kind. A recursive walker is the wrong tool here: `x[i, ]` holds an
  # empty symbol and binding it raises "argument is missing".
  anywhere <- 0L
  positioned <- 0L
  defined <- 0L
  for (fp in files) {
    ex <- parse(fp, keep.source = FALSE)
    for (i in seq_along(ex)) {
      e <- ex[[i]]
      is_def <- is.call(e) && identical(as.character(e[[1L]]), "<-") &&
        is.name(e[[2L]]) &&
        identical(as.character(e[[2L]]), "frm_check_dots")
      if (is_def) {
        # the helper's own definition mentions its name once, as the
        # assignment target, and is not a call site
        defined <- defined + 1L
        next
      }
      anywhere <- anywhere + sum(all.names(e) == "frm_check_dots")
      if (!is.call(e) || !identical(as.character(e[[1L]]), "<-")) next
      if (!is.name(e[[2L]])) next
      rhs <- e[[3L]]
      if (!is.call(rhs) ||
            !identical(as.character(rhs[[1L]]), "function")) next
      if (!"..." %in% names(rhs[[2L]])) next
      b <- rhs[[3L]]
      stmts <- if (is.call(b) && identical(as.character(b[[1L]]), "{")) {
        as.list(b)[-1L]
      } else {
        list(b)
      }
      for (st in stmts) {
        if (is.call(st) && is.name(st[[1L]]) &&
              identical(as.character(st[[1L]]), "frm_check_dots")) {
          positioned <- positioned + 1L
        }
      }
    }
  }
  # non-vacuous: the sources really do contain call sites to check
  expect_gt(anywhere, 50L)
  expect_identical(positioned, anywhere)
  # and the helper was defined exactly once, so subtracting that one
  # mention cannot be hiding a call site
  expect_identical(defined, 1L)
})

test_that("the two refusal helpers really are unconditional", {
  # 38 methods are exempt from the swallow guard because their whole
  # body is one call to one of these. If either grows a branch, those
  # 38 become exempt SWALLOWERS with nothing failing.
  for (nm in c("fit_no_draws", "multiple_no_draws")) {
    fn <- getFromNamespace(nm, "frmtmb")
    b <- ar_body(fn)
    expect_true(is.call(b) && identical(as.character(b[[1L]]), "{"),
                label = paste(nm, "has a braced body"))
    expect_identical(length(b), 2L, label = paste(nm, "has one statement"))
    expect_identical(as.character(b[[2L]][[1L]]), "stop",
                     label = paste(nm, "statement is stop()"))
  }
  # the absent case: a helper with a branch must NOT satisfy this
  branchy <- function(fn) {
    b <- ar_body(fn)
    is.call(b) && identical(as.character(b[[1L]]), "{") &&
      length(b) == 2L && identical(as.character(b[[2L]][[1L]]), "stop")
  }
  expect_false(branchy(function(x) {
    if (x) stop("a")
    invisible(NULL)
  }))
  # and covr's wrapper removal does not launder a branch: the same
  # body in the shape covr writes still fails, as does a lone branch
  expect_true(branchy(function(x) {
    if (TRUE) {
      covr:::count("k")
      stop("a")
    }
  }))
  expect_false(branchy(function(x) {
    if (TRUE) {
      covr:::count("k")
      if (x) stop("a")
    }
  }))
  expect_false(ar_refuses_always(function(x, ...) {
    if (TRUE) {
      covr:::count("k")
      if (x) fit_no_draws("f")
    }
  }))
})

test_that("no registered S3 method has a `...` it never touches", {
  ns <- asNamespace("frmtmb")
  s3 <- parseNamespaceFile("frmtmb", dirname(system.file(package = "frmtmb")))
  reg <- s3$S3methods
  nm <- ifelse(is.na(reg[, 3]), paste(reg[, 1], reg[, 2], sep = "."),
               reg[, 3])
  nm <- setdiff(unique(nm), ar_dots_exempt)
  bad <- Filter(function(f) {
    obj <- tryCatch(get(f, envir = ns), error = function(e) NULL)
    !is.null(obj) && ar_swallows(obj)
  }, nm)
  expect_identical(bad, character(0))
})
