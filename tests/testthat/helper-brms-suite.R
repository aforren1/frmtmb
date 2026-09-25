# Harness for the brms suite port (item 2.6b): brms 2.23.0's own bin-1
# assertions, run against frmtmb by the generated files
# test-brms-suite-*.R. dev/brmsport-gen.R writes those files and
# dev/brmsport-ledger.R turns a recorded run into dev/brmsport-ledger.tsv.
#
# Every brms assertion is wrapped in brms_port() with the verdict the
# ledger gives it. A "pass" must hold. Anything else (a defect, a
# deliberate divergence, pending 2.6d or 2.6e, cannot transfer) must
# NOT hold, so the tier stays green while it cannot be made to pass by
# loosening the port, and it goes red the day the verdict goes stale.

# The tier is gated with the other brms fit tiers. It compiles nothing
# on its own, but its fixtures fit models and the sample side draws.
skip_unless_brms_suite <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("brms")
  if (!identical(Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "true")) {
    testthat::skip("set FRMTMB_BRMS_FIT_TESTS=true to run the brms suite")
  }
}

# An error of this kind means an assertion never reached frmtmb code, so
# an expect_error() that catches it is not a pass. The famlink lane's
# ledger found ten such passes on base.
brms_vacuous_pattern <- paste0(
  "could not find function|object '[^']*' not found|",
  "of mode 'function' was not found|",
  "there is no package called|is not an exported object from|",
  "is missing, with no default"
)

# frmtmb refusing an ARGUMENT's name or saying it cannot honor one says
# nothing about the value brms refuses, whatever brms's pattern matches
# in the message (a regex often names the argument too). An assertion
# whose own pattern is about such a refusal is exempt.
brms_argument_refusal <- paste0(
  "has no argument|unused argument|cannot honor|",
  "Cannot interpret [^ ]+ argument"
)

brms_about_argument_refusal <- function(call) {
  args <- as.list(call)[-1]
  chr <- Filter(is.character, args[-1])
  length(chr) > 0L && any(vapply(chr, function(x) {
    grepl(brms_argument_refusal, x)
  }, TRUE))
}

# A `$` read that does not reach the named element of a plain list: a
# PARTIAL match (fit$data is fit$data2), or an absent name read as an
# intermediate step (the chain goes on reading NULL). A final absent
# name is not flagged, since brms may mean it (bprior$x is NULL). Only
# plain lists are checked: an object with its own `$` method answers
# names() for something else. The object prefix is evaluated again,
# so this runs only on the few assertion forms that need it.
brms_hollow_read <- function(e, env, intermediate = FALSE) {
  if (!is.call(e)) return(FALSE)
  head <- as.character(e[[1]])[1]
  if (identical(head, "$") && length(e) == 3L) {
    nm <- as.character(e[[3]])
    obj <- tryCatch(eval(e[[2]], env), error = function(err) NULL)
    plain <- is.list(obj) && !is.data.frame(obj) &&
      !any(vapply(class(obj), function(k) {
        !is.null(utils::getS3method("$", k, optional = TRUE))
      }, TRUE))
    nms <- names(obj)
    if (is.null(nms)) nms <- character()
    if (plain && !nm %in% nms) {
      partial <- any(startsWith(nms, nm))
      if (partial || intermediate) return(TRUE)
    }
    return(brms_hollow_read(e[[2]], env, intermediate = TRUE))
  }
  args <- as.list(e)[-1]
  for (a in args) {
    if (is.name(a) && !nzchar(as.character(a))) next
    if (brms_hollow_read(a, env, intermediate = FALSE)) return(TRUE)
  }
  FALSE
}

# expect_null(), expect_true(is.null()), expect_length() and a NULL
# comparison held on a hollow `$` read. Two NULLs from a genuine source
# (names() of unnamed vectors) are not flagged.
brms_hollow_null <- function(call, val, env) {
  if (!is.call(call)) return(FALSE)
  head <- as.character(call[[1]])[1]
  target <- switch(
    head,
    expect_null = call[[2]],
    expect_length = call[[2]],
    expect_true = if (is.call(call[[2]]) &&
                      identical(call[[2]][[1]], as.name("is.null"))) {
      call[[2]][[2]]
    },
    expect_equal = , expect_identical = , expect_equivalent =
      if (length(call) >= 3L && is.null(val) &&
          !identical(call[[3]], quote(NULL))) {
        call("list", call[[2]], call[[3]])
      },
    NULL
  )
  !is.null(target) && brms_hollow_read(target, env)
}

brms_port_state <- new.env(parent = emptyenv())

brms_port_clean <- function(x) {
  x <- paste(x, collapse = " ")
  x <- gsub("[\r\n\t]+", " ", x)
  substr(x, 1L, 400L)
}

# Evaluates one assertion and reports whether it HELD. Successes are
# muffled so the reporter counts the harness's own expectation once,
# and warnings and messages that brms's assertion did not itself
# capture are muffled as testthat would merely report them.
brms_port_run <- function(call, env) {
  held <- TRUE
  msg <- ""
  val <- tryCatch(
    withCallingHandlers(
      eval(call, env),
      expectation_success = function(e) {
        invokeRestart("muffle_expectation")
      },
      warning = function(w) invokeRestart("muffleWarning"),
      message = function(m) invokeRestart("muffleMessage")
    ),
    expectation_failure = function(e) {
      held <<- FALSE
      msg <<- conditionMessage(e)
      NULL
    },
    skip = function(e) {
      held <<- FALSE
      msg <<- paste("SKIP:", conditionMessage(e))
      NULL
    },
    error = function(e) {
      held <<- FALSE
      msg <<- paste("ERROR:", conditionMessage(e))
      NULL
    }
  )
  # testthat 3e returns the condition an expect_error(), expect_warning()
  # or expect_message() caught, which is how a vacuous pass is seen
  # without evaluating the call a second time
  raw_held <- held
  caught <- if (inherits(val, "condition")) conditionMessage(val) else ""
  vacuous <- grepl(brms_vacuous_pattern, if (held) caught else msg)
  why <- caught
  if (held && grepl(brms_argument_refusal, caught) &&
      !brms_about_argument_refusal(call)) {
    vacuous <- TRUE
  }
  if (held && brms_hollow_null(call, val, env)) {
    vacuous <- TRUE
    why <- "a NULL read through a partial or intermediate $ match"
  }
  if (held && vacuous) {
    held <- FALSE
    msg <- paste("VACUOUS:", why)
  }
  list(held = held, vacuous = vacuous, msg = brms_port_clean(msg),
       caught = brms_port_clean(caught), raw_held = raw_held)
}

brms_port_record <- function(kind, id, verdict, res) {
  path <- Sys.getenv("FRMTMB_BRMSPORT_RECORD")
  pkg <- Sys.getenv("FRMTMB_BRMSPORT_PKG", "frmtmb")
  # raw_held: whether testthat alone would have counted it a pass
  raw <- if (is.null(res$raw_held)) res$held else res$raw_held
  row <- c(kind, pkg, id, verdict, res$held, res$vacuous, res$msg,
           res$caught, raw)
  cat(paste(row, collapse = "\t"), "\n", sep = "", file = path,
      append = TRUE)
}

# The names a call assigns with `<-` or `=`, at any depth, so that a
# failed setup line or a failed assertion that assigns can mark them.
# A complex target (x$a, x[[i]], names(x)) assigns its root object only
# in part, so `complex = TRUE` asks for those roots separately.
brms_assign_root <- function(target) {
  while (is.call(target)) target <- target[[2]]
  if (is.name(target)) as.character(target) else character()
}

brms_assigned <- function(e, complex = FALSE) {
  if (!is.call(e)) return(character())
  out <- character()
  if (!complex && identical(as.character(e[[1]])[1], "assign") &&
      length(e) >= 3L && is.character(e[[2]])) {
    out <- e[[2]]
  }
  if (as.character(e[[1]])[1] %in% c("<-", "=") && length(e) == 3L) {
    if (is.name(e[[2]]) && !complex) out <- as.character(e[[2]])
    if (is.call(e[[2]]) && complex) out <- brms_assign_root(e[[2]])
  }
  args <- as.list(e)[-1]
  for (i in seq_along(args)) {
    # the empty argument of x[i, ] cannot be bound to a name
    if (is.name(args[[i]]) && !nzchar(as.character(args[[i]]))) next
    out <- c(out, brms_assigned(args[[i]], complex))
  }
  unique(out)
}

# An object whose last assignment FAILED still holds whatever an earlier
# line put there, so an assertion reading it would test the earlier
# object. brms's block would have stopped at the failure instead. Such a
# name is stale until a later line assigns it successfully.
brms_mark <- function(call, ok, msg) {
  env <- parent.frame(2L)
  st <- brms_port_state$stale
  if (is.null(st) || !identical(brms_port_state$env, env)) {
    st <- list()
    brms_port_state$env <- env
  }
  for (nm in brms_assigned(call)) {
    st[[nm]] <- if (ok) NULL else msg
  }
  # a failed x$a <- leaves x half assigned; a successful one does not
  # repair an x that was already stale
  if (!ok) {
    for (nm in brms_assigned(call, complex = TRUE)) st[[nm]] <- msg
  }
  brms_port_state$stale <- st
}

brms_stale_reads <- function(call, env) {
  st <- if (identical(brms_port_state$env, env)) brms_port_state$stale
  hit <- intersect(all.vars(call), names(st))
  # a name the call itself assigns first is not read stale
  setdiff(hit, brms_assigned(call))
}

brms_port <- function(id, verdict, reason, code) {
  call <- substitute(code)
  env <- parent.frame()
  stale <- brms_stale_reads(call, env)
  res <- brms_port_run(call, env)
  if (length(stale)) {
    res$held <- FALSE
    res$vacuous <- TRUE
    res$msg <- brms_port_clean(sprintf(
      "STALE: reads '%s', whose last assignment failed: %s",
      stale[1], brms_port_state$stale[[stale[1]]]))
  }
  brms_mark(call, res$held, res$msg)
  if (nzchar(Sys.getenv("FRMTMB_BRMSPORT_RECORD"))) {
    brms_port_record("assert", id, verdict, res)
    testthat::succeed()
    return(invisible(res))
  }
  if (verdict %in% c("pass", "hollow")) {
    # a hollow row held on something that is not frmtmb's answer; it is
    # asserted to keep holding so that a change to it is seen
    if (res$held) {
      testthat::succeed()
    } else {
      testthat::fail(sprintf(
        "brms %s is recorded as %s and no longer holds: %s", id, verdict,
        res$msg))
    }
  } else if (res$held) {
    testthat::fail(sprintf(paste(
      "brms %s now HOLDS but is recorded as '%s' (%s). Change its row in",
      "dev/brmsport-verdicts.tsv and regenerate."), id, verdict, reason))
  } else {
    testthat::succeed()
  }
  invisible(res)
}

# brms's refusal assertion rewritten to frmtmb's own message for the same
# call: the regexp becomes `pattern` and `fixed` is dropped. Nothing else
# in the call changes, so the refusal must fire on brms's own case.
brms_own_call <- function(call, pattern) {
  args <- as.list(call)[-1]
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  args <- args[nms != "fixed"]
  nms <- nms[nms != "fixed"]
  rx <- which(nms == "regexp")
  if (!length(rx)) rx <- which(nms == "")[2]
  if (is.na(rx[1])) {
    args <- c(args[1], list(pattern), args[-1])
  } else {
    args[[rx[1]]] <- pattern
  }
  as.call(c(call[[1]], args))
}

# The user's rule of 2026-09-17: a refusal frmtmb makes for the same case,
# in its own words, is a PASS. brms's assertion as written must still NOT
# hold (else the row is a plain pass and this verdict is stale), and the
# same call with frmtmb's own pattern must hold, non-vacuously, so a
# different error (the object-not-found of brm:81) cannot satisfy it.
# Messages a specific own-words pattern must NOT match: the empty string
# and refusals of unrelated cases. A pattern that matches one of these
# (".", "error", "not") would accept any refusal at all.
brms_own_controls <- c(
  "",
  "Error",
  "invalid argument",
  "is not supported",
  "must be",
  paste("x is not a supported family. Supported families are: gaussian,",
        "poisson, binomial"),
  "The model formula could not be evaluated: object not found"
)

brms_port_own <- function(id, pattern, reason, code) {
  call <- substitute(code)
  env <- parent.frame()
  stopifnot(is.call(call), as.character(call[[1]]) %in%
              c("expect_error", "expect_warning", "expect_message"))
  stale <- brms_stale_reads(call, env)
  as_brms <- brms_port_run(call, env)
  own <- brms_port_run(brms_own_call(call, pattern), env)
  res <- own
  res$raw_held <- as_brms$raw_held
  generic <- vapply(brms_own_controls, function(x) grepl(pattern, x),
                    TRUE)
  if (any(generic)) {
    res$held <- FALSE
    res$vacuous <- TRUE
    res$msg <- brms_port_clean(sprintf(
      "OWN-WORDS PATTERN NOT SPECIFIC: '%s' matches the control '%s'",
      pattern, brms_own_controls[which(generic)[1]]))
  } else if (length(stale)) {
    res$held <- FALSE
    res$vacuous <- TRUE
    res$msg <- brms_port_clean(sprintf(
      "STALE: reads '%s', whose last assignment failed: %s",
      stale[1], brms_port_state$stale[[stale[1]]]))
  } else if (as_brms$held) {
    res$held <- FALSE
    res$msg <- "STALE OWN-WORDS: brms's assertion holds as written"
  } else if (!own$held) {
    res$msg <- brms_port_clean(paste("OWN-WORDS FAILED:", own$msg))
  }
  brms_mark(call, res$held, res$msg)
  if (nzchar(Sys.getenv("FRMTMB_BRMSPORT_RECORD"))) {
    brms_port_record("assert", id, "own", res)
    testthat::succeed()
  } else if (res$held) {
    testthat::succeed()
  } else {
    testthat::fail(sprintf(
      "brms %s is recorded as a pass in frmtmb's own words (%s): %s",
      id, pattern, res$msg))
  }
  invisible(res)
}

# An assertion whose text reaches brms's own namespace (brms:::inv_link,
# brms:::family_names). Evaluating it would test brms, and a pass would
# be brms's, so it is recorded and never run.
brms_port_not_run <- function(id, verdict, reason, why) {
  res <- list(held = FALSE, vacuous = FALSE,
              msg = paste("NOT RUN:", why), caught = "")
  if (nzchar(Sys.getenv("FRMTMB_BRMSPORT_RECORD"))) {
    brms_port_record("assert", id, verdict, res)
  } else if (identical(verdict, "pass")) {
    testthat::fail(sprintf("brms %s cannot be a pass: %s", id, why))
  }
  testthat::succeed()
  invisible(res)
}

# A setup line of brms's block. One that fails in frmtmb must not take
# the block's later, independent assertions with it; the assertions that
# do depend on it then fail on the missing object, which the ledger
# reports with this line's error beside it.
brms_setup <- function(id, code) {
  call <- substitute(code)
  res <- list(held = TRUE, vacuous = FALSE, msg = "", caught = "")
  tryCatch(
    withCallingHandlers(
      code,
      warning = function(w) invokeRestart("muffleWarning"),
      message = function(m) invokeRestart("muffleMessage")
    ),
    error = function(e) {
      res$held <<- FALSE
      res$msg <<- brms_port_clean(conditionMessage(e))
    }
  )
  brms_mark(call, res$held, res$msg)
  if (nzchar(Sys.getenv("FRMTMB_BRMSPORT_RECORD"))) {
    brms_port_record("setup", id, "", res)
  }
  invisible(NULL)
}

# A setup line that reaches brms's own namespace. It is not run, and the
# names it would assign are stale, so nothing downstream reads brms's
# answer as frmtmb's.
brms_setup_not_run <- function(id, names) {
  msg <- "NOT RUN: the setup line reaches brms's own namespace"
  env <- parent.frame()
  st <- if (identical(brms_port_state$env, env)) brms_port_state$stale
  if (is.null(st)) st <- list()
  for (nm in strsplit(names, " ", fixed = TRUE)[[1]]) st[[nm]] <- msg
  brms_port_state$env <- env
  brms_port_state$stale <- st
  if (nzchar(Sys.getenv("FRMTMB_BRMSPORT_RECORD"))) {
    brms_port_record("setup", id, "", list(held = FALSE, vacuous = FALSE,
                                           msg = msg, caught = ""))
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------
# Shims. Each maps a brms entry point onto the frmtmb one that answers
# the same question, and nothing else: formulas, arguments and expected
# messages stay brms's. A generated file binds them at its top, so they
# never leak into another test file.
# ---------------------------------------------------------------------

# brms arguments that steer Stan, a backend or a file cache. frm() has no
# counterpart, and passing them on would make every call fail on the
# argument rather than on what brms's assertion is about.
brms_sampling_args <- c(
  "file", "file_refit", "backend", "chains", "iter", "warmup", "cores",
  "seed", "refresh", "silent", "testmode", "save_pars", "sample_prior",
  "algorithm", "stan_model_args", "mock_fit", "rename", "threads",
  "internal"
)

# Works on the unevaluated arguments, so a dropped save_pars() call is
# never evaluated (frmtmb has no save_pars()).
brms_forward <- function(fun, fixed, dots_call, env) {
  nm <- names(dots_call)
  if (is.null(nm)) nm <- rep("", length(dots_call))
  keep <- dots_call[!nm %in% brms_sampling_args]
  vals <- lapply(keep, eval, envir = env)
  do.call(fun, c(fixed, vals))
}

# brm() is frm(). brms refuses before sampling and frm() refuses before
# or while fitting, so the whole fit runs.
brms_shim_brm <- function(formula, data = NULL, family = NULL, ...) {
  brms_forward(frmtmb::frm,
               list(formula = formula, data = data, family = family),
               as.list(substitute(list(...)))[-1], parent.frame())
}

# standata() is frm(dry_run = "frame"), where frmtmb assembles and
# validates the response, the addition terms and the design, which is
# what brms's standata() does before it names things for Stan. The
# result is read through brms_standata_view(), which gives the frame's
# contents brms's Stan data names.
brms_shim_standata <- function(formula, data = NULL, family = NULL, ...) {
  fr <- brms_forward(frmtmb::frm,
                     list(formula = formula, data = data, family = family,
                          dry_run = "frame"),
                     as.list(substitute(list(...)))[-1], parent.frame())
  brms_standata_view(fr)
}

# brms names the intercept column "Intercept" because Stan identifiers
# cannot carry parentheses.
brms_x_names <- function(X) {
  cn <- colnames(X)
  if (!is.null(cn)) colnames(X)[cn == "(Intercept)"] <- "Intercept"
  X
}

brms_standata_view <- function(fr) {
  out <- list()
  resps <- names(fr$y)
  mv <- length(resps) > 1L
  sfx <- function(nm, r) if (mv) paste0(nm, "_", r) else nm
  for (r in resps) {
    y <- fr$y[[r]]
    out[[sfx("Y", r)]] <- as.array(y)
    out[[sfx("N", r)]] <- length(y)
    av <- fr$aterm_values[[r]]
    # addition terms: brms's names, the same per-row vectors
    # a constant term such as trials(10) is stored once; brms repeats it
    for (nm in intersect(names(av), c("se", "weights", "trials"))) {
      out[[sfx(nm, r)]] <- as.array(rep_len(av[[nm]], length(y)))
    }
    if (!is.null(av[["cens"]])) {
      out[[sfx("cens", r)]] <- as.array(av[["cens"]])
      if (!is.null(av[["cens_y2"]])) {
        # brms's rcens holds the upper bound on interval-censored rows
        # and 0 elsewhere
        out[[sfx("rcens", r)]] <-
          as.array(ifelse(av[["cens"]] == 2, av[["cens_y2"]], 0))
      }
    }
    for (lp in fr$linpreds) {
      if (!identical(lp$resp, r) || is.null(lp$X)) next
      nm <- if (identical(lp$dpar, "mu")) "X" else paste0("X_", lp$dpar)
      out[[sfx(nm, r)]] <- brms_x_names(lp$X)
    }
  }
  tau <- fr$par_template$tau_raw
  if (!is.null(tau)) out$nthres <- length(tau)
  out
}

brms_shim_rename_pars <- function(x) x

# brms:::validate_newdata() is what brms's prediction methods call on
# newdata; fitted() is frmtmb's prediction method that takes brms's
# arguments today (predict() moves under item 2.6d).
brms_shim_validate_newdata <- function(newdata, object, ...) {
  invisible(stats::fitted(object, newdata = newdata, ...))
}

# update() with brms's sampling switches removed, as brm() drops them.
brms_shim_update <- function(object, ...) {
  brms_forward(stats::update, list(object),
               as.list(substitute(list(...)))[-1], parent.frame())
}

# ---------------------------------------------------------------------
# Fixtures standing in for brms:::brmsfit_example1..6. The DATA are
# brms's own, read from its namespace, so every count the suite asserts
# (40 rows, 10 patients, 4 visits, 8 subjects) is the one brms's fits
# have. The FORMULAS are brms's except where frmtmb cannot express a
# term, and `changed` says what moved and why; an assertion that reaches
# a changed part is classified in the ledger, never loosened here.
# ---------------------------------------------------------------------
brms_fixture_data <- function(k) {
  fit <- get(paste0("brmsfit_example", k), envir = asNamespace("brms"))
  d <- as.data.frame(fit$data)
  attr(d, "terms") <- NULL
  d
}

brms_fixture_spec <- function(k) {
  switch(
    k,
    list(
      brms = paste("count ~ Trt * Age + mo(Exp) + s(Age) + volume +",
                   "offset(Age) + (1 + Trt | visit) + arma(visit, patient),",
                   "sigma ~ Trt; student"),
      formula = bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume +
                     offset(Age) + (1 + Trt | visit) +
                     arma(visit, patient, cov = TRUE),
                   sigma ~ Trt),
      family = student(),
      changed = paste("arma() gains cov = TRUE: frmtmb has only the",
                      "residual-covariance ARMA (?frmtmb-autocor)")),
    list(
      brms = paste("count | weights(AgeSD) ~ 1/(1 + exp(-a)) *",
                   "exp(b * Trt), a ~ Age + (1 | ID1 | patient),",
                   "b ~ Age + (1 | ID1 | patient); gamma(identity)"),
      formula = bf(count | weights(AgeSD) ~ 1 / (1 + exp(-a)) *
                     exp(b * Trt),
                   a ~ Age + (1 | ID1 | patient),
                   b ~ Age + (1 | ID1 | patient), nl = TRUE),
      family = Gamma("identity"),
      data = function(d) {
        d$Trt <- as.numeric(as.character(d$Trt))
        d
      },
      changed = paste("Trt enters the nonlinear body as its numeric codes",
                      "0 and 1, the covariate brms's own Stan data C_1",
                      "carries; frm() dies on the factor with an RTMB",
                      "advector error (dev/brmsport-probe-nl.R)")),
    list(
      brms = paste("count ~ Trt * me(Age, AgeSD) +",
                   "(1 + mmc(Age, volume) | mm(patient, visit)); gaussian"),
      formula = bf(count ~ Trt * Age +
                     (1 + mmc(Age, volume) | mm(patient, visit))),
      family = gaussian(),
      changed = paste("me(Age, AgeSD) becomes Age: frmtmb's me()",
                      "takes numeric interaction multipliers only, and",
                      "brms's Trt is a factor (dev/me-findings.md)")),
    list(
      brms = "rating ~ x1 + cs(x2) + (cs(x2) || subject), disc ~ 1; sratio",
      formula = bf(rating ~ x1 + cs(x2) + (1 + x2 || subject)),
      family = sratio(),
      changed = paste("disc ~ 1 is dropped, since frmtmb's ordinal",
                      "families have no disc parameter; (cs(x2) || subject)",
                      "becomes (1 + x2 || subject), since cs() is not",
                      "available inside a group-level term")),
    list(
      brms = paste("count ~ Age + (1 | gr(patient, by = gender)),",
                   "mu2 ~ Age; mixture(gaussian, exponential)"),
      formula = bf(count ~ Age + (1 | patient), mu2 ~ Age),
      family = mixture(gaussian(), exponential()),
      changed = paste("gr(patient, by = gender) becomes patient: frmtmb's",
                      "gr() takes no by (dev/brms-suite-audit.md",
                      "section 5)")),
    list(
      brms = paste("volume ~ Trt + gp(Age, by = Trt, gr = TRUE); gaussian",
                   "and count ~ Trt + Age; poisson; no rescor"),
      formula = bf(volume ~ Trt + gp(Age)) +
        gaussian() + bf(count ~ Trt + Age) + poisson() + set_rescor(FALSE),
      family = NULL,
      changed = paste("gp(Age, by = Trt, gr = TRUE) becomes gp(Age):",
                      "frmtmb's gp() takes no by or gr"))
  )
}

# One fit per fixture per process: the six are shared by every block of
# a file, as brms's are.
brms_fixture <- function(k) {
  key <- paste0("fit", k)
  hit <- brms_port_state$fixtures[[key]]
  if (!is.null(hit)) return(hit)
  spec <- brms_fixture_spec(k)
  warns <- character()
  d <- brms_fixture_data(k)
  if (!is.null(spec$data)) d <- spec$data(d)
  # update() re-evaluates the stored call, so the call must carry the
  # formula and data themselves rather than names local to this function
  args <- list(formula = spec$formula, data = d)
  if (!is.null(spec$family)) args$family <- spec$family
  fit <- withCallingHandlers(
    do.call(frmtmb::frm, args),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  attr(fit, "brmsport_warnings") <- warns
  if (is.null(brms_port_state$fixtures)) brms_port_state$fixtures <- list()
  brms_port_state$fixtures[[key]] <- fit
  fit
}

