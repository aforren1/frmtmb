# Classed conditions: the brms contract that every refusal can be caught
# by class.
#
# brms raises every refusal through stop2(), which gives the condition
# the class `brms_error`, so `tryCatch(brms_error = )` catches any of
# them. frmtmb and its extensions raise through these three functions
# instead of stop(), warning() and message(), for the same reason. They
# take the same arguments as the base functions they replace, so a call
# site changes its name and nothing else, and the message text and the
# recorded call are what the base function would have produced.
#
# dev/conditions-findings.md records why the conditions are built with
# base R rather than rlang, and tests/testthat/test-conditions.R keeps a
# bare stop(), warning() or message() from coming back.

#' Classed conditions
#'
#' @description
#' Every error, warning and message that frmtmb raises itself carries a
#' class that names frmtmb, so a caller can catch all of them without
#' matching message text:
#'
#' | raised as | class vector |
#' |---|---|
#' | error | `c("frmtmb_error", "error", "condition")` |
#' | warning | `c("frmtmb_warning", "warning", "condition")` |
#' | message | `c("frmtmb_message", "message", "condition")` |
#'
#' A condition that an extension package raises has one more class in
#' front, made from the package name with each dot changed to an
#' underscore. For example, an error from frmtmb.eam has the class
#' `c("frmtmb_eam_error", "frmtmb_error", "error", "condition")`. Thus
#' `tryCatch(frmtmb_error = )` catches a refusal from frmtmb and from
#' every extension, and `tryCatch(frmtmb_eam_error = )` catches the
#' refusals that are frmtmb.eam's by the rules below.
#'
#' The subclass names one package, found in this order:
#' \enumerate{
#'   \item The package that the raising function gives as `package =`.
#'     frmtmb gives it in two cases. The first is a refusal whose text an
#'     extension gave as data, such as the `refusals` of
#'     [frmtmb_structure()] or the `sim_refusal` of [frmtmb_family()].
#'     The second is a refusal about a declaration of a family: an
#'     addition term that the family needs or does not read, or a
#'     simulator, mean, variance, unit deviance, latent state or
#'     likelihood factor that the family does not supply. In both cases
#'     the package is the one that built the family or the structure.
#'     frmtmb.eam gives it for a refusal of its non-decision-time seam
#'     (`frmtmb.eam::ndt_bound()` and the functions near it), where it
#'     is the package that called the seam. Thus a refusal about the
#'     bound of `frmtmb.learn::rlddm()` is a `frmtmb_learn_error`.
#'   \item Otherwise, the package of the function that raises the
#'     condition.
#' }
#' A family or a structure belongs to the package whose function called
#' [frmtmb_family()] or [frmtmb_structure()] to build it. No subclass is
#' added for frmtmb itself or for code outside a package. So a refusal
#' that frmtmb raises for any other reason has no subclass, also when it
#' is about a model with a family from an extension. Examples are a
#' refusal of an argument of [predict()], of a prior, of `newdata`, or
#' of a value outside the range of a link function. An extension
#' function that the user calls, such as `frmtmb.eam::ndt_time()` on a
#' fit of another package, raises with its own subclass.
#'
#' This is the contract that brms keeps with its class `brms_error`.
#' frmtmb builds the conditions with base R, so the class vector does
#' not include `rlang_error`.
#'
#' Some conditions do not have these classes:
#' \itemize{
#'   \item An error from another package, such as RTMB, TMB, Matrix or
#'     base R, keeps its own class. The exception is an error that
#'     occurs while [frm()] optimizes a model: frmtmb catches it and
#'     raises it again as a `frmtmb_fit_error`, which is also a
#'     `frmtmb_error`.
#'   \item An argument that the function does not have, such as
#'     `frm(priors = )`, is refused by R itself with
#'     "unused argument".
#'   \item The refusal that `model.matrix()` gives for a multivariate
#'     fit also has the class `simpleError`, after `frmtmb_error`. The
#'     insight package examines that class to replace a failed standard
#'     error with `NA`.
#' }
#'
#' @section For extension authors:
#' `frm_stop()`, `frm_warning()` and `frm_message()` are the functions
#' that frmtmb uses to raise its conditions. Use them in an extension
#' in place of [stop()], [warning()] and [message()]. They take the
#' same arguments as the base functions and give the same message text
#' and the same call, and they add the classes in the table above. The
#' class of the extension is found from the namespace of the function
#' that calls them, so an extension does not supply it. A function of
#' base R between the two, such as [lapply()] when a helper is given as
#' its `FUN`, is skipped.
#'
#' Given one condition object, the helpers signal that object, as the
#' base functions do, with the frmtmb classes added in front of its own
#' classes.
#'
#' `class` adds classes in front of the others, for a condition that a
#' caller must be able to catch by itself. `package` names the package
#' whose subclass the condition gets, in place of the package of the
#' calling function. Use it where one package raises a refusal that
#' another package wrote. `package = "frmtmb"` gives no subclass.
#' `frm_family_package()` gives the package that built a family, the
#' value to use for a refusal about that family.
#'
#' `frm_match_arg()` does what [match.arg()] does, and a value that
#' matches no choice is a `frmtmb_error` that names the argument, the
#' value and the choices. The subclass is that of the function that
#' calls `frm_match_arg()`.
#'
#' These differences from the base functions remain:
#' \itemize{
#'   \item With `call. = TRUE`, an error that nothing catches in a
#'     session that is not interactive can print one more line,
#'     `Calls:`, that includes the name of the helper. This line is R's
#'     short traceback. The message and the call are the same as from
#'     [stop()]. frmtmb itself uses `call. = FALSE`.
#'   \item `noBreaks.` has no effect.
#'   \item `frm_message()` records the call to itself, where
#'     [message()] records the call to [message()]. R does not print the
#'     call of a message.
#' }
#'
#' @param ... Zero or more objects that are pasted together, as in
#'   [stop()], or one condition object.
#' @param call. Logical. If `TRUE`, the call of the function that raises
#'   the condition is part of the condition, as in [stop()].
#' @param immediate. Logical. If `TRUE` and `getOption("warn")` is 0,
#'   the warning prints at once, as in [warning()].
#' @param noBreaks. Ignored. It is accepted so that a call to
#'   [warning()] can use the same arguments.
#' @param domain Passed to [.makeMessage()], for translation.
#' @param appendLF Logical. If `TRUE`, a newline is added to the message,
#'   as in [message()].
#' @param class A character vector of more classes, put in front of the
#'   classes that frmtmb adds.
#' @param package `NULL`, or the name of the package whose subclass the
#'   condition gets.
#' @param arg,choices,several.ok As in [match.arg()].
#' @param family A family object from [frmtmb_family()] or
#'   [frmtmb_structure()].
#' @return `frm_stop()` does not return. `frm_warning()` returns its
#'   message invisibly. `frm_message()` returns `NULL` invisibly.
#'   `frm_match_arg()` returns the matched choices. `frm_family_package()`
#'   returns the name of the package that built `family`, or `"frmtmb"`.
#' @seealso [frmtmb-extension-api] for the rest of the interface an
#'   extension uses.
#' @examples
#' d <- data.frame(y = c(1, 2, 3))
#' e <- tryCatch(frm(y ~ 1, data = d, family = "not_a_family"),
#'               frmtmb_error = function(e) e)
#' class(e)
#' conditionMessage(e)
#'
#' f <- function() frm_warning("a warning from f()")
#' w <- tryCatch(f(), frmtmb_warning = function(w) w)
#' class(w)
#' conditionCall(w)
#'
#' g <- function(type = c("response", "link")) frm_match_arg(type)
#' g("resp")
#' e <- tryCatch(g("bogus"), frmtmb_error = function(e) e)
#' conditionMessage(e)
#' @name frmtmb-conditions
#' @aliases frm_stop frm_warning frm_message frm_match_arg frm_family_package
NULL

#' The package that owns the code evaluated in `env`, or `NULL` for code
#' outside a namespace.
#'
#' base R's own closures are skipped, because a helper passed as a value
#' (`lapply(x, frm_stop)`) is called from inside `lapply()`, whose frame
#' belongs to base, while the refusal belongs to the function that
#' passed it. The search only runs on that path; a direct call stops at
#' the first `topenv()`.
#'
#' @noRd
frm_env_package <- function(env) {
  top <- topenv(env)
  if (identical(top, .BaseNamespaceEnv)) {
    frames <- sys.frames()
    parents <- sys.parents()
    hit <- which(vapply(frames, identical, NA, env))
    n <- if (length(hit)) hit[length(hit)] else 0L
    while (n > 0L && identical(topenv(frames[[n]]), .BaseNamespaceEnv)) {
      n <- parents[n]
    }
    top <- if (n > 0L) topenv(frames[[n]]) else globalenv()
  }
  if (isNamespace(top)) unname(getNamespaceName(top))
}

#' The frmtmb classes of a condition of one kind. The package is found
#' from where the raising function was defined rather than passed at
#' each call site, so that a call site stays a rename of the base call;
#' `package` overrides it where one package raises another's text.
#'
#' @noRd
frm_condition_class <- function(kind, env, package = NULL) {
  pkg <- if (is.null(package)) frm_env_package(env) else package
  own <- if (length(pkg) && !identical(pkg, "frmtmb")) {
    paste0(gsub(".", "_", pkg, fixed = TRUE), "_", kind)
  }
  c(own, paste0("frmtmb_", kind), kind, "condition")
}

#' A condition object given to a helper keeps its own message, call and
#' classes, as it does under the base function, and gains the frmtmb
#' classes in front of them.
#'
#' @noRd
frm_reclass <- function(cnd, kind, class, env, package) {
  add <- c(class, frm_condition_class(kind, env, package))
  # the last two are `kind` and "condition", which the object has
  add <- add[seq_len(length(add) - 2L)]
  class(cnd) <- c(setdiff(add, class(cnd)), class(cnd))
  cnd
}

#' Whether `...` is one condition object, the shape base R signals as is.
#'
#' @noRd
frm_is_condition_arg <- function(args) {
  length(args) == 1L && inherits(args[[1L]], "condition")
}

#' @rdname frmtmb-conditions
#' @export
frm_stop <- function(..., call. = TRUE, domain = NULL, class = NULL,
                     package = NULL) {
  args <- list(...)
  cnd <- if (frm_is_condition_arg(args)) {
    frm_reclass(args[[1L]], "error", class, parent.frame(), package)
  } else {
    # stop() records the call of the function that called it; one frame
    # further up from here is that same function
    structure(
      class = c(class,
                frm_condition_class("error", parent.frame(), package)),
      list(message = .makeMessage(..., domain = domain),
           call = if (call.) sys.call(-1L)))
  }
  stop(cnd)
}

#' @rdname frmtmb-conditions
#' @export
frm_warning <- function(..., call. = TRUE, immediate. = FALSE,
                        noBreaks. = FALSE, domain = NULL, class = NULL,
                        package = NULL) {
  args <- list(...)
  cnd <- if (frm_is_condition_arg(args)) {
    frm_reclass(args[[1L]], "warning", class, parent.frame(), package)
  } else {
    structure(
      class = c(class,
                frm_condition_class("warning", parent.frame(), package)),
      list(message = .makeMessage(..., domain = domain),
           call = if (call.) sys.call(-1L)))
  }
  # warning() ignores `immediate.` for a condition object. With warn = 0
  # an immediate warning prints the way warn = 1 prints every warning, so
  # that setting is held for the one signal.
  if (isTRUE(immediate.) && identical(as.numeric(getOption("warn", 0)), 0)) {
    old <- options(warn = 1)
    on.exit(options(old), add = TRUE)
  }
  warning(cnd)
}

#' @rdname frmtmb-conditions
#' @export
frm_message <- function(..., domain = NULL, appendLF = TRUE,
                        class = NULL, package = NULL) {
  args <- list(...)
  cnd <- if (frm_is_condition_arg(args)) {
    frm_reclass(args[[1L]], "message", class, parent.frame(), package)
  } else {
    # message() records its own call rather than its caller's, so this
    # does the same
    structure(
      class = c(class,
                frm_condition_class("message", parent.frame(), package)),
      list(message = .makeMessage(..., domain = domain,
                                  appendLF = appendLF),
           call = sys.call()))
  }
  message(cnd)
}

#' @rdname frmtmb-conditions
#' @export
frm_family_package <- function(family) {
  attr(family, "frmtmb_package", exact = TRUE) %||% "frmtmb"
}

#' @rdname frmtmb-conditions
#' @export
frm_match_arg <- function(arg, choices, several.ok = FALSE) {
  # match.arg()'s own reading of the caller: its formals give the
  # choices, and its call is the one the refusal records
  sysp <- sys.parent()
  caller <- parent.frame()
  if (missing(choices)) {
    formal_args <- formals(sys.function(sysp))
    choices <- eval(formal_args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysp))
  }
  if (is.null(arg)) return(choices[1L])
  name <- deparse1(substitute(arg))
  refuse <- function(what) {
    frm_stop(errorCondition(paste0("`", name, "` ", what),
                            call = sys.call(sysp)),
             package = frm_env_package(caller) %||% "frmtmb")
  }
  one_of <- paste0("one of ", paste0("\"", choices, "\"", collapse = ", "))
  if (!is.character(arg)) {
    refuse(paste0("must be ", one_of, ", not ", arg_desc(arg)))
  }
  if (!several.ok) {
    if (identical(arg, choices)) return(arg[1L])
    if (length(arg) > 1L) {
      refuse(paste0("must be ", one_of, ", not ", arg_desc(arg)))
    }
  } else if (!length(arg)) {
    refuse(paste0("must name at least ", one_of))
  }
  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L)) {
    refuse(paste0("must be ", one_of, ", not ", arg_desc(arg)))
  }
  choices[i[i > 0L]]
}
