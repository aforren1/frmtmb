# A fit can outlive the session that made it. saveRDS() stores a
# family's closures with a reference to the namespace that built them,
# and readRDS() in a session where that package is not installed puts
# the global environment in its place, so every internal the closures
# call is then "not found", a message that names a function of another
# package and not the package to install. These helpers stop first, by
# the package's name, at the points where core calls a model's own
# functions: frame assembly (frm(), update(), refits), the objective
# (every retape, which summary(), vcov() and confint() reach through
# sdr_of()), the response mean, the simulators and the curve functions.
#
# They are written for any part of a model an extension supplies, not
# for one family: `require_extension()` takes a package and a phrase,
# and `require_model_packages()` reads the packages off the families of
# a fit, frame or spec. A covariance structure that moves to an
# extension can call `require_extension()` with its own phrase.
#
# A check fires only for the part it guards, so a method that needs
# nothing of an absent package still answers. Once one fires, the
# refusal names every part of the model this session cannot use
# (families, covariance structures, residual correlation terms): a fit
# with a family from one extension and a structure from another named
# one package at a time before, so the user found the second package
# only after installing the first.

#' Whether `pkg` is core, or an extension this session can load.
#'
#' @noRd
extension_usable <- function(pkg) {
  is.null(pkg) || !nzchar(pkg) || identical(pkg, "frmtmb") ||
    isNamespaceLoaded(pkg) || requireNamespace(pkg, quietly = TRUE)
}

#' An error handler for `tryCatch()` that falls back to `value`, or
#' calls it with the condition when it is a function, except for the
#' refusal for a part of the model this session cannot use, which it
#' signals again.
#'
#' A handler that reads any error as "no covariance" or "no answer"
#' turned that refusal into an NA standard error (`fixef()`) or a
#' second, wrong message (`predict()`, `influence()`). Every handler
#' whose expression can reach the refusal uses this one
#' (`dev/smallfix-handlers.R` lists them).
#'
#' @noRd
keep_missing <- function(value = NULL) {
  force(value)
  function(e) {
    if (inherits(e, "frmtmb_extension_missing")) frm_stop(e)
    if (is.function(value)) value(e) else value
  }
}

#' Fit-level checks for a method of a fit: every family and every
#' structure, so a refusal names every part this session cannot use
#' before a method reaches one part on its own.
#'
#' @noRd
require_fit_packages <- function(fit) {
  require_model_packages(fit)
  require_structure_packages(fit)
}

#' Stop, naming `pkg`, when it is neither core nor loadable here.
#'
#' `what` completes "<what> is from the extension package <pkg>", for
#' example "The cox family of this model".
#'
#' @noRd
require_extension <- function(pkg, what) {
  if (extension_usable(pkg)) return(invisible(TRUE))
  stop_extension_parts(list(list(kind = "missing", pkg = pkg, what = what)))
}

#' The packages a family's own functions come from: the package that
#' built it (`frm_family_package()`), and for a mixture, the packages
#' that built its components.
#'
#' @noRd
family_packages <- function(fam) {
  unique(c(frm_family_package(fam),
           attr(fam, "frmtmb_component_packages", exact = TRUE)))
}

#' The extension that now holds a family an older frmtmb built, or
#' `NULL`. Such a family has a name in `extension_family_homes`, records
#' frmtmb as its package, and has its density in frmtmb's namespace: a
#' `cox()` fit saved by frmtmb 0.73.0 or earlier. A family the extension
#' built records the extension, and a user's family of the same name has
#' its density elsewhere, so neither is matched; a `wiener()` fit an
#' older frmtmb.eam built records frmtmb.eam and still works.
#'
#' @noRd
old_core_family_home <- function(fam) {
  nm <- fam[["family"]]
  if (!is.character(nm) || length(nm) != 1L ||
        !nm %in% names(extension_family_homes) ||
        !identical(frm_family_package(fam), "frmtmb")) {
    return(NULL)
  }
  f <- fam[["lpdf"]]
  if (!is.function(f) ||
        !identical(topenv(environment(f)), asNamespace("frmtmb"))) {
    return(NULL)
  }
  unname(extension_family_homes[[nm]])
}

#' The parts of a fit, a frame, a spec or a single family that this
#' session cannot use because of a family. One list per part: `kind`
#' ("old" for a family an older frmtmb built, "missing" for a package
#' that cannot be loaded), `pkg` and `what`.
#'
#' @noRd
family_missing_parts <- function(x) {
  fams <- if (inherits(x, "frmtmb_family")) {
    list(x)
  } else {
    spec <- if (!is.null(x[["spec"]])) x[["spec"]] else x
    lapply(spec[["responses"]], function(r) r[["family"]])
  }
  out <- list()
  for (fam in fams) {
    if (!inherits(fam, "frmtmb_family")) next
    what <- paste0("The ", fam[["family"]], " family of this model")
    home <- old_core_family_home(fam)
    if (!is.null(home)) {
      out[[length(out) + 1L]] <- list(kind = "old", pkg = home, what = what)
      next
    }
    for (p in family_packages(fam)) {
      if (extension_usable(p)) next
      out[[length(out) + 1L]] <- list(kind = "missing", pkg = p, what = what)
    }
  }
  out
}

#' The same for the random-effect blocks of a fit or a frame. A block's
#' structure must be registered in this session, which loading the
#' package that owns it does.
#'
#' The owner is what the block recorded when the frame was built
#' (`frmtmb_register_covstruct(package =)`), or for a structure core
#' lists as an extension's (`extension_special_homes`), that extension.
#' A structure registered in a session rather than by a package has no
#' package to load ("session"). One whose package is loaded but does not
#' register it is from another version of that package ("stale").
#'
#' @noRd
structure_missing_parts <- function(x) {
  frame <- if (!is.null(x[["frame"]])) x[["frame"]] else x
  out <- list()
  for (bk in frame[["re_blocks"]] %||% list()) {
    nm <- bk[["covstruct"]]
    if (!is.null(frmtmb_covstruct_registry$reg[[nm]])) next
    pkg <- bk[["covstruct_package"]] %||% extension_special_home(nm)
    kind <- if (is.null(pkg)) {
      "session"
    } else if (!extension_usable(pkg)) {
      "missing"
    } else if (is.null(frmtmb_covstruct_registry$reg[[nm]])) {
      # loading the package registers what it owns, so a name still
      # absent after the load is not this version's
      "stale"
    } else {
      next
    }
    out[[length(out) + 1L]] <- list(
      kind = kind, pkg = pkg,
      what = paste0("The ", nm, "() structure of this model"))
  }
  out
}

#' The parts of a fit or a frame that this session cannot use because
#' of a covariance structure or a residual correlation term
#' (`autocor_missing_parts()`, R/autocor-seam.R): both are seams an
#' extension registers, and a refusal names them with the families.
#'
#' @noRd
seam_missing_parts <- function(x) {
  c(structure_missing_parts(x), autocor_missing_parts(x))
}

#' Check every family of a fit, a frame, a spec or a single family.
#'
#' @noRd
require_model_packages <- function(x) {
  parts <- family_missing_parts(x)
  if (!length(parts)) return(invisible(TRUE))
  stop_extension_parts(c(parts, seam_missing_parts(x)))
}

#' Check every random-effect block and every residual correlation term
#' of a fit or a frame. Called where the families are checked (the
#' objective, the standard-error machinery) and at the head of the
#' methods that read a block's entry directly.
#'
#' @noRd
require_structure_packages <- function(x) {
  parts <- seam_missing_parts(x)
  if (!length(parts)) return(invisible(TRUE))
  stop_extension_parts(c(family_missing_parts(x), parts))
}

#' One refusal for every part of a model this session cannot use, and
#' one instruction that installs and loads every package they need. A
#' single part reads as each per-part refusal did before.
#'
#' @noRd
stop_extension_parts <- function(parts) {
  key <- vapply(parts, function(p) paste(p$kind, p$pkg %||% "", p$what), "")
  parts <- parts[!duplicated(key)]
  of <- function(k) Filter(function(p) identical(p$kind, k), parts)
  pkgs <- function(ps) unique(vapply(ps, function(p) p$pkg, ""))
  low <- function(s) sub("^The ", "the ", s)
  old <- of("old")
  miss <- of("missing")
  s <- character(0)
  for (p in old) {
    s <- c(s, paste0(p$what, " was built by an older frmtmb, when the ",
                     "family was in frmtmb. It is in the extension package ",
                     p$pkg, " now, and the functions this fit calls are no ",
                     "longer in frmtmb."))
  }
  mp <- pkgs(miss)
  if (length(mp) == 1L) {
    w <- vapply(miss, function(p) p$what, "")
    s <- c(s, paste0(comma_and(c(w[1L], low(w[-1L]))),
                     if (length(w) > 1L) " are" else " is",
                     " from the extension package ", mp, ", which is not ",
                     "installed or cannot be loaded in this session, so its ",
                     "functions are not available."))
  } else if (length(mp) > 1L) {
    per <- vapply(mp, function(m) {
      w <- vapply(Filter(function(p) identical(p$pkg, m), miss),
                  function(p) low(p$what), "")
      paste0(m, " (", comma_and(w), ")")
    }, "")
    s <- c(s, paste0("This model needs ", length(mp), " extension packages ",
                     "that are not installed or cannot be loaded in this ",
                     "session, so their functions are not available: ",
                     comma_and(per), "."))
  }
  for (p in of("session")) {
    s <- c(s, paste0(p$what, " was registered in an R session, not by a ",
                     "package, and is not registered in this one. Register ",
                     "it again with ",
                     p[["register"]] %||% "frmtmb_register_covstruct()",
                     " before using the fit."))
  }
  for (p in of("stale")) {
    s <- c(s, paste0(p$what, " is from the extension package ", p$pkg,
                     ", which is loaded but does not register it. Install ",
                     "the version of ", p$pkg, " the fit was made with."))
  }
  # An old family's home is asked whether it is installed, and whether
  # the installation is whole (an installed package has
  # Meta/package.rds), not loaded: loading it here would register its
  # families in a session that only read a fit. A home installed but
  # not whole cannot be loaded, so it is to be reinstalled, not loaded
  homes <- pkgs(old)
  state <- vapply(homes, function(h) {
    if (isNamespaceLoaded(h)) return("ok")
    d <- system.file(package = h)
    if (!nzchar(d)) "absent" else if (!file.exists(file.path(d, "Meta",
      "package.rds"))) "broken" else "ok"
  }, "")
  broken <- homes[state == "broken"]
  if (length(broken)) {
    s <- c(s, paste0(comma_and(broken), if (length(broken) > 1L) " are" else
      " is", " installed but cannot be loaded in this session."))
  }
  # the install list and the library() list in one order
  order_ <- unique(c(homes, mp))
  inst <- order_[order_ %in% c(mp, homes[state == "absent"])]
  reinst <- order_[order_ %in% broken]
  libs <- paste0("library(", order_, ")", collapse = "; ")
  # "it" only where the sentence naming the package is the one before
  more <- length(of("session")) + length(of("stale")) > 0L
  if (length(homes)) {
    verb <- c(if (length(inst)) paste0("install ", comma_and(inst)),
              if (length(reinst)) paste0("reinstall ", comma_and(reinst)))
    s <- c(s, if (length(verb)) {
      paste0(toupper(substring(verb[1L], 1L, 1L)),
             substring(paste(verb, collapse = " and "), 2L),
             " and refit the model with ", libs, " loaded")
    } else paste0("Refit the model with ", libs, " loaded"))
  } else if (length(mp)) {
    s <- c(s, paste0("Install ", if (more) comma_and(mp) else
      if (length(mp) > 1L) "them" else "it", " and run ", libs))
  } else {
    s[length(s)] <- sub("[.]$", "", s[length(s)])
  }
  msg <- paste(s, collapse = " ")
  frm_stop(msg, class = "frmtmb_extension_missing", call. = FALSE)
}
