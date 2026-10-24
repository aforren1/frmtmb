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

#' Stop, naming `pkg`, when it is neither core nor loadable here.
#'
#' `what` completes "<what> is from the extension package <pkg>", for
#' example "The cox family of this model".
#'
#' @noRd
require_extension <- function(pkg, what) {
  if (is.null(pkg) || !nzchar(pkg) || identical(pkg, "frmtmb") ||
        isNamespaceLoaded(pkg) ||
        requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }
  frm_stop(what, " is from the extension package ", pkg, ", which is not ",
           "installed or cannot be loaded in this session, so its ",
           "functions are not available. Install it and run library(",
           pkg, ")", class = "frmtmb_extension_missing", call. = FALSE)
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

#' Check every family of a fit, a frame, a spec or a single family with
#' `require_extension()`.
#'
#' @noRd
require_model_packages <- function(x) {
  fams <- if (inherits(x, "frmtmb_family")) {
    list(x)
  } else {
    spec <- if (!is.null(x[["spec"]])) x[["spec"]] else x
    lapply(spec[["responses"]], function(r) r[["family"]])
  }
  for (fam in fams) {
    if (!inherits(fam, "frmtmb_family")) next
    home <- old_core_family_home(fam)
    if (!is.null(home)) {
      frm_stop("The ", fam[["family"]], " family of this model was built ",
               "by an older frmtmb, when the family was in frmtmb. It is ",
               "in the extension package ", home, " now, and the ",
               "functions this fit calls are no longer in frmtmb. Refit ",
               "the model with library(", home, ") loaded",
               class = "frmtmb_extension_missing", call. = FALSE)
    }
    for (p in family_packages(fam)) {
      require_extension(p, paste0("The ", fam[["family"]],
                                  " family of this model"))
    }
  }
  invisible(TRUE)
}

#' Check every random-effect block of a fit or a frame with
#' `require_extension()`: the block's structure must be registered in
#' this session, which loading the package that owns it does.
#'
#' The owner is what the block recorded when the frame was built
#' (`frmtmb_register_covstruct(package =)`), or for a structure core
#' lists as an extension's (`extension_special_homes`), that extension.
#' A structure registered in a session rather than by a package has no
#' package to load, and is refused with that said. Called where the
#' families are checked (the objective, the standard-error machinery)
#' and at the head of the methods that read a block's entry directly.
#'
#' @noRd
require_structure_packages <- function(x) {
  frame <- if (!is.null(x[["frame"]])) x[["frame"]] else x
  for (bk in frame[["re_blocks"]] %||% list()) {
    nm <- bk[["covstruct"]]
    if (!is.null(frmtmb_covstruct_registry$reg[[nm]])) next
    pkg <- bk[["covstruct_package"]] %||% extension_special_home(nm)
    what <- paste0("The ", nm, "() structure of this model")
    require_extension(pkg, what)
    if (!is.null(frmtmb_covstruct_registry$reg[[nm]])) next
    if (is.null(pkg)) {
      frm_stop(what, " was registered in an R session, not by a package, ",
               "and is not registered in this one. Register it again with ",
               "frmtmb_register_covstruct() before using the fit",
               class = "frmtmb_extension_missing", call. = FALSE)
    }
    frm_stop(what, " is from the extension package ", pkg, ", which is ",
             "loaded but does not register it. Install the version of ",
             pkg, " the fit was made with", class = "frmtmb_extension_missing",
             call. = FALSE)
  }
  invisible(TRUE)
}
