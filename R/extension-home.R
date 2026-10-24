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
    for (p in family_packages(fam)) {
      require_extension(p, paste0("The ", fam[["family"]],
                                  " family of this model"))
    }
  }
  invisible(TRUE)
}
