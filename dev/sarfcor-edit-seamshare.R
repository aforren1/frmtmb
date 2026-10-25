# One-off source edit for lane sarfcor, part A: the checks the two
# registration seams share become helpers with one message template
# each (test-message-uniqueness.R). Run on the worktree and on the
# stage-A copy: Rscript dev/sarfcor-edit-seamshare.R <root>
root <- commandArgs(TRUE)[1]
edit <- function(p, pairs) {
  p <- file.path(root, p)
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
pkg_old <- function(fn) paste0(
"  arg <- function(nm) paste0(\"", fn, "(\", nm, \" =)\")
  if (identical(package, \"frmtmb\")) package <- NULL
  if (!is.null(package) && (!is.character(package) ||
                              length(package) != 1L || is.na(package) ||
                              !nzchar(package))) {
    frm_stop(arg(\"package\"), \" must be one package name, or NULL\",
             call. = FALSE)
  }")
pkg_new <- function(fn) paste0(
"  arg <- function(nm) paste0(\"", fn, "(\", nm, \" =)\")
  package <- seam_check_package(package, arg(\"package\"))")
name_old <- "  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(name) || !identical(name, make.names(name))) {
    frm_stop(arg(\"name\"), \" must be one syntactic name, written as it \",
             \"appears in a formula and without parentheses\", call. = FALSE)
  }"
name_new <- "  seam_check_name(name, arg(\"name\"))"
fn_old <- function(null_ok) paste0(
"  fn_or <- function(x, nm, argn, null_ok = ", null_ok, ") {
    if (is.null(x) && null_ok) return(NULL)
    if (!is.function(x) || length(formals(x)) < argn) {
      frm_stop(arg(nm), \" must be a function of \", argn, \" argument\",
               if (argn > 1L) \"s\", call. = FALSE)
    }
    x
  }")
fn_new <- function(null_ok) paste0(
"  fn_or <- function(x, nm, argn, null_ok = ", null_ok, ") {
    seam_fn(x, arg(nm), argn, null_ok)
  }")
edit("R/covstruct-seam.R", list(
  c(pkg_old("frmtmb_register_covstruct"),
    pkg_new("frmtmb_register_covstruct")),
  c(name_old, name_new),
  c(fn_old("TRUE"), fn_new("TRUE")),
  c("#' The package a registration came from, or `NULL` for a call from",
"#' The checks the registration seams share, one message template each:
#' the owning package, the name as a formula writes it, and a slot that
#' must be a function of so many arguments.
#'
#' @noRd
seam_check_package <- function(package, what) {
  if (identical(package, \"frmtmb\")) return(NULL)
  if (!is.null(package) && (!is.character(package) ||
                              length(package) != 1L || is.na(package) ||
                              !nzchar(package))) {
    frm_stop(what, \" must be one package name, or NULL\", call. = FALSE)
  }
  package
}

#' @noRd
seam_check_name <- function(name, what) {
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(name) || !identical(name, make.names(name))) {
    frm_stop(what, \" must be one syntactic name, written as it \",
             \"appears in a formula and without parentheses\", call. = FALSE)
  }
  invisible(name)
}

#' @noRd
seam_fn <- function(x, what, argn, null_ok) {
  if (is.null(x) && null_ok) return(NULL)
  if (!is.function(x) || length(formals(x)) < argn) {
    frm_stop(what, \" must be a function of \", argn, \" argument\",
             if (argn > 1L) \"s\", call. = FALSE)
  }
  x
}

#' The package a registration came from, or `NULL` for a call from")))
edit("R/autocor-seam.R", list(
  c(pkg_old("frmtmb_register_autocor"), pkg_new("frmtmb_register_autocor")),
  c(name_old, name_new),
  c(fn_old("FALSE"), fn_new("FALSE"))))
cat("edited", root, "\n")
