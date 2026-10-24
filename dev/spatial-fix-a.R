# Fix round 1, part A edits to R/covstruct-seam.R (A1, A3, A4, A5, A8),
# applied to the tree named as the argument (the worktree, or the
# reconstructed AB stage tree), so that the two stay in step.
root <- commandArgs(trailingOnly = TRUE)[1]
rep_file <- function(p, old, new) {
  p <- file.path(root, p)
  s <- paste(readLines(p, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, s, fixed = TRUE)[[1]])
  if (regexpr(old, s, fixed = TRUE) < 0 || n != 1) {
    stop(p, ": count ", n, ": ", substr(old, 1, 70))
  }
  writeLines(sub(old, new, s, fixed = TRUE), p)
}
p <- "R/covstruct-seam.R"
rep_file(p, '#\' The names a registered structure may not take: the core\'s own
#\' formula vocabulary and its internal structure names. A structure is a
#\' formula special, so taking one of these would give one spelling two
#\' meanings in one formula. `car` and `spde` are brms and frmtmb
#\' spellings that an extension supplies, so they are not reserved.
#\'
#\' @noRd
covstruct_reserved_names <- function() {
  setdiff(unique(c(brms_whole_term_specials, "mo", "mi", "me", "gr", "mm",
                   "mmc", "offset", "I", "propto", "rr", "us_t", "diag_t",
                   "smooth", "s", "t2", "te", "ti", core_aterms)),
          names(extension_special_homes))
}', '#\' The names a registered structure may not take: the names the rest of
#\' frmtmb, brms and the other extensions own. A structure is a formula
#\' special, and the parser aliases every call to a special\'s name in a
#\' formula, so taking one of these would give one spelling two meanings
#\' in one formula: a structure named `poly` stopped `poly(x, 2)` from
#\' working anywhere.
#\'
#\' Reserved: core\'s formula vocabulary and internal structure names, the
#\' whole brms vocabulary, core\'s families, the specials, addition terms
#\' and families core lists as other extensions\' (only their home package
#\' registers those, `covstruct_home_owner()`), and the addition terms
#\' registered in this session. `covstruct_formula_function()` adds the
#\' functions of the packages a formula calls into.
#\'
#\' @noRd
covstruct_reserved_names <- function() {
  unique(c(brms_whole_term_specials, "mo", "mi", "me", "gr", "mm",
           "mmc", "offset", "I", "propto", "rr", "us_t", "diag_t",
           "smooth", "s", "t2", "te", "ti", core_aterms,
           unlist(brms_vocabulary, use.names = FALSE),
           names(family_registry), names(extension_family_homes),
           names(extension_aterm_homes), names(extension_special_homes),
           names(frmtmb_aterm_registry$reg)))
}

#\' The package whose function `name` is, when it is a function of a
#\' package formulas call into (base, stats, splines, utils, methods,
#\' mgcv), or `NULL`.
#\'
#\' @noRd
covstruct_formula_function <- function(name) {
  for (ns in c("base", "stats", "splines", "utils", "methods", "mgcv")) {
    env <- tryCatch(asNamespace(ns), error = function(e) NULL)
    if (!is.null(env) && exists(name, envir = env, inherits = FALSE,
                                mode = "function")) {
      return(ns)
    }
  }
  NULL
}

#\' The package that owns a reserved name core lists as an extension\'s,
#\' or `NULL` when no extension owns it.
#\'
#\' @noRd
covstruct_home_owner <- function(name) {
  homes <- c(extension_special_homes, extension_aterm_homes,
             extension_family_homes)
  if (name %in% names(homes)) unname(homes[[name]]) else NULL
}')
rep_file(p, '#\' `solve()`. A test that registers calls the suite\'s
#\' `local_frmtmb_registries()` first, which restores every registry
#\' when the test ends.', '#\' `solve()`. A test that registers calls [frmtmb_local_registries()]
#\' first, which restores every registry when the test ends.')
rep_file(p, '#\' @param name The structure\'s name, as a formula writes it, without
#\'   parentheses. It must be a syntactic name that is not one of
#\'   frmtmb\'s own structures or formula specials. Registering a name the
#\'   same package registered before replaces it, so a reload is not an
#\'   error.', '#\' @param name The structure\'s name, as a formula writes it, without
#\'   parentheses. It must be a syntactic name that nothing else in a
#\'   formula means: not one of frmtmb\'s own structures, specials,
#\'   families or addition terms, not a brms name, not a name core lists
#\'   as another extension\'s (only that package registers it), and not a
#\'   function of base, stats, splines, utils, methods or mgcv.
#\' @param replace A registered name is replaced when the package that
#\'   registered it registers it again, so a reload is not an error. A
#\'   name another package, or the session, registered is replaced only
#\'   with `replace = TRUE`; without it the call refuses.
#\' @param package The package that owns the structure. A fit records it
#\'   for each block, so that a fit read where that package is not
#\'   installed is refused by the package\'s name before any method uses
#\'   the structure. The default is the package whose namespace calls,
#\'   and `NULL` outside a package.')
rep_file(p, '                                      expand = NULL, term = NULL,
                                      from_natural = NULL, chol = NULL) {
  arg <- function(nm) paste0("frmtmb_register_covstruct(", nm, " =)")', '                                      expand = NULL, term = NULL,
                                      from_natural = NULL, chol = NULL,
                                      replace = FALSE,
                                      package = compat_registrant_name(
                                        parent.frame())) {
  arg <- function(nm) paste0("frmtmb_register_covstruct(", nm, " =)")
  if (identical(package, "frmtmb")) package <- NULL
  if (!is.null(package) && (!is.character(package) ||
                              length(package) != 1L || is.na(package) ||
                              !nzchar(package))) {
    frm_stop(arg("package"), " must be one package name, or NULL",
             call. = FALSE)
  }
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    frm_stop(arg("replace"), " must be TRUE or FALSE", call. = FALSE)
  }')
rep_file(p, '  if (name %in% covstruct_reserved_names()) {
    frm_stop("frmtmb_register_covstruct() cannot register \'", name,
             "\', because ", name, "() already means something else in a ",
             "formula. Choose another name.", call. = FALSE)
  }', '  home <- covstruct_home_owner(name)
  if (!is.null(home) && !identical(package, home)) {
    frm_stop("frmtmb_register_covstruct() cannot register \'", name,
             "\': it names a term of the extension package ", home,
             ", and only ", home, " registers it", call. = FALSE)
  }
  if (is.null(home) && name %in% covstruct_reserved_names()) {
    frm_stop("frmtmb_register_covstruct() cannot register \'", name,
             "\', because ", name, "() already means something else in a ",
             "formula. Choose another name.", call. = FALSE)
  }
  fpkg <- covstruct_formula_function(name)
  if (!is.null(fpkg)) {
    frm_stop("frmtmb_register_covstruct() cannot register \'", name,
             "\', because ", name, "() is a function of ", fpkg,
             " that formulas call, and a structure of that name would ",
             "take it from every formula. Choose another name.",
             call. = FALSE)
  }
  old <- frmtmb_covstruct_registry$reg[[name]]
  if (!is.null(old) && !replace && !identical(old[["package"]], package)) {
    frm_stop("frmtmb_register_covstruct(): \'", name, "\' is already ",
             "registered by ",
             if (is.null(old[["package"]])) "this session" else {
               paste0("the package ", old[["package"]])
             }, ". Only that package replaces it on its own; pass ",
             "replace = TRUE to replace it anyway", call. = FALSE)
  }')
rep_file(p, '                term = term, from_natural = from_natural, chol_L = chol,
                package = compat_registrant_name(parent.frame()))
  entry <- entry[!vapply(entry, is.null, TRUE)]', '                term = term, from_natural = from_natural, chol_L = chol,
                package = package)
  entry <- entry[!vapply(entry, is.null, TRUE)]
  bad <- cs_check_entry(name, entry)
  if (!is.null(bad)) {
    frm_stop("frmtmb_register_covstruct() would register an incomplete ",
             "entry: ", bad, call. = FALSE)
  }')
cat("fix-a applied to", root, "\n")
