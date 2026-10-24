# Fix round 1: quadrature and importance refuse a registered structure in
# frm_compat() as frm() does (A6), and contributed rule sets resolve in
# the order of their packages' names, not of library() calls (C2), on
# the tree named as the argument.
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
p <- "R/compat.R"
rep_file(p, '  ## quadrature --------------------------------------------------------
', '  ## quadrature --------------------------------------------------------
  # A structure another package registers is on neither of the groups
  # below, and frm() refuses it by name; this kind rule says so where
  # the groups, more specific, do not reach.
  r("quadrature", "kind:covstruct", "refused",
    "Refused: quadrature marginalizes one scalar random intercept at a time, and frm() refuses every block that is not a dimension-1 us, diag or homdiag term (or a Student-t latent over one), which a structure another package registers is not.")
')
rep_file(p, '  r("importance", "group:importance_blocks", "conditional",', '  # the kind rule, below the groups, is for a structure another package
  # registers: the scope check knows frmtmb\'s own structures only
  r("importance", "kind:covstruct", "refused",
    "Refused by name at the scope check: the correction recovers each level\'s prior density from the block density, which it does for frmtmb\'s own structures that are gaussian in a level\'s coefficients and independent between levels, and for no structure another package registers.")
  r("importance", "group:importance_blocks", "conditional",')
rep_file(p, '  arg <- "frmtmb_register_compat()"
  who <- compat_registrant(parent.frame())', '  arg <- "frmtmb_register_compat()"
  who <- compat_registrant(parent.frame())
  owner <- compat_registrant_name(parent.frame())')
rep_file(p, '  if (!is.null(rules)) {
    frmtmb_compat_contrib$rules <-
      c(frmtmb_compat_contrib$rules, list(rules))
  }
  invisible(NULL)
}', '  if (!is.null(rules)) {
    # the owner orders the contributions (frmtmb_compat_rules_tbl())
    attr(rules, "frmtmb_package") <- owner %||% ""
    frmtmb_compat_contrib$rules <-
      c(frmtmb_compat_contrib$rules, list(rules))
  }
  invisible(NULL)
}')
rep_file(p, 'frmtmb_compat_rules_tbl <- function() {
  hand <- compat_hand_rules_tbl()
  out <- do.call(rbind, c(
    list(hand, compat_aterm_rules(compat_core_family_accepts())),
    lapply(frmtmb_compat_contrib$rules, function(mk) mk())))', 'frmtmb_compat_rules_tbl <- function() {
  hand <- compat_hand_rules_tbl()
  out <- do.call(rbind, c(
    list(hand, compat_aterm_rules(compat_core_family_accepts())),
    lapply(compat_contrib_rules_ordered(), function(mk) mk())))')
rep_file(p, '# The hand-written rules do not change within a session:', '#\' The contributed rule sets in the order they resolve: by the name of
#\' the package that registered each (a set registered outside a package
#\' first), and within one package in registration order.
#\'
#\' Two contributed rules of one signature tie, and the later wins. In
#\' load order that made frm_compat() depend on the order of the
#\' library() calls: with frmtmb.spatial loaded before frmtmb.sample, two
#\' pairs showed the other package\'s note. Ordered by package, the
#\' answer is a property of the packages installed, not of the session.
#\'
#\' @noRd
compat_contrib_rules_ordered <- function() {
  rl <- frmtmb_compat_contrib$rules
  if (length(rl) < 2L) return(rl)
  pk <- vapply(rl, function(f) {
    attr(f, "frmtmb_package", exact = TRUE) %||% ""
  }, "")
  rl[order(pk, seq_along(pk))]
}

# The hand-written rules do not change within a session:')
rep_file(p, '#\' Contributions are APPENDED, which is the only safe direction. Rules
#\' of equal specificity resolve later-wins, so a contributed rule may
#\' override a core default and a core default can never silently
#\' override a contributed one.', '#\' Contributions are APPENDED, which is the only safe direction. Rules
#\' of equal specificity resolve later-wins, so a contributed rule may
#\' override a core default and a core default can never silently
#\' override a contributed one. Between packages, contributions resolve
#\' in the order of the packages\' names, so the table does not depend
#\' on the order in which they were loaded.')
cat("fix-a6 applied to", root, "\n")
