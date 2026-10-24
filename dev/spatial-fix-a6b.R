# Fix round 1, A6 (second form): the kind rules fix-a6 wrote win no pair
# (every core structure sits in a group, which is more specific), so
# the rows move to frmtmb_register_covstruct(), which writes them for
# each structure it registers. The notes become constants that the
# group rules share. Also makes the replace = refusal a unique template
# and renames a local in R/registries.R that the hazard scanner bans.
# Run on the tree named as the argument.
root <- commandArgs(trailingOnly = TRUE)[1]
rd <- function(f) readLines(file.path(root, f), warn = FALSE)
wr <- function(x, f) writeLines(x, file.path(root, f))
one <- function(x, pat) {
  i <- grep(pat, x, fixed = TRUE)
  if (length(i) != 1L) stop("expected one match for: ", pat)
  i
}

## R/compat.R ---------------------------------------------------------
x <- rd("R/compat.R")
# drop the two kind rules and the comments that introduce them
i <- one(x, 'r("importance", "kind:covstruct", "refused",')
stopifnot(grepl("^  # ", x[i - 1L]), grepl("^  # ", x[i - 2L]))
n_scope <- sub('^\\s*"(.*)"\\)$', "\\1", x[i + 1L])
x <- x[-c(i - 2L, i - 1L, i, i + 1L)]
i <- one(x, 'r("quadrature", "kind:covstruct", "refused",')
stopifnot(all(grepl("^  # ", x[i - (1:3)])))
x <- x[-c(i - 3L, i - 2L, i - 1L, i, i + 1L)]
# the group notes become constants the seam shares
i <- one(x, 'r("importance", "group:crosslevel_blocks", "refused",')
n_cross <- sub('^\\s*"(.*)"\\)$', "\\1", x[i + 1L])
x[i + 1L] <- "    compat_note_crosslevel_importance)"
i <- one(x, 'r("quadrature", "group:wide_blocks", "refused",')
n_wide <- sub('^\\s*"(.*)"\\)$', "\\1", x[i + 1L])
x[i + 1L] <- "    compat_note_wide_quadrature)"
j <- one(x, "compat_hand_rules_tbl <- function() {")
# back up to the start of the roxygen block above the function
k <- j - 1L
while (k > 1L && grepl("^#'", x[k - 1L])) k <- k - 1L
consts <- c(
  "# The notes of the rows that refuse a covariance structure under",
  "# importance = and quadrature. The groups below give them to frmtmb's",
  "# own structures and frmtmb_register_covstruct() to each structure it",
  "# registers, so a structure moved to an extension reports what it",
  "# reported in core.",
  paste0('compat_note_crosslevel_importance <- "', n_cross, '"'),
  paste0('compat_note_scope_importance <- "', n_scope, '"'),
  paste0('compat_note_wide_quadrature <- "', n_wide, '"'),
  "")
x <- append(x, consts, after = k - 1L)
wr(x, "R/compat.R")

## R/covstruct-seam.R -------------------------------------------------
x <- rd("R/covstruct-seam.R")
i <- one(x, "  frmtmb_covstruct_registry$reg[[name]] <- entry")
add <- c(
  "  # Quadrature and importance = refuse every registered structure, and",
  "  # frm_compat() has to say so. The rows are keyed by the structure, so",
  "  # a replacement replaces them rather than leaving a row that wins",
  "  # nothing.",
  "  rules <- covstruct_compat_rules(name, levels_coupled)",
  "  attr(rules, \"frmtmb_package\") <- package %||% \"\"",
  "  frmtmb_compat_contrib$rules[[paste0(\"covstruct:\", name)]] <- rules")
x <- append(x, add, after = i - 1L)
# the package's shared TRUE/FALSE refusal, not a second copy of it
i <- one(x, '    frm_stop(arg("replace"), " must be TRUE or FALSE", call. = FALSE)')
stopifnot(grepl("^  if [(]!is.logical[(]replace[)]", x[i - 1L]),
          x[i + 1L] == "  }")
x <- c(x[seq_len(i - 2L)], '  check_flag(replace, arg("replace"))',
       x[-seq_len(i + 1L)])
x <- c(x, "",
  "#' The frm_compat() rows frmtmb writes for a registered structure.",
  "#'",
  "#' @noRd",
  "covstruct_compat_rules <- function(name, levels_coupled) {",
  "  force(name)",
  "  force(levels_coupled)",
  "  function() {",
  "    b <- compat_rule_builder()",
  "    b$r(\"importance\", name, \"refused\",",
  "        if (isTRUE(levels_coupled)) compat_note_crosslevel_importance",
  "        else compat_note_scope_importance)",
  "    b$r(\"quadrature\", name, \"refused\", compat_note_wide_quadrature)",
  "    b$rules()",
  "  }",
  "}")
wr(x, "R/covstruct-seam.R")

## R/registries.R -----------------------------------------------------
x <- rd("R/registries.R")
n0 <- sum(grepl("\\baterms\\b", x))
x <- gsub("\\baterms\\b", "atreg", x)
cat("registries.R: renamed", n0, "lines\n")
wr(x, "R/registries.R")
cat("fix-a6b applied to", root, "\n")
