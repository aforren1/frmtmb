# Write extensions/frmtmb.spatial/R/zzz.R with the compatibility notes
# copied byte for byte from core's R/compat.R at the base commit, so the
# rows frm_compat() reports with the extension loaded are the rows core
# reported (rule 4 of a move). Run once, from the worktree root.
x <- system2("git", c("show", "19b162b9:R/compat.R"), stdout = TRUE)
note_after <- function(pat) {
  i <- grep(pat, x, fixed = TRUE)
  stopifnot(length(i) == 1L)
  sub("^\\s*", "", x[i + 1L])
}
n_car <- note_after('  r("car", "*", "conditional",')
n_car <- sub(",$", "", n_car)
n_spde <- note_after('  r("spde", "*", "conditional",')
n_spde <- sub(",$", "", n_spde)
i <- grep('for (cs in frmtmb_compat_groups_lst$gmrf) {', x, fixed = TRUE)
n_sim <- sub("^\\s*", "", x[i + 2L])
n_pred <- sub("^\\s*", "", x[i + 4L])
rows <- function(nm, own) {
  c(sprintf('  r("%s", "*", "conditional",', nm),
    paste0("    ", own, ","), "    override = TRUE)",
    sprintf('  r("%s", "simulate", "works",', nm), paste0("    ", n_sim),
    sprintf('  r("%s", "predict", "conditional",', nm), paste0("    ", n_pred))
}
args <- commandArgs(trailingOnly = TRUE)
which <- if (length(args)) args else c("car", "spde")
own <- c(car = n_car, spde = n_spde)
out <- c(
  "# What this package tells frmtmb at load time: the structures it",
  "# supplies, through frmtmb_register_covstruct(), and their rows in",
  "# frm_compat(), through frmtmb_register_compat(). The structures are",
  "# registered FIRST, because registering one adds it to the",
  "# compatibility vocabulary that the rules then name.",
  "#",
  "# Every note below is core's note for the same pair before the move,",
  "# copied from frmtmb's R/compat.R (dev/spatial-gen-zzz.R in the frmtmb",
  "# repository writes this file), so frm_compat() reports the same rows",
  "# with this package loaded as core reported with the structures in",
  "# core. The rows for importance = and quadrature are not here:",
  "# frmtmb_register_covstruct() writes them for every structure it",
  "# registers, with the notes core's groups gave these two.",
  "",
  "#' @noRd",
  ".onLoad <- function(libname, pkgname) {",
  paste0("  register_", which, "()"),
  "  frmtmb_register_compat(rules = spatial_compat_rules)",
  "  invisible()",
  "}",
  "",
  "#' The compatibility rules for the structures this package supplies.",
  "#'",
  "#' @noRd",
  "spatial_compat_rules <- function() {",
  "  b <- compat_rule_builder()",
  "  r <- b$r",
  unlist(lapply(which, function(nm) rows(nm, own[[nm]]))),
  "  b$rules()",
  "}")
writeLines(out, "extensions/frmtmb.spatial/R/zzz.R")
cat("wrote zzz.R for", which, "\n")
