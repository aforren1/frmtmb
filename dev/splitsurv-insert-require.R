ins <- function(f, anchor, line, after = TRUE) {
  x <- readLines(f)
  i <- which(x == anchor)
  if (length(i) != 1L) stop(f, ": not once: ", anchor)
  x <- if (after) append(x, line, after = i) else append(x, line, after = i - 1L)
  writeLines(x, f)
}
ins("R/frame.R", "                           refit_na_action = NULL) {",
    c("  # a saved fit's family from a package this session lacks",
      "  require_model_packages(spec)"))
ins("R/fit.R", "                            \"vcov(), confint(), frm_linpred(se.fit =))\"))",
    "  require_model_packages(fit)")
ins("R/objective.R", "build_objective <- function(frame) {",
    "  require_model_packages(frame)")
ins("R/families.R", "response_mean <- function(fam, dpars, aterms) {",
    "  require_model_packages(fam)")
ins("R/families.R", "sim_draw <- function(ctx) {",
    "  require_model_packages(ctx[[\"family\"]])")
ins("R/curve.R", "  if (!inherits(object, \"frmtmb_fit\")) return(invisible(NULL))",
    "  require_model_packages(object)")
x <- readLines("R/families.R")
i <- which(x %in% c("  fam[[\"component_families\"]] <- vapply(comps, `[[`, \"\", \"family\")",
                    "  fam[[\"component_families\"]] <- vapply(cf, `[[`, \"\", \"family\")"))
stopifnot(length(i) == 2L)
v <- c("comps", "cf")
for (k in rev(seq_along(i))) {
  x <- append(x, c(
    "  # the packages whose functions the components carry, which a saved",
    "  # fit needs in the session that reads it (require_model_packages())",
    paste0("  attr(fam, \"frmtmb_component_packages\") <-"),
    paste0("    unique(vapply(", if (k == 1) "comps" else "cf",
           ", frm_family_package, \"\"))")), after = i[k])
}
writeLines(x, "R/families.R")
