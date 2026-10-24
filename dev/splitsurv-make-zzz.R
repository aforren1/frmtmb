# Build frmtmb.survival's R/zzz.R from frmtmb.spline's (the Royston and
# Parmar rows) and the cox() rows cut from core's R/compat.R. Run once.
sp <- readLines("extensions/frmtmb.spline/R/zzz.R")
a <- grep("^  ## ---- the family ---", sp)
b <- grep("^  ## ---- curve inference ---", sp)
c1 <- grep("^  r[(]\"rp_floored\", \"cens[(][)]\", \"works\",", sp)
e <- grep("^  b[$]rules[(][)]$", sp)
stopifnot(length(a) == 1L, length(b) == 1L, length(c1) == 1L,
          length(e) == 1L)
rp_rules <- sp[(a + 1L):(b - 1L)]
rp_floored_rules <- sp[c1:(e - 1L)]
cox <- readLines("dev/splitsurv-out/cox-compat-cut.R")
head <- c(
  "# What this package tells frmtmb about itself at load time.",
  "#",
  "# The registration uses the seam frmtmb exports for the purpose",
  "# (?frmtmb::`frmtmb-sampling-api`, section \"The compatibility registry\").",
  "# Registering from .onLoad() rather than at top level is what a",
  "# contributor outside the package must do: by then every namespace is",
  "# sealed, so the collation-order question that governs frmtmb's own",
  "# in-package contributors does not arise.",
  "#",
  "# The rows are the ones core and frmtmb.spline registered for these",
  "# features before the two families moved here, word for word, so that",
  "# frm_compat() reports the same pairs with this package loaded as it",
  "# did before (dev/splitsurv-compat-dump.R compares the two sessions).",
  "",
  "#' @noRd",
  ".onLoad <- function(libname, pkgname) {",
  "  frmtmb_register_compat(",
  "    features = c(cox = \"family\", royston_parmar = \"family\",",
  "                 rp_floored = \"method\"),",
  "    rules = sv_compat_rules)",
  "  invisible()",
  "}",
  "",
  "#' The compatibility rules for what this package supplies.",
  "#'",
  "#' Every row was run. `untested` is used where the pair was not",
  "#' exercised, which is the honest third state: an absent guard and a",
  "#' passing guard look the same from outside.",
  "#'",
  "#' @noRd",
  "sv_compat_rules <- function() {",
  "  b <- compat_rule_builder()",
  "  r <- b$r",
  "",
  "  ## ---- cox() ------------------------------------------------------")
out <- c(head, cox, "",
         "  ## ---- royston_parmar() --------------------------------------",
         rp_rules, rp_floored_rules, "  b$rules()", "}")
writeLines(out, "extensions/frmtmb.survival/R/zzz.R")
cat(length(out), "lines\n")
