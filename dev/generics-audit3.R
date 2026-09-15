# The names frmtmb masks that are NOT S3 generics.
#
# Found by accident: dev/generics-scale-brms.R called
# `brm(family = lognormal())` after attaching both packages and brms
# refused the object, because frmtmb's `lognormal()` had masked brms's.
# Items 2.5a to 2.5c are about GENERICS, where a shared generic fixes
# the collision.  A family constructor cannot be shared that way: the
# two return genuinely different objects.  This script measures how
# many such names there are, so the residue is on the record with a
# number rather than an impression.
LIB <- "C:/Users/adf44/source/r/generics-lib"
UL  <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
SYS <- file.path(R.home(), "library")
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib", UL, SYS))

frm <- parseNamespaceFile("frmtmb", dirname(find.package("frmtmb")))$exports
bexp <- parseNamespaceFile("brms", dirname(find.package("brms")))$exports
both <- sort(intersect(frm, bexp))

ns <- asNamespace("frmtmb")
is_gen <- function(nm) {
  f <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
  if (!is.function(f)) return(FALSE)
  any(grepl("UseMethod", deparse(body(f)), fixed = TRUE))
}
gen <- vapply(both, is_gen, NA)
# a name frmtmb re-exports from the owner is not a rival at all
shared <- vapply(both, function(nm) {
  f <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
  if (!is.function(f)) return(FALSE)
  e <- environmentName(topenv(environment(f)))
  nzchar(e) && e != "frmtmb"
}, NA)

cat("```\n")
cat("== names frmtmb and brms both export, dev/generics-audit3.R ==\n")
cat(sprintf("brms exports                       %d\n", length(bexp)))
cat(sprintf("frmtmb exports                     %d\n", length(frm)))
cat(sprintf("both                               %d\n", length(both)))
cat(sprintf("  of which S3 generics             %d\n", sum(gen)))
cat(sprintf("  of which already handed to owner %d\n", sum(shared)))
cat(sprintf("  of which NOT generics            %d\n", sum(!gen)))
cat("\nthe non-generic collisions:\n")
nm <- both[!gen]
for (i in seq(1, length(nm), by = 5)) {
  cat("  ", paste(nm[i:min(i + 4, length(nm))], collapse = ", "), "\n")
}
cat("```\n")
