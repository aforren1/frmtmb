# Fix round 1's "seen to fail on the build before the fix": switch the
# four fix-round guards off in the source (or back on), so the lane
# library can hold the pre-fix behaviour for one test run.
#   Rscript dev/famsimplex-prefix.R off|on
# "off" saves the four files under dev/famsimplex-prefix-bak/ first;
# "on" restores them byte for byte from there.
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex"
files <- c("R/parse.R", "R/diagnostics.R", "R/conditional-effects.R",
           "R/predict.R")
bak <- file.path(wt, "dev", "famsimplex-prefix-bak")
if (identical(arm, "on")) {
  for (f in files) {
    stopifnot(file.copy(file.path(bak, basename(f)), file.path(wt, f),
                        overwrite = TRUE))
  }
  cat("restored\n")
  quit(save = "no")
}
dir.create(bak, showWarnings = FALSE)
for (f in files) {
  stopifnot(file.copy(file.path(wt, f), file.path(bak, basename(f)),
                      overwrite = TRUE))
}
swap <- function(f, from, to) {
  p <- file.path(wt, f)
  s <- readLines(p)
  i <- which(s == from)
  stopifnot(length(i) >= 1L)
  s[i] <- to
  writeLines(s, p)
  cat(f, ":", length(i), "line(s)\n")
}
swap("R/parse.R", "  if (length(primaries) > 1L) {",
     "  if (FALSE && length(primaries) > 1L) {")
swap("R/diagnostics.R", "  if (fam_is_simplex(rspec$family)) {",
     "  if (FALSE) {")
swap("R/conditional-effects.R",
     "  if (identical(method, \"predict\") && fam_is_simplex(rspec$family)) {",
     "  if (FALSE) {")
# residuals_one()'s early refusal; residual_values() keeps its own, as
# the pre-fix build had it
p <- file.path(wt, "R/predict.R")
s <- readLines(p)
i <- which(s == "  if (fam_is_simplex(fam)) simplex_residuals_refusal(fam)")
stopifnot(length(i) == 2L)
s[i[1L]] <- "  if (FALSE) simplex_residuals_refusal(fam)"
writeLines(s, p)
cat("R/predict.R : 1 line\n")
