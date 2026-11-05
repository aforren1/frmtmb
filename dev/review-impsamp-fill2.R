# Correct the Jacobian wording now that all six named tests have run.
p <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/reviews/2026-11-02-impsamp.md"
s <- paste(readLines(p), collapse = "\n")
sub1 <- function(a, b) {
  if (!grepl(a, s, fixed = TRUE)) stop("no: ", a)
  s <<- sub(a, b, s, fixed = TRUE)
}
sub1("reverting the total-derivative rate fails no test; half the Jacobian passes all but one",
     "reverting the total-derivative rate fails no test; no test compares the objective off the anchor with an exact value")
sub1("textbook, curvature, numDeriv, anchor tests all pass; only the gaussian test errors",
     "textbook, curvature, numDeriv, anchor tests pass; gaussian and scalar-intercept tests error")
sub1("gaussian errors, textbook fails 6 |", "gaussian and scalar-intercept tests error, textbook fails 6 |")
sub1("`half_jacobian` is a wrong density at every theta
but the anchor, and the only test that sees it is an exact-identity one:
nothing compares the objective off the anchor with an exact value.",
     "`half_jacobian` is a wrong density at every theta
but the anchor. Two of six named tests see it, both by an error at the end of the
fit and neither by comparing a value: the textbook test passes 15 of 15. Nothing
compares the objective off the anchor with an exact value.")
writeLines(s, p)
cat("ok\n")
