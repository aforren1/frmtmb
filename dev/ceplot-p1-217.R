# Lane ceplot: record the user's decision on row :217 (2026-09-30) in
# the findings and NEWS.
rep1 <- function(f, old, new) {
  s <- paste(readLines(f), collapse = "\n")
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop(basename(f), ": expected one match, got ", n)
  writeLines(strsplit(sub(old, new, s, fixed = TRUE), "\n", fixed = TRUE)[[1]], f)
}
fd <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-findings.md"
rep1(fd, "per effect; frmtmb once per call. Row :217 asserts the warning on the
DEFAULT call, so it does not hold; flipping it means changing the
default display, which is a decision for the user (section 8).",
"per effect; frmtmb once per call. Row :217 asserts the warning on the
DEFAULT call, so it does not hold. The user decided on 2026-09-30 to
keep frmtmb's per-category default, the display brms's own warning
recommends, as a deliberate divergence from brms; the row is
re-verdicted at consolidation.")
rep1(fd, "- The ordinal default stays the per-category display, and brms's
  \"treated as continuous\" warning comes with `categorical = FALSE`
  (1.5). brms's own warning calls its default likely invalid, which is
  the project's tiebreaker case; the user can overrule it.",
"- The ordinal default stays the per-category display, and brms's
  \"treated as continuous\" warning comes with `categorical = FALSE`
  (1.5). Reason: brms's own warning calls its default display likely
  invalid for ordinal families and asks for `categorical = TRUE`, which
  is frmtmb's default. The user decided this on 2026-09-30 (row :217,
  re-verdicted at consolidation).")
rep1(fd, "- Row :217 needs the default display changed; not done (section 8).",
"- Row :217 does not hold, by the user's decision of 2026-09-30 (a
  deliberate divergence, section 6).")
rep1(fd, "- Row :217: should `conditional_effects()` on an ordinal fit default to
  brms's expected-category display (with brms's warning), as brms does,
  instead of the per-category display that brms's warning recommends?
  This lane kept the per-category default.",
"- Nothing. Row :217 was decided by the user on 2026-09-30: the
  per-category default stays, as a deliberate divergence from brms.")
fn <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/NEWS.md"
rep1(fn, "  default per-category display gives none. brms gives the warning on
  its default call, whose display is the expected category number;
  frmtmb's default display is the per-category one that the warning
  asks for.",
"  default per-category display gives none. This is a deliberate
  divergence from brms, decided by the user: brms's default display is
  the expected category number and warns on every default call, and
  its own warning calls that display likely invalid for ordinal
  families and asks for `categorical = TRUE`, which is frmtmb's
  default.")
cat("edited\n")
