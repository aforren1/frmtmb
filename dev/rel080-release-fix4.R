# The short check of the release fixes (dev/reviews/2026-11-07-release.md,
# "Short check of the release fixes", SHIP): its backlog edit and its NIT,
# after dev/rel080-release-fix3.R. NEWS is edited by dev/rel080-news.R.
#
#   Rscript dev/rel080-release-fix4.R
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel080"
rd <- function(p) paste(readLines(file.path(REL, p), warn = FALSE),
                        collapse = "\n")
wr <- function(x, p) {
  con <- file(file.path(REL, p), open = "wb")
  on.exit(close(con))
  writeLines(strsplit(x, "\n", fixed = TRUE)[[1L]], con, sep = "\r\n")
}
edit <- function(p, old, new, what) {
  x <- rd(p)
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(p, ": ", what, ": the old text occurs ", n, " times")
  wr(sub(old, new, x, fixed = TRUE), p)
  cat("EDIT", p, "|", what, "\n")
}
edit("dev/test-backlog.md",
paste0("  standard-error warning (release review of 0.80.0, item 1). 0.80.0 keeps\n"),
paste0("  standard-error warning (release review of 0.80.0, item 1). The refit of\n",
"  a saddle restart does the same: bernoulli `sig4` seed 55 of the review's\n",
"  fresh2 rows reads -189.696351, 2.8e-4 below the best start, silent, where\n",
"  0.79.1 warned of its gradient (release review, short check). 0.80.0 keeps\n"),
"short check: the ridge item names the saddle-restart refit")
edit("dev/test-backlog.md",
"### From the reviews, filed at the release\n\n",
paste0("### From the reviews, filed at the release\n\n",
"- **NIT. `design_center()` centers `0 + w1 * w2` at offsets of 1e5 and\n",
"  1e6**, a design that does not span the constant, because its absolute\n",
"  tolerance (1e-8 per row on the raw columns) counts the constant as in\n",
"  that span. The rank decision and the fit are right on all 6 fits\n",
"  (seeds 1 to 3, against `lm()`; release review of 0.80.0, short check).\n"),
"short check: the NIT on 0 + w1 * w2")
cat("RELEASE FIX 4 DONE\n")
