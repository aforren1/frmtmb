# Review finding 10: the old subsite's articles index and authors page
# get redirects, its NEWS page points at the page that holds its
# history, and the comment above the block says which lines are hand
# written. Run once.
f <- "_pkgdown.yml"
x <- readLines(f)
one <- function(s) {
  i <- which(x == s)
  if (length(i) != 1L) stop("not once: ", s)
  i
}
i <- one("  - [\"frmtmb.spline/news/index.html\", \"frmtmb.survival/news/index.html\"]")
x[i] <- "  - [\"frmtmb.spline/news/index.html\", \"frmtmb.survival/articles/spline-history.html\"]"
x <- append(x, c(
  "  - [\"frmtmb.spline/articles/index.html\", \"articles/index.html\"]",
  "  - [\"frmtmb.spline/authors.html\", \"frmtmb.survival/authors.html\"]"),
  after = i)
i <- one("# dev/docsci-gen-redirects.R, which follows each chain to a live page.")
x <- append(x, c(
  "# The first block below is not generated: dev/splitsurv-retire-pkgdown.R",
  "# and dev/splitsurv-edit-pkgdown.R wrote it when frmtmb.spline was",
  "# retired and cox_baseline() moved to frmtmb.survival. The pages of the",
  "# retired subsite go to the pages that now hold the same functions,",
  "# articles and history."), after = i)
writeLines(x, f)
