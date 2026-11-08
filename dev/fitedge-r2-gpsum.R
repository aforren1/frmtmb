# Fix round, B1: the fix build (lane2, dev/fitedge-r2-gpfresh.R) against
# rellib-r18 and the reviewed build on the review's two 96-row sets. The
# reference is the best of the reference arm and every arm, as in
# dev/review-fitedge-gpsum.R. Short: more than 1e-3 nats below it; worse:
# more than 1e-3 below rellib-r18.
o <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/review-fitedge-out"
lo <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-out"
r2 <- file.path(lo, "r2")
rd <- function(f) {
  x <- utils::read.delim(f, colClasses = c(wclass = "character"))
  names(x)[1] <- "design"
  x
}
tol <- 1e-3
for (set in c("fresh", "lane96")) {
  if (set == "fresh") {
    b <- rd(file.path(o, "gpf-fresh-base.tsv"))
    l <- rd(file.path(o, "gpf-fresh-lane.tsv"))
    r <- rd(file.path(o, "gpf-fresh-ref.tsv"))
  } else {
    b <- rd(file.path(lo, "gpbody-fit-base.tsv"))
    l <- rd(file.path(lo, "gpbody-fit-lane.tsv"))
    r <- rd(file.path(lo, "gpbody-ref-lane.tsv"))
  }
  n2 <- rd(file.path(r2, sprintf("gpf-%s-lane2.tsv", set)))
  k <- c("design", "curve", "n", "seed")
  stopifnot(identical(b[k], n2[k]), identical(b[k], l[k]), identical(b[k], r[k]))
  ref <- pmax(r$ll, b$ll, l$ll, n2$ll, na.rm = TRUE)
  silent <- function(x) is.na(x$wclass) | x$wclass == ""
  line <- function(lab, i) {
    sh <- function(x) sum(i & x$ll < ref - tol)
    si <- function(x) sum(i & x$ll < ref - tol & silent(x))
    wo <- function(x) min(c(0, (x$ll - ref)[i]))
    cat(sprintf(paste0("R2 %s %-4s rows %2d | short base %2d (silent %2d) reviewed %2d ",
                       "fix %2d (silent %2d) | worst base %8.3f fix %8.3f | fix worse ",
                       "than base %d (max %.3f), better %d | time base %.1f fix %.1f ",
                       "ratio %.2f\n"),
                set, lab, sum(i), sh(b), si(b), sh(l), sh(n2), si(n2), wo(b), wo(n2),
                sum(i & n2$ll < b$ll - tol),
                max(c(0, (b$ll - n2$ll)[i & n2$ll < b$ll - tol])),
                sum(i & n2$ll > b$ll + tol), sum(b$sec[i]), sum(n2$sec[i]),
                sum(n2$sec[i]) / sum(b$sec[i])))
  }
  line("all", rep(TRUE, nrow(b)))
  for (d in unique(b$design)) line(d, b$design == d)
  for (i in which(n2$ll < b$ll - tol)) {
    cat(sprintf("R2W %s %s %s %d %d base %.4f fix %.4f\n", set, b$design[i],
                b$curve[i], b$n[i], b$seed[i], b$ll[i], n2$ll[i]))
  }
  for (i in which(n2$ll < ref - tol)) {
    cat(sprintf("R2S %s %s %s %d %d fix %.4f ref %.4f short %.4f base %.4f warn '%s'\n",
                set, b$design[i], b$curve[i], b$n[i], b$seed[i], n2$ll[i], ref[i],
                n2$ll[i] - ref[i], b$ll[i], n2$wclass[i]))
  }
}
