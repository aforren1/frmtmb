# Review of lane fitedge, item 1: summary of review-fitedge-gpfresh.R. Args:
# set (fresh | lane96). The reference is the best of ref and every arm.
# Short: more than 1e-3 nats below it. Worse: more than 1e-3 below base.
set <- commandArgs(TRUE)[1]
# a second argument names the lane arm's tsv (the fix round's rerun)
lane_file <- commandArgs(TRUE)[2]
o <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/review-fitedge-out"
lo <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-out"
rd <- function(f) utils::read.delim(f, colClasses = c(theta = "character",
                                                       wclass = "character"))
if (set == "fresh") {
  b <- rd(file.path(o, "gpf-fresh-base.tsv"))
  l <- rd(file.path(o, "gpf-fresh-lane.tsv"))
  p <- rd(file.path(o, "gpf-fresh-pick.tsv"))
  r <- rd(file.path(o, "gpf-fresh-ref.tsv"))
} else {
  rd2 <- function(f) {
    x <- utils::read.delim(f, colClasses = c(wclass = "character"))
    names(x)[1] <- "design"
    x
  }
  b <- rd2(file.path(lo, "gpbody-fit-base.tsv"))
  l <- rd2(file.path(lo, "gpbody-fit-lane.tsv"))
  r <- rd2(file.path(lo, "gpbody-ref-lane.tsv"))
  p <- rd(file.path(o, "gpf-lane96-pick.tsv"))
}
if (!is.na(lane_file)) l <- rd(lane_file)
k <- c("design", "curve", "n", "seed")
stopifnot(identical(b[k], l[k]), identical(b[k], r[k]), identical(b[k], p[k]))
ref <- pmax(r$ll, b$ll, l$ll, p$ll, na.rm = TRUE)
two <- pmax(b$ll, l$ll)
tol <- 1e-3
silent <- function(x) is.na(x$wclass) | x$wclass == ""
line <- function(lab, i) {
  sh <- function(v) sum(v[i] < ref[i] - tol)
  wo <- function(v) min(c(0, (v - ref)[i]))
  cat(sprintf(paste0("GS %s %-6s rows %3d | short base %2d (silent %2d) lane %2d (silent %2d) ",
                     "pick %2d best-of-two %2d | worst base %8.3f lane %8.3f two %8.3f | ",
                     "lane worse than base %d (max %.3f) better %d | time base %.1f lane %.1f\n"),
              set, lab, sum(i), sh(b$ll), sum(i & b$ll < ref - tol & silent(b)),
              sh(l$ll), sum(i & l$ll < ref - tol & silent(l)), sh(p$ll), sh(two),
              wo(b$ll), wo(l$ll), wo(two), sum(i & l$ll < b$ll - tol),
              max(c(0, (b$ll - l$ll)[i & l$ll < b$ll - tol])),
              sum(i & l$ll > b$ll + tol), sum(b$sec[i]), sum(l$sec[i])))
}
line("all", rep(TRUE, nrow(b)))
for (d in unique(b$design)) line(d, b$design == d)
w <- which(l$ll < b$ll - tol)
for (i in w) {
  cat(sprintf("GSW %s %s %s %d %d base %.4f lane %.4f ref %.4f lane-base %.4f lane warn '%s'\n",
              set, b$design[i], b$curve[i], b$n[i], b$seed[i], b$ll[i], l$ll[i],
              ref[i], l$ll[i] - b$ll[i], l$wclass[i]))
}
s <- which(l$ll < ref - tol)
for (i in s) {
  cat(sprintf("GSS %s %s %s %d %d lane %.4f ref %.4f short %.4f base %.4f warn '%s'\n",
              set, b$design[i], b$curve[i], b$n[i], b$seed[i], l$ll[i], ref[i],
              l$ll[i] - ref[i], b$ll[i], l$wclass[i]))
}
