# Lane gpstop: summary of dev/gpstop-gpfresh.R on the two 96-row sets. The
# reference per row is the best of the ref arm and every arm. Short: more
# than 1e-3 nats below it; silent: no warning. Worse or better: more than
# 1e-3 below or above rellib-r20 (base). Args: lane arm name (default lane).
lane_arm <- commandArgs(TRUE)[1]
if (is.na(lane_arm)) lane_arm <- "lane"
o <- "C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-out"
rd <- function(f) {
  x <- utils::read.delim(file.path(o, f),
                         colClasses = c(wclass = "character",
                                        theta = "character"))
  names(x)[1] <- "design"
  x
}
tol <- 1e-3
for (set in c("fresh", "lane96")) {
  b <- rd(sprintf("gpf-%s-base.tsv", set))
  l <- rd(sprintf("gpf-%s-%s.tsv", set, lane_arm))
  r <- rd(sprintf("gpf-%s-ref.tsv", set))
  k <- c("design", "curve", "n", "seed")
  stopifnot(identical(b[k], l[k]), identical(b[k], r[k]))
  ref <- pmax(r$ll, b$ll, l$ll, na.rm = TRUE)
  silent <- function(x) is.na(x$wclass) | x$wclass == ""
  line <- function(lab, i) {
    sh <- function(x) sum(i & x$ll < ref - tol)
    si <- function(x) sum(i & x$ll < ref - tol & silent(x))
    wo <- function(x) min(c(0, (x$ll - ref)[i]))
    cat(sprintf(paste0("GPS %s %-4s rows %2d | short base %2d (silent %2d) ",
                       "lane %2d (silent %2d) | worst base %8.3f lane %8.3f | ",
                       "lane worse than base %d (max %.3f), better %d | ",
                       "time base %.1f lane %.1f\n"),
                set, lab, sum(i), sh(b), si(b), sh(l), si(l), wo(b), wo(l),
                sum(i & l$ll < b$ll - tol),
                max(c(0, (b$ll - l$ll)[i & l$ll < b$ll - tol])),
                sum(i & l$ll > b$ll + tol), sum(b$sec[i]), sum(l$sec[i])))
  }
  line("all", rep(TRUE, nrow(b)))
  for (d in unique(b$design)) line(d, b$design == d)
  for (i in which(l$ll < b$ll - tol)) {
    cat(sprintf("GPSW %s %s %s %d %d base %.4f lane %.4f\n", set, b$design[i],
                b$curve[i], b$n[i], b$seed[i], b$ll[i], l$ll[i]))
  }
  for (i in which(l$ll > b$ll + tol)) {
    cat(sprintf("GPSB %s %s %s %d %d base %.4f lane %.4f gain %.4f\n", set,
                b$design[i], b$curve[i], b$n[i], b$seed[i], b$ll[i], l$ll[i],
                l$ll[i] - b$ll[i]))
  }
  for (i in which(l$ll < ref - tol)) {
    cat(sprintf(paste0("GPSS %s %s %s %d %d lane %.4f ref %.4f short %.4f ",
                       "base %.4f warn '%s'\n"),
                set, b$design[i], b$curve[i], b$n[i], b$seed[i], l$ll[i],
                ref[i], l$ll[i] - ref[i], b$ll[i], l$wclass[i]))
  }
}
