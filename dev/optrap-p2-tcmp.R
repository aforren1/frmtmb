# Lane optrap, punch round 2: the review's heavy-tail cases (student
# fits of t(1) and t(2) data at y x 1 and 1e3, dev/optrap-rev2-b1.R) per
# seed, this build (tag p2, dev/optrap-log/rev2/) against rellib-r7 and
# punch round 1's build 26 (the reviewer's logs, read only). Short means
# more than 1e-4 below the best any of them reached at either scale.
for (cs in c("t1", "t2")) {
  rd <- function(f) read.delim(f, quote = "")
  b <- list(p2 = rd(sprintf("dev/optrap-log/rev2/b1-%s-p2.tsv", cs)),
            b26 = rd(sprintf("dev/optrap-rev2-log/b1-%s-lane.tsv", cs)),
            r7 = rd(sprintf("dev/optrap-rev2-log/b1-%s-r7.tsv", cs)))
  best <- tapply(pmax(b$p2$ll, b$b26$ll, b$r7$ll), b$p2$seed, max)
  for (sc in c(1, 1000)) {
    out <- vapply(names(b), function(k) {
      x <- b[[k]]
      i <- x$setting == sc
      gp <- x$ll[i] - best[as.character(x$seed[i])]
      sh <- gp < -1e-4
      sprintf("%s %d (min %.3g, silent %d)", k, sum(sh), min(gp),
              sum(sh & x$nwarn[i] == 0))
    }, "")
    cat(cs, "scale", sc, "|", paste(out, collapse = " | "), "\n")
  }
}
