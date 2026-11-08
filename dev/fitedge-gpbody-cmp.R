# Item 1: row by row, base against lane (dev/fitedge-gpbody-sweep.R logs), and
# time by body. Worse means more than 1e-3 nats below base.
o <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-out"
rd <- function(f) read.delim(file.path(o, f), stringsAsFactors = FALSE)
k <- c("body", "curve", "n", "seed")
m <- merge(rd("gpbody-fit-base.tsv"), rd("gpbody-fit-lane.tsv"), by = k,
           suffixes = c(".b", ".l"))
m <- merge(m, setNames(rd("gpbody-ref-lane.tsv")[, c(k, "ll")], c(k, "ref")), by = k)
for (b in unique(m$body)) {
  s <- m[m$body == b, ]
  cat(sprintf("CMP body %s time base %.1f s lane %.1f s ratio %.2f | lane better %d worse %d\n",
              b, sum(s$sec.b), sum(s$sec.l), sum(s$sec.l) / sum(s$sec.b),
              sum(s$ll.l > s$ll.b + 1e-3), sum(s$ll.l < s$ll.b - 1e-3)))
}
w <- m[m$ll.l < m$ll.b - 1e-3, ]
for (i in seq_len(nrow(w))) {
  cat(sprintf("WORSE %s %s n %d seed %d base %.4f lane %.4f ref %.4f\n", w$body[i],
              w$curve[i], w$n[i], w$seed[i], w$ll.b[i], w$ll.l[i], w$ref[i]))
}
# the body against the linear spelling's own default, where one exists
li <- rd("gpbody-lin-lane.tsv")
ml <- merge(m, setNames(li[, c(k, "ll", "sec")], c(k, "lin", "lin_sec")), by = k)
for (b in unique(ml$body)) {
  s <- ml[ml$body == b, ]
  cat(sprintf("LIN body %s rows %d body equals linear (1e-6 rel) %d, body above %d, below %d; time body %.1f s linear %.1f s\n",
              b, nrow(s), sum(abs(s$ll.l - s$lin) <= 1e-6 * abs(s$lin)),
              sum(s$ll.l > s$lin + 1e-6 * abs(s$lin)),
              sum(s$ll.l < s$lin - 1e-6 * abs(s$lin)), sum(s$sec.l), sum(s$lin_sec)))
  for (i in which(abs(s$ll.l - s$lin) > 1e-6 * abs(s$lin))) {
    cat(sprintf("LINDIFF %s %s n %d seed %d body %.4f linear %.4f base body %.4f ref %.4f\n",
                b, s$curve[i], s$n[i], s$seed[i], s$ll.l[i], s$lin[i], s$ll.b[i], s$ref[i]))
  }
}
