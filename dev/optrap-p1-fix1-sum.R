# Lane optrap, punch round 1: dev/optrap-rev-fix1.R's starts per build
# (default plus 40 jittered, sd 0.3; w1: sd 1): how many hold the ledger
# row :326, how many end at code 0, and the lowest objective.
for (t in c("r7", "lane", "p1i", "r7w1", "lanew1", "p1iw1")) {
  f <- sprintf("dev/optrap-rev-log/fix1-%s.tsv", t)
  if (!file.exists(f)) next
  x <- read.delim(f)
  cat(sprintf("%-7s n %d | holds %d | code 0 %d | best %.4f | default %.10f code %d holds %s\n",
              t, nrow(x), sum(x$h326), sum(x$code == 0), min(x$obj),
              x$obj[1], x$code[1], x$h326[1]))
}
