# Reviewer of lane optrap, item 7: per seed, lane (and lane with skip1)
# against rellib-r7 on the one-level design.
r <- read.delim("dev/optrap-rev-log/onelevel-r7.tsv")
for (lab in c("lane", "laneskip1")) {
  l <- read.delim(sprintf("dev/optrap-rev-log/onelevel-%s.tsv", lab))
  m <- merge(r, l, by = "seed", suffixes = c(".r7", ".l"))
  m$dobj <- m$obj.l - m$obj.r7
  cat(lab, ": code != 0 on", sum(m$code.l != 0), "| objective lane - r7:",
      "max over code-0 seeds", signif(max(abs(m$dobj[m$code.l == 0])), 3),
      "| on code != 0 seeds", paste(signif(m$dobj[m$code.l != 0], 3),
                                    collapse = " "), "\n")
  x <- m[m$code.l != 0, ]
  print(data.frame(seed = x$seed, dobj = signif(x$dobj, 3),
                   db0 = signif(x$b0.l - x$b0.r7, 3),
                   dlogsig = signif(x$logsig.l - x$logsig.r7, 3),
                   logsd.r7 = signif(x$logsd.r7, 3),
                   logsd.l = signif(x$logsd.l, 3)))
}
