# Summary of dev/predraw2-ncheck.R: the check's chance of firing, by how
# many test draws it reads, on fits whose 1000-draw median excess is
# below 0.5 (it should not fire) and above 1.5 (it should).
# Usage: Rscript dev/predraw2-ncheck-sum.R <tsv>...
for (f in commandArgs(TRUE)) {
  x <- utils::read.delim(f)
  lo <- x$med1000 < 0.5
  hi <- x$med1000 > 1.5
  cat(sprintf(paste("%-34s fits %2d | median < 0.5: n %2d, mean P(fire)",
                    "n20 %.4f n40 %.4f seq %.4f, max n20 %.4f seq %.4f |",
                    "median > 1.5: n %d, min P(fire) n20 %s seq %s\n"),
              basename(f), nrow(x), sum(lo), mean(x$p20[lo]),
              mean(x$p40[lo]), mean(x$pseq[lo]), max(x$p20[lo]),
              max(x$pseq[lo]), sum(hi),
              if (any(hi)) sprintf("%.4f", min(x$p20[hi])) else "-",
              if (any(hi)) sprintf("%.4f", min(x$pseq[hi])) else "-"))
}
