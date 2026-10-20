# Lane illse2, item 1: summarize dev/illse2-log/p1-<arm>-<re>-<deg>.txt
# into one line per design and route: the largest |ratio - 1| against
# the exact reference over seeds 1 to 3, or the error count.
#   Rscript dev/illse2-p1sum.R [suffix]
a <- commandArgs(TRUE)
sfx <- if (length(a)) a[1] else ""
dir <- "C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-log"
for (arm in c("base", "lane")) for (re in 0:1) for (deg in 5:8) {
  f <- file.path(dir, sprintf("p1-%s-%d-%d%s.txt", arm, re, deg, sfx))
  if (!file.exists(f)) next
  x <- readLines(f)
  one <- function(tag) {
    l <- grep(tag, x, value = TRUE)
    if (!length(l)) return("-")
    nerr <- sum(grepl("ERROR", l))
    ok <- l[!grepl("ERROR", l)]
    dv <- regmatches(ok, regexpr("max dev [0-9.e+-]+", ok))
    st <- regmatches(ok, regexpr("stat [0-9.]+", ok))
    v <- c(as.numeric(sub("max dev ", "", dv)),
           abs(as.numeric(sub("stat ", "", st)) - 1))
    paste0(if (length(v)) sprintf("%.1e", max(v)) else "",
           if (nerr) sprintf("%d of %d error", nerr, length(l)) else "")
  }
  cat(sprintf("%-4s %-6s deg %d | cooks %-12s | D1 all %-12s | D2 all %s\n",
              arm, if (re) "(1|g)" else "no RE", deg, one("^  cooks"),
              one("^  D1 all"), one("^  D2 all")))
}
