# Review of lane arstep: join fresh-<part>-base.txt and fresh-<part>-lane.txt
# by row. Reference = best of base default, lane default and the lane's
# multi-start references. Short = more than 1e-3 below it; silent = conv 0
# and no warning. Worse = lane below base by more than 1e-3 (and by 1e-6).
# Usage: Rscript review-arstep-sum.R <part>
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
part <- commandArgs(TRUE)[1]
rd <- function(arm) {
  f <- file.path(RV_OUT, sprintf("fresh-%s-%s.txt", part, arm))
  L <- readLines(f)
  cat("log", basename(f), ":", grep("^lib:", L, value = TRUE), "\n")
  fr <- grep("^FR ", L, value = TRUE)
  fr <- fr[!grepl(" ERROR ", fr)]
  sp <- strsplit(fr, " +")
  x <- data.frame(i = as.integer(vapply(sp, `[`, "", 3)),
                  lbl = vapply(sp, function(s) paste(s[5:9], collapse = " "), ""),
                  ll = as.numeric(vapply(sp, function(s) s[which(s == "ll") + 1], "")),
                  conv = as.integer(vapply(sp, function(s) s[which(s == "conv") + 1], "")),
                  w = as.integer(vapply(sp, function(s) s[which(s == "w") + 1], "")),
                  wall = as.numeric(vapply(sp, function(s) s[which(s == "wall") + 1], "")),
                  ac = vapply(sp, function(s) s[which(s == "ac") + 1], ""))
  rr <- grep("^FRR ", L, value = TRUE)
  if (length(rr)) {
    s2 <- strsplit(rr, " +")
    ref <- setNames(as.numeric(vapply(s2, function(s) s[which(s == "best") + 1], "")),
                    vapply(s2, `[`, "", 3))
    x$ref <- ref[as.character(x$i)]
  } else x$ref <- NA
  err <- grep(" ERROR ", grep("^FR ", L, value = TRUE), value = TRUE)
  if (length(err)) cat(paste("ERR", err), sep = "\n")
  x
}
b <- rd("base"); l <- rd("lane")
j <- merge(b, l, by = c("i", "lbl"), suffixes = c(".b", ".l"))
j$best <- pmax(j$ll.b, j$ll.l, j$ref.l, j$ref.b, na.rm = TRUE)
j$gap.b <- j$ll.b - j$best; j$gap.l <- j$ll.l - j$best
sh <- function(g, conv, w) sum(g < -1e-3)
for (grp in c("all", unique(sub(" .*", "", sub("^[A-Z] ", "", j$lbl))))) {
  k <- if (grp == "all") rep(TRUE, nrow(j)) else startsWith(j$lbl, paste0(part, " ", grp, " "))
  if (!any(k)) next
  jj <- j[k, ]
  cat(sprintf("FSUM %s %-8s rows %3d | base short %2d silent %2d worst %8.3g | lane short %2d silent %2d worst %8.3g | lane worse than base (>1e-3) %d (>1e-6) %d, better %d | wall base %.1f lane %.1f\n",
              part, grp, nrow(jj),
              sum(jj$gap.b < -1e-3), sum(jj$gap.b < -1e-3 & jj$conv.b == 0 & jj$w.b == 0), min(jj$gap.b),
              sum(jj$gap.l < -1e-3), sum(jj$gap.l < -1e-3 & jj$conv.l == 0 & jj$w.l == 0), min(jj$gap.l),
              sum(jj$ll.l < jj$ll.b - 1e-3), sum(jj$ll.l < jj$ll.b - 1e-6),
              sum(jj$ll.l > jj$ll.b + 1e-3), sum(jj$wall.b), sum(jj$wall.l)))
}
cat("rows:\n")
for (r in seq_len(nrow(j))) {
  x <- j[r, ]
  tag <- c(if (x$ll.l < x$ll.b - 1e-3) "WORSE", if (x$ll.l > x$ll.b + 1e-3) "better",
           if (x$gap.l < -1e-3) "lane-short", if (x$gap.b < -1e-3) "base-short",
           if (x$w.l != x$w.b) sprintf("warn %d->%d", x$w.b, x$w.l))
  if (length(tag)) {
    cat(sprintf("FROW %3d %-34s base %.4f lane %.4f best %.4f | ac b %s l %s | conv %d/%d w %d/%d | %s\n",
                x$i, x$lbl, x$ll.b, x$ll.l, x$best, x$ac.b, x$ac.l, x$conv.b, x$conv.l,
                x$w.b, x$w.l, paste(tag, collapse = ",")))
  }
}
cat(sprintf("rows joined %d (base %d, lane %d)\n", nrow(j), nrow(b), nrow(l)))
