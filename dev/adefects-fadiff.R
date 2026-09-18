# Diff the two arms of dev/adefects-falsealarm.R: which harvested call
# or hypothesis string changed its frmtmb verdict between the base build
# and this lane's. That difference, not the raw disagreement count, is
# what this lane's refusals cost.
#
#   Rscript dev/adefects-fadiff.R > dev/adefects-log/fadiff.txt 2>&1
#
# quote = "" is not optional: the harvest writes with quote = FALSE and
# some calls contain a double quote (`ma("log")`), which read.delim's
# default quoting swallows, shifting the columns and inventing a
# verdict change. That happened on the first run of this script.
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
b <- rd("dev/adefects-log/falsealarm-base.tsv")
l <- rd("dev/adefects-log/falsealarm-lane.tsv")
key <- function(d) paste(d$kind, d$text, sep = "\r")
stopifnot(!anyDuplicated(key(b)), !anyDuplicated(key(l)))
m <- merge(data.frame(k = key(b), kind = b$kind, text = b$text,
                      src = b$src, brms = b$brms, base = b$frm),
           data.frame(k = key(l), lkind = l$kind, ltext = l$text,
                      lsrc = l$src, lbrms = l$brms, lane = l$frm),
           by = "k", all = TRUE)
# a row present in one arm only has NA on the other arm's columns, so
# every field is read from whichever side has it
fill <- is.na(m$kind)
m$kind[fill] <- m$lkind[fill]
m$text[fill] <- m$ltext[fill]
m$src[fill] <- m$lsrc[fill]
m$brms[fill] <- m$lbrms[fill]
both <- !is.na(m$base) & !is.na(m$lane)
ch <- m[both & m$base != m$lane, , drop = FALSE]
cat("rows in base:", nrow(b), " rows in lane:", nrow(l),
    " rows in both arms:", sum(both), "\n")
cat("rows whose frmtmb verdict CHANGED:", nrow(ch), "\n")
for (i in seq_len(nrow(ch))) {
  cat(sprintf("  %-10s brms=%-5s base=%-5s lane=%-5s %s  [%s]\n",
              ch$kind[i], ch$brms[i], ch$base[i], ch$lane[i],
              substr(ch$text[i], 1, 60), ch$src[i]))
}
onlyb <- m[is.na(m$lane), , drop = FALSE]
onlyl <- m[is.na(m$base), , drop = FALSE]
cat("only in the base corpus:", nrow(onlyb),
    "(the spellings this lane rewrote)\n")
for (i in seq_len(nrow(onlyb))) {
  cat(sprintf("   - %-8s brms=%-5s base=%-5s %s\n", onlyb$kind[i],
              onlyb$brms[i], onlyb$base[i], substr(onlyb$text[i], 1, 60)))
}
cat("only in the lane corpus:", nrow(onlyl),
    "(the replacements, and this lane's own test file)\n")
for (i in seq_len(nrow(onlyl))) {
  cat(sprintf("   + %-8s brms=%-5s lane=%-5s %s\n", onlyl$kind[i],
              onlyl$brms[i], onlyl$lane[i], substr(onlyl$text[i], 1, 60)))
}
