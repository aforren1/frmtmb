# Lane postfix: the hand-translated vignette rows whose outcome differs
# between rellib-r8 and the lane build (dev/postfix-bv.sh before/after).
#   Rscript dev/postfix-bv-diff.R
rd <- function(arm) {
  fs <- list.files(file.path("dev/postfix-bv-out", arm), "[.]csv$",
                   full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    # the n-th row of each label, so a repeated label keeps its place
    d$key <- paste(d$vignette, d$label, ave(seq_along(d$label), d$vignette,
                                            d$label, FUN = seq_along),
                   sep = " | ")
    d
  }))
}
b <- rd("before")
a <- rd("after")
cat("rows: before", nrow(b), " after", nrow(a), "\n")
m <- merge(b[c("key", "ok", "msg")], a[c("key", "ok", "msg")], by = "key",
           suffixes = c(".before", ".after"), all = TRUE)
ch <- m[(m$ok.before != m$ok.after) %in% TRUE |
          is.na(m$ok.before) | is.na(m$ok.after), ]
cat("rows whose outcome changed:", nrow(ch), "\n")
cat("  FALSE -> TRUE:", sum(!ch$ok.before & ch$ok.after, na.rm = TRUE),
    "\n")
cat("  TRUE -> FALSE:", sum(ch$ok.before & !ch$ok.after, na.rm = TRUE),
    "\n\n")
for (i in seq_len(nrow(ch))) {
  cat(sprintf("%s\n    before ok=%s %s\n    after  ok=%s %s\n", ch$key[i],
              ch$ok.before[i], substr(ch$msg.before[i], 1, 150),
              ch$ok.after[i], substr(ch$msg.after[i], 1, 150)))
}
# failing rows whose message changed (a refusal that now names the
# reason, or a different stop)
st <- m[m$ok.before %in% FALSE & m$ok.after %in% FALSE &
          m$msg.before != m$msg.after, ]
cat("\nrows failing on both builds with a new message:", nrow(st), "\n")
for (i in seq_len(nrow(st))) {
  cat(sprintf("%s\n    before %s\n    after  %s\n", st$key[i],
              substr(st$msg.before[i], 1, 150),
              substr(st$msg.after[i], 1, 150)))
}
