# Move a function's mutable locals into an environment, by parse data.
#
# Every SYMBOL token named in `syms` between lines `from` and `to` of
# `file` becomes `<env>$<sym>`, and every `<<-` in that range becomes
# `<-`. Parse data, not a regular expression, so a name inside a string
# or a comment, an argument name (`f(M = 1)`), a called function and the
# right side of `x$M` are left as they are. Line endings are kept.
#
# Usage: Rscript dev/superassign-rename.R <file> <from> <to> <env> <sym>...

a <- commandArgs(trailingOnly = TRUE)
file <- a[1L]
from <- as.integer(a[2L])
to <- as.integer(a[3L])
env <- a[4L]
syms <- a[-(1:4)]

raw <- readBin(file, "raw", file.info(file)$size)
txt <- rawToChar(raw)
crlf <- grepl("\r\n", txt, fixed = TRUE)
lines <- strsplit(gsub("\r\n", "\n", txt, fixed = TRUE), "\n",
                  fixed = TRUE)[[1L]]

pd <- utils::getParseData(parse(text = lines, keep.source = TRUE))
pd <- pd[pd$terminal, , drop = FALSE]
pd <- pd[order(pd$line1, pd$col1), , drop = FALSE]
prev <- c("", pd$token[-nrow(pd)])
inr <- pd$line1 >= from & pd$line1 <= to
hit_sym <- inr & pd$token == "SYMBOL" & pd$text %in% syms &
  !prev %in% c("'$'", "'@'")
hit_ga <- inr & pd$token == "LEFT_ASSIGN" & pd$text == "<<-"

edits <- rbind(
  data.frame(line = pd$line1[hit_sym], col = pd$col1[hit_sym],
             len = nchar(pd$text[hit_sym]),
             new = paste0(env, "$", pd$text[hit_sym])),
  data.frame(line = pd$line1[hit_ga], col = pd$col1[hit_ga], len = 3L,
             new = rep("<-", sum(hit_ga))))
# right to left within a line, so earlier columns stay valid
edits <- edits[order(edits$line, -edits$col), , drop = FALSE]
for (i in seq_len(nrow(edits))) {
  l <- lines[edits$line[i]]
  c0 <- edits$col[i]
  lines[edits$line[i]] <- paste0(substr(l, 1L, c0 - 1L), edits$new[i],
                                 substr(l, c0 + edits$len[i], nchar(l)))
}
out <- paste0(paste(lines, collapse = if (crlf) "\r\n" else "\n"),
              if (crlf) "\r\n" else "\n")
writeBin(charToRaw(out), file)
cat(sprintf("%s: %d symbols, %d `<<-`\n", file, sum(hit_sym), sum(hit_ga)))
