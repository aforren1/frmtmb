# Realign the continuation lines of calls the sp_ -> curve_ rename
# lengthened or shortened, so arguments stay under their open paren.
#   Rscript dev/splitsurv-realign.R <file> [<file> ...]
for (f in commandArgs(trailingOnly = TRUE)) {
  x <- readLines(f)
  n <- 0L
  for (i in seq_len(length(x) - 1L)) {
    s <- x[i]
    if (grepl("^ *#", s)) next
    m <- gregexpr("curve_[a-z_]+\\(", s)[[1]]
    if (m[1] < 0) next
    for (k in rev(seq_along(m))) {
      open <- m[k] + attr(m, "match.length")[k] - 1L
      rest <- substr(s, open, nchar(s))
      ch <- strsplit(rest, "")[[1]]
      bal <- cumsum((ch == "(") - (ch == ")"))
      if (tail(bal, 1) > 0 && all(bal > 0)) {
        j <- i + 1L
        depth <- tail(bal, 1)
        while (j <= length(x) && depth > 0) {
          ind <- attr(regexpr("^ *", x[j]), "match.length")
          if (abs(ind - open) <= 4L && ind != open) {
            x[j] <- paste0(strrep(" ", open), sub("^ *", "", x[j]))
            n <- n + 1L
          }
          cj <- strsplit(x[j], "")[[1]]
          depth <- depth + sum(cj == "(") - sum(cj == ")")
          j <- j + 1L
        }
        break
      }
    }
  }
  writeLines(x, f)
  cat(f, n, "lines realigned\n")
}
