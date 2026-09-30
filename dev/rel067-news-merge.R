# Resolve a two-sided NEWS.md conflict block by section: each side's
# text is split at its "## " headings, and the sections are joined
# ours then theirs, one blank line between bullets. The lead text
# before a side's first heading is kept, ours then theirs, for a hand
# edit afterwards.
args <- commandArgs(trailingOnly = TRUE)
f <- args[1]
x <- readLines(f, warn = FALSE)
a <- grep("^<<<<<<< ", x); m <- grep("^=======$", x); b <- grep("^>>>>>>> ", x)
stopifnot(length(a) == 1, length(m) == 1, length(b) == 1)
split_side <- function(v) {
  h <- grep("^## ", v)
  lead <- if (length(h)) v[seq_len(h[1] - 1)] else v
  secs <- list()
  for (i in seq_along(h)) {
    end <- if (i < length(h)) h[i + 1] - 1 else length(v)
    body <- v[seq.int(h[i] + 1, length.out = max(0, end - h[i]))]
    secs[[sub("^## ", "", v[h[i]])]] <- body
  }
  list(lead = lead, secs = secs)
}
bullets <- function(body) {
  s <- grep("^\\* ", body)
  if (!length(s)) return(list())
  lapply(seq_along(s), function(i) {
    end <- if (i < length(s)) s[i + 1] - 1 else length(body)
    y <- body[s[i]:end]
    while (length(y) && y[length(y)] == "") y <- y[-length(y)]
    y
  })
}
# the part above the conflict belongs to the section it opens in
pre <- x[seq_len(a - 1)]
open_head <- tail(grep("^## ", pre), 1)
lead_head <- if (length(open_head)) sub("^## ", "", pre[open_head]) else NA
if (length(open_head)) pre <- pre[seq_len(open_head - 1)]
ours <- x[(a + 1):(m - 1)]; theirs <- x[(m + 1):(b - 1)]
if (!is.na(lead_head)) {
  ours <- c(paste("##", lead_head), ours)
  theirs <- c(paste("##", lead_head), theirs)
}
o <- split_side(ours); t <- split_side(theirs)
keys <- unique(c(names(o$secs), names(t$secs)))
order_first <- c("Breaking changes", "New features", "Formula grammar",
                 "New arguments and answers", "Bug fixes")
keys <- c(intersect(order_first, keys), setdiff(keys, order_first))
out <- c(o$lead, if (length(t$lead) && any(nzchar(t$lead))) c("", t$lead))
for (k in keys) {
  bl <- c(bullets(o$secs[[k]]), bullets(t$secs[[k]]))
  out <- c(out, paste("##", k), "")
  for (y in bl) out <- c(out, y, "")
}
while (length(pre) && pre[length(pre)] == "") pre <- pre[-length(pre)]
res <- c(pre, "", out, x[(b + 1):length(x)])
res <- res[!(c(FALSE, res[-1] == "" & res[-length(res)] == ""))]
writeLines(res, f)
cat(f, ": sections", paste(keys, collapse = " | "), "\n")
