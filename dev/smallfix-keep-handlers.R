# Rewrites the error handlers that dev/smallfix-handlers.R lists so they
# signal the refusal for a part of the model again (keep_missing(),
# R/extension-home.R). Only handlers inside the listed functions are
# touched, and only the simple fallback forms; anything else is printed
# for a hand edit.
#   Rscript dev/smallfix-keep-handlers.R <handlers.txt>
# Run once; it edits R/ in place.
a <- commandArgs(trailingOnly = TRUE)
x <- readLines(a[1])
x <- x[grepl("^ R/", x)]
lst <- unique(data.frame(file = sub("^ (R/\\S+).*", "\\1", x),
                         fun = sub("^ R/\\S+\\s+(\\S+).*", "\\1", x)))
forms <- c(
  "error = function\\(e\\) NULL" = "error = keep_missing(NULL)",
  "error = function\\(err\\) NULL" = "error = keep_missing(NULL)",
  "error = function\\(e\\) integer\\(0\\)" = "error = keep_missing(integer(0))",
  "error = function\\(e\\) NA_real_" = "error = keep_missing(NA_real_)",
  "error = function\\(e\\) e([),]|$)" = "error = keep_missing(identity)\\1",
  "error = identity([),]|$)" = "error = keep_missing(identity)\\1",
  "error = function\\(e\\) fit\\$cache\\$sdr" =
    "error = keep_missing(function(e) fit$cache$sdr)")
for (f in unique(lst$file)) {
  src <- readLines(f)
  ex <- parse(f, keep.source = TRUE)
  refs <- attr(ex, "srcref")
  funs <- lst$fun[lst$file == f]
  n0 <- 0L
  for (k in seq_along(ex)) {
    e <- ex[[k]]
    if (!(is.call(e) && identical(e[[1L]], as.name("<-")) &&
          deparse(e[[2L]]) %in% funs)) next
    rows <- refs[[k]][1L]:refs[[k]][3L]
    for (i in rows) {
      for (p in names(forms)) {
        if (grepl(p, src[i])) {
          src[i] <- gsub(p, forms[[p]], src[i])
          n0 <- n0 + 1L
        }
      }
    }
  }
  writeLines(src, f)
  cat(sprintf("%-26s %2d handler lines rewritten\n", f, n0))
}
