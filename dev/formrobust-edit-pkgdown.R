# One-off edit of _pkgdown.yml: the new topics. Kept as the record.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/_pkgdown.yml"
x <- readLines(p)
i <- which(x == "  - frmtmb-autocor"); stopifnot(length(i) == 1L)
x <- append(x, "  - autocor-terms", i)
j <- which(x == "  - autocor_matrix"); stopifnot(length(j) == 1L)
x <- append(x, "  - autocor", j)
crlf <- any(grepl("\r", readChar(p, 2000, useBytes = TRUE)))
con <- file(p, "wb"); writeLines(x, con, sep = if (crlf) "\r\n" else "\n")
close(con)
