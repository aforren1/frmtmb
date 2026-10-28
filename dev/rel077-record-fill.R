# Paste the generated blocks into dev/round-20261028.md: each marker
# line is replaced by the file it names, between BEGIN and END
# GENERATED comments, so the record's counts come from the logs.
#   Rscript dev/rel077-record-fill.R
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel077")
f <- "dev/round-20261028.md"
x <- sub("\r$", "", readLines(f, warn = FALSE))
blocks <- list(
  "PENDING-COUNTS" = c("dev/rel077-log/counts-plain.md",
                       "dev/rel077-log/counts-ref.md"),
  "PENDING-ACCOUNT-BLOCK" = "dev/rel077-log/account.txt",
  "PENDING-MAP" = "dev/rel077-news-map.md")
for (k in names(blocks)) {
  i <- which(x == k)
  if (length(i) != 1L) stop(k, ": expected one marker, found ", length(i))
  body <- unlist(lapply(blocks[[k]], function(p) {
    y <- sub("\r$", "", readLines(p, warn = FALSE))
    while (length(y) && !nzchar(y[length(y)])) y <- y[-length(y)]
    c(paste0("<!-- BEGIN GENERATED: ", p, " -->"), "", y, "",
      "<!-- END GENERATED -->", "")
  }))
  x <- c(x[seq_len(i - 1L)], body[-length(body)], x[-seq_len(i)])
  cat("filled", k, "\n")
}
con <- file(f, "wb"); writeLines(x, con, sep = "\n"); close(con)
