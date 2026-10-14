# Reviewer listing: the bin-2T passes and bin-2P passes with their
# assertion text, from the lane's bin-2 ledger.
b <- read.delim("dev/brmsport-bin2-ledger.tsv", quote = "",
                stringsAsFactors = FALSE)
x <- b[b$bin == "2T" & b$outcome == "pass", ]
for (i in seq_len(nrow(x))) {
  cat(x$id[i], "|", x$held[i], "|", substr(x$assertion[i], 1, 110), "\n")
}
y <- b[b$bin == "2P" & b$outcome == "pass", ]
kind <- ifelse(grepl("length|dim", y$assertion), "length/dim", "value")
cat("\n2P passes:", nrow(y), "by kind:\n"); print(table(kind))
cat("\n2P value passes:\n")
for (i in which(kind == "value")) cat(y$id[i], "|", y$assertion[i], "\n")
