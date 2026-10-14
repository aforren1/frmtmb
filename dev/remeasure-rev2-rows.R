# Reviewer round 2: the ordinal rows and the reclassified rows in the
# lane's bin-2 ledger
b <- read.delim("dev/brmsport-bin2-ledger.tsv", quote = "",
                stringsAsFactors = FALSE)
pick <- grepl("^log_lik:4[67][0-9]$|^posterior_predict:3[0-2][0-9]$", b$id)
print(b[pick, c("id", "n", "outcome", "class")], row.names = FALSE)
print(b[b$id %in% c("log_lik:32", "log_lik:302", "log_lik:267",
                    "distributions:238", "distributions:446"),
        c("id", "n", "outcome", "class")], row.names = FALSE)
