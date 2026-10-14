# Reviewer round 3: the ledger's ported flag against outcome and class
b <- read.delim("dev/brmsport-bin2-ledger.tsv", quote = "",
                stringsAsFactors = FALSE)
p <- b[b$bin == "2P", ]
print(tapply(p$n, list(paste(p$outcome, p$class), p$ported), sum))
cat("\nnot-ported rows:\n")
x <- p[p$class == "not-ported", c("id", "n", "ported", "block", "reason")]
x$block <- substr(x$block, 1, 50); x$reason <- substr(x$reason, 1, 110)
print(x, row.names = FALSE)
cat("\nbrms-backends and a sample of brms-internal rows:\n")
y <- p[p$class %in% c("brms-backends", "brms-internal"),
       c("id", "n", "class", "block")]
y$block <- substr(y$block, 1, 70)
print(y, row.names = FALSE)
cat("\nG1 passes (distributions:*):\n")
z <- p[grepl("^distributions", p$id) & p$outcome == "pass",
       c("id", "n", "kind", "assertion")]
z$assertion <- substr(z$assertion, 1, 80)
print(z, row.names = FALSE)
