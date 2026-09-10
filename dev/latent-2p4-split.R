# Lane `latent`, item 2.4: split the 200 replicates into the ones where
# lca() and poLCA reached the same optimum and the ones where they did
# not, so the recovery table can be read both ways.
#
# A recovery table that silently averages over a local optimum is not a
# recovery table; a recovery table that silently drops it is worse. Both
# files are written and dev/latent-findings.md reports both.
#
#   Rscript dev/latent-2p4-split.R

source("dev/latent-env.R")
f <- "dev/latent-2p4-lca.tsv"
d <- utils::read.delim(f, check.names = FALSE)
g <- d$ll_frm - d$ll_polca
bad <- abs(g) > 1e-6 * abs(d$ll_frm)
cat("replicates:", nrow(d), "  same optimum:", sum(!bad),
    "  different:", sum(bad), "\n")
cat("gap on the different ones:",
    formatC(range(g[bad]), digits = 6, format = "f"), "\n")
utils::write.table(d[!bad, , drop = FALSE],
                   "dev/latent-2p4-lca-agree.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)
utils::write.table(d[bad, , drop = FALSE],
                   "dev/latent-2p4-lca-disagree.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)
cat("wrote dev/latent-2p4-lca-agree.tsv and",
    "dev/latent-2p4-lca-disagree.tsv\n")
