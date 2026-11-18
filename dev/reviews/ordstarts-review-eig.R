O <- "dev/ordstarts-out/"
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE, quote = "", comment.char = "")
X <- do.call(rbind, lapply(Sys.glob(paste0(O, "main-*.tsv")), rd))
X$key <- paste(X$family, X$re, X$n, X$seed)
A <- X[X$arm == "default", ]; B <- X[X$arm == "starts10", ][match(A$key, X$key[X$arm=="starts10"]), ]
ok <- A$err == "" | is.na(A$err)
cat("mineig A distribution (non-error, code 0):\n")
s <- ok & A$code %in% 0
print(table(cut(A$mineig[s], c(-Inf, -1e-3, 0, 1e-6, 1e-4, 1e-2, Inf))))
print(table(cut(A$mineig[s], c(-Inf, -1e-3, 0, 1e-6, 1e-4, 1e-2, Inf)), A$nlost[s]))
cat("maxse vs mineig, A code0 by design (median mineig):\n")
print(tapply(A$mineig[s], paste(A$family, A$re, A$n)[s], function(v) signif(quantile(v, c(0, .25, .5), na.rm=TRUE), 2)))
d <- B$ll - A$ll; up <- d > 1e-3 & !is.na(d)
print(B[up, c("key", "mineig", "maxse", "nlost", "tse1", "tse2", "se_deferred")])
