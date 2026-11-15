# Lane silent2: how many control fits the 8-ulp "collapsed" test changes
# (lane against base, columns collapsed1..3 of dev/silent2-ordmeasure.R)
sets <- c("forty", "sweep", "steep", "wide")
rd <- function(arm) do.call(rbind, lapply(sets, function(s) {
  x <- utils::read.delim(sprintf("dev/silent2-out/ord-%s-%s.tsv", s, arm))
  x[, c("case", "collapsed1", "collapsed2", "collapsed3")]
}))
B <- rd("base")
L <- rd("lane")
M <- merge(B, L, by = "case", suffixes = c(".b", ".l"))
anyc <- function(sfx) {
  apply(M[, grep(paste0("^collapsed[123][.]", sfx, "$"), names(M))], 1,
        function(r) any(as.logical(r) %in% TRUE))
}
b <- anyc("b")
l <- anyc("l")
cat("fits", nrow(M), "; collapsed on base", sum(b), "; on lane", sum(l),
    "; newly collapsed", sum(l & !b), "\n")
if (any(l & !b)) print(M$case[l & !b])
