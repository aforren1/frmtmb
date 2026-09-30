set.seed(1); a <- sample.int(12, 1); k <- sample.int(36, 1); cat("g index", a, "g:h index", k, "\n")
fl <- as.vector(t(outer(1:12, 1:3, paste, sep = ":")))
cat("factor-order level k:", levels(interaction(factor(rep(1:12,3)), factor(rep(1:3, each=12)), sep=":", lex.order=TRUE))[k], "\n")
cat("alphabetical level k:", sort(unique(paste(rep(1:12,3), rep(1:3, each=12), sep="_")))[k], "\n")