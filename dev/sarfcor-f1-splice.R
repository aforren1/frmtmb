# Fix round 1: replace the runtime section of R/autocor-seam.R (from the
# block-entry helper to the fcor() section) with dev/sarfcor-f1-block.txt,
# and let fcor_unwhiten() take a matrix of innovations.
p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/R/autocor-seam.R"
x <- readLines(p)
a <- grep("^#' A block's entry, refused by its package's name", x)
b <- grep("^# -+ fcor\\(\\)$", x)
stopifnot(length(a) == 1, length(b) == 1, a < b)
blk <- readLines("C:/Users/adf44/source/r/frmtmb-wt-sarfcor/dev/sarfcor-f1-block.txt")
x <- c(x[seq_len(a - 1)], blk, x[b:length(x)])
i <- which(x == "  mu + sigma * as.vector(aux[[\"L\"]] %*% u)")
stopifnot(length(i) == 1)
x[i] <- "  mu + sigma * (aux[[\"L\"]] %*% u)"
writeLines(x, p)
