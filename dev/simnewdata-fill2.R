# Punch round 1: paste the round's logs into the findings, verbatim.
#   Rscript dev/simnewdata-fill2.R <punch section file>
f <- "dev/simnewdata-findings.md"
x <- readLines(f)
s <- readLines(commandArgs(trailingOnly = TRUE)[1])
L <- function(p) readLines(file.path("dev/simnewdata-log", p), warn = FALSE)
blk <- function(v) paste0("      ", v)
fills <- list(
  "CACHE-BLOCK" = blk(tail(L("cache-compare.txt"), 1)),
  "TIMING-BLOCK" = blk(c(grep("per-replicate", L("timing-before.txt"), value = TRUE),
                         grep("per-replicate", L("timing-after.txt"), value = TRUE))),
  "RPLOTS2-BLOCK" = paste0("    ", L("rplots-lane.txt")),
  "BOOTCMP-BLOCK" = blk(L("boot-compare.txt")),
  "BOOTLANE-BLOCK" = blk(grep("^lane:", L("boot-lane.txt"), value = TRUE)))
for (k in names(fills)) {
  i <- which(s == k); stopifnot(length(i) == 1L)
  s <- c(s[seq_len(i - 1L)], fills[[k]], s[-seq_len(i)])
}
i <- which(x == "## What I did not do, and why"); stopifnot(length(i) == 1L)
x <- c(x[seq_len(i - 1L)], s, "", x[i:length(x)])
con <- file(f, "wb"); writeLines(x, con, sep = "\n"); close(con)
