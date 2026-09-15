## The follow-up numbers dev/coh-summarize.R does not print: what the
## three non-convergent replicates do to the headline coverage, the
## paired win-loss between the correct model and the survey's, and the
## variance components the correct model recovers.
##
## Run: Rscript dev/coh-extra.R dev/coh-recovery-main.tsv
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- "dev/coh-partial.tsv"
src <- readLines("dev/coh-summarize.R")
stop_at <- grep("^truth <- ", src)[1L]
eval(parse(text = paste(src[seq_len(stop_at - 1L)], collapse = "\n")))

d <- read_tsv_kv(args[1L])
n_rung <- length(unique(d$rung))
keep <- names(which(table(d$seed) == n_rung))
d <- d[as.character(d$seed) %in% keep, ]
f <- d[d$rung == "full", ]
i <- d[d$rung == "id", ]
cat("complete replicates:", length(keep), "\n")
cat("full: covered", sum(f$covers), "of", nrow(f), "=",
    signif(mean(f$covers), 4), "\n")
g <- f[f$conv == 0, ]
cat("full, convergence 0 only: covered", sum(g$covers), "of", nrow(g),
    "=", signif(mean(g$covers), 4), "\n")
nc <- f[f$conv != 0, ]
if (nrow(nc)) {
  cat("the", nrow(nc), "non-convergent replicates: seeds",
      paste(nc$seed, collapse = " "), "; covers",
      paste(nc$covers, collapse = " "), "; max gradient",
      paste(signif(nc$maxgrad, 3), collapse = " "), "\n")
}
cat("mean width: id", signif(mean(i$width), 4), " full",
    signif(mean(f$width), 4), " ratio",
    signif(mean(i$width) / mean(f$width), 4), "\n")
cat("full sd_id", signif(mean(f$sd_id), 4), "+/-",
    signif(stats::sd(f$sd_id), 3), " sd_idcond",
    signif(mean(f$sd_idcond), 4), "+/-",
    signif(stats::sd(f$sd_idcond), 3), "\n")
m <- merge(i[, c("seed", "width", "covers")],
           f[, c("seed", "width", "covers")], by = "seed",
           suffixes = c("_id", "_full"))
cat("id missed and full covered:", sum(!m$covers_id & m$covers_full),
    "of", nrow(m), "\n")
cat("id covered and full missed:", sum(m$covers_id & !m$covers_full),
    "\n")
