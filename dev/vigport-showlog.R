# Print the logged output of the bv() rows whose drift is given, so a
# RUNS-NOW or FAILS-NOW row can be read rather than counted.
#   Rscript dev/vigport-showlog.R <bv_out_dir> <drift> [max_lines]
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
want <- args[2]
mx <- if (length(args) > 2) as.integer(args[3]) else 25L
d <- utils::read.csv(file.path(out, "drift.csv"), stringsAsFactors = FALSE)
d <- d[d$drift == want, ]
for (i in seq_len(nrow(d))) {
  lg <- readLines(file.path(out, paste0("bv-", d$vignette[i], ".log")),
                  warn = FALSE)
  hdr <- which(startsWith(lg, paste0("--- [", d$kind[i], "] ", d$label[i])))
  cat("\n=====", d$vignette[i], "/", d$label[i], "\n")
  for (h in hdr) {
    end <- h + which(grepl("^\\[(ok|err) ", lg[(h + 1):length(lg)]))[1]
    body <- lg[(h + 1):end]
    if (length(body) > mx) body <- c(head(body, mx), "  ...")
    cat(body, sep = "\n")
  }
}
