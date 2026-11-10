# Summarize a dev/translator-audit.R outdir: how many pairs have latent
# blocks, and every block L1 cannot see at its points. Prints; writes
# nothing.
#
#   Rscript dev/translator-audit-summary.R <outdir>
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
fs <- list.files(outdir, pattern = "^p[0-9]+[.]tsv$", full.names = TRUE)
rows <- lapply(fs, function(p) read.delim(p, stringsAsFactors = FALSE))
allc <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(r) {
  for (m in setdiff(allc, names(r))) r[[m]] <- NA
  r[, allc]
})
a <- do.call(rbind, rows)
logs <- list.files(file.path(outdir, "logs"), pattern = "[.]log$")
ids <- sub("[.]log$", "", logs)
cat("pairs launched", length(ids), "; with a result", length(fs),
    "; no result:", paste(setdiff(ids, sub("[.]tsv$", "", basename(fs))),
                          collapse = " "), "\n")
nob <- a[!nzchar(a$block) | is.na(a$block), ]
cat("pairs without a latent block or not translated:", length(unique(nob$id)),
    "\n")
why <- table(sub(":.*$", "", nob$note))
for (k in names(why)) cat("  ", why[[k]], k, "\n")
odd <- nob[!grepl("^no latent block|^no brms spelling", nob$note), ]
if (nrow(odd)) {
  cat("  of which errors:\n")
  for (i in seq_len(nrow(odd))) {
    cat("   ", odd$id[i], odd$pair[i], ":", substr(odd$note[i], 1, 110), "\n")
  }
}
lb <- a[nzchar(a$block) & !is.na(a$block), ]
cat("pairs with latent blocks:", length(unique(lb$id)), "; blocks:",
    nrow(lb), "; untestable (one value):", sum(is.na(lb$blind)), "\n")
bl <- lb[lb$blind %in% TRUE, ]
cat("BLIND blocks:", nrow(bl), "in", length(unique(bl$id)), "pairs\n")
op <- options(width = 300)
if (nrow(bl)) {
  print(data.frame(id = bl$id, pair = bl$pair, block = bl$block,
                   points = bl$points,
                   max_delta = signif(bl$max_delta, 3),
                   tol = signif(bl$tol, 3), sds = substr(bl$sds, 1, 70)),
        row.names = FALSE, right = FALSE)
}
# the closest calls: seen, but by less than 100 times the tolerance
near <- lb[lb$blind %in% FALSE & lb$max_delta < 100 * lb$tol, ]
cat("\nseen by less than 100 x tol:", nrow(near), "blocks\n")
if (nrow(near)) {
  print(data.frame(id = near$id, block = near$block,
                   ratio = signif(near$max_delta / near$tol, 3),
                   sds = substr(near$sds, 1, 70)), row.names = FALSE,
        right = FALSE)
}
options(op)
