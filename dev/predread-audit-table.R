# The audit matrix as markdown, generated from the logs of dev/predread-audit.R
# and dev/predread-audit-sample.R: A agrees, x differs, - no reading (a side
# refuses, or the reader has no analog on this shape). Writes the block pasted
# into dev/predread-findings.md.
#   Rscript dev/predread-audit-table.R > dev/predread-log/audit-table.md
tab <- function(f) {
  r <- read.delim(f, stringsAsFactors = FALSE, quote = "")
  r$mark <- ifelse(r$ok %in% TRUE, "A", ifelse(r$ok %in% FALSE, "x", "-"))
  r
}
show <- function(r) {
  sh <- unique(r$shape)
  rd <- unique(r$reader)
  cat("| reader |", paste(sh, collapse = " | "), "|\n")
  cat("|---|", paste(rep("---", length(sh)), collapse = "|"), "|\n", sep = "")
  for (k in rd) {
    v <- vapply(sh, function(s) {
      m <- r$mark[r$shape == s & r$reader == k]
      if (length(m)) m[1L] else ""
    }, "")
    cat("| ", gsub("[|]", "/", k), " | ", paste(v, collapse = " | "), " |\n",
        sep = "")
  }
}
cnt <- function(r, lab) {
  cat(sprintf("%s: %d cells; %d agree, %d differ, %d no reading.\n\n", lab,
              nrow(r), sum(r$ok %in% TRUE), sum(r$ok %in% FALSE),
              sum(is.na(r$ok))))
}
for (a in list(c("core", "dev/predread-log/audit-%s.tsv"),
               c("frmtmb.sample", "dev/predread-log/audit-sample-%s.tsv"))) {
  cat("###", a[1L], "after the change\n\n")
  af <- tab(sprintf(a[2L], "after"))
  be <- tab(sprintf(a[2L], "before"))
  cnt(be, "before (rellib-r16)")
  cnt(af, "after")
  show(af)
  m <- merge(be[c("shape", "reader", "mark", "metric")],
             af[c("shape", "reader", "mark", "metric")],
             by = c("shape", "reader"), suffixes = c("_was", "_now"))
  ch <- m[m$mark_was != m$mark_now, ]
  cat("\nCells that changed (before to after; metric is the relative",
      "difference, or the ratio of a z to its bound):\n\n")
  cat("| shape | reader | was | now | metric was | metric now |\n|---|---|---|---|---|---|\n")
  for (i in seq_len(nrow(ch))) {
    cat(sprintf("| %s | %s | %s | %s | %s | %s |\n", ch$shape[i], ch$reader[i],
                ch$mark_was[i], ch$mark_now[i], format(ch$metric_was[i]),
                format(ch$metric_now[i])))
  }
  cat("\nCells that still differ:\n\n")
  d <- af[af$ok %in% FALSE, ]
  for (i in seq_len(nrow(d))) {
    cat(sprintf("- %s, %s: metric %s (%s)\n", d$shape[i], d$reader[i],
                format(d$metric[i]), d$note[i]))
  }
  cat("\n")
}
