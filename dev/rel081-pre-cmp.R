# Compare pairs remeasured on the merged build with each lane's own
# outdir, column by column (seconds and the build's own path excepted).
#
#   Rscript dev/rel081-pre-cmp.R <outdir> <id>=<lane outdir> ...
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
rd <- function(p) utils::read.delim(p, colClasses = "character",
                                    na.strings = character())
skip <- c("seconds", "frmtmb", "frm_warn")
for (a in args[-1]) {
  id <- sub("=.*", "", a)
  lane <- sub("^[^=]*=", "", a)
  f <- file.path(out, paste0(id, ".tsv"))
  if (!file.exists(f)) { cat(id, "NO RESULT\n"); next }
  m <- rd(f)
  l <- rd(file.path(lane, paste0(id, ".tsv")))
  cols <- setdiff(intersect(names(m), names(l)), skip)
  dif <- cols[vapply(cols, function(k) !identical(m[[k]], l[[k]]), NA)]
  cat(sprintf("%s merged: %s %s const %s grad %s at %s live %s | lane: %s %s const %s grad %s | differ: %s%s\n",
              id, m$verdict, m$level, m$l1_const, m$l1_grad, m$l1_at,
              if ("live_status" %in% names(m)) m$live_status else "-",
              l$verdict, l$level, l$l1_const, l$l1_grad,
              if (length(dif)) paste(dif, collapse = ",") else "none",
              if (length(setdiff(names(m), names(l))))
                paste0(" | only merged: ",
                       paste(setdiff(names(m), names(l)), collapse = ","))
              else ""))
}
