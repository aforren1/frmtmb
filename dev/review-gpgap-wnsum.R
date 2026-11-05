# Short check of fix round 1, item 3: the WN lines as a confusion table.
#   Rscript dev/review-gpgap-wnsum.R
L <- unlist(lapply(Sys.glob("dev/review-gpgap-out/wn-*.txt"), readLines, warn = FALSE))
L <- L[startsWith(L, "WN ")]
bad <- L[grepl("ERROR|REF FAILED", L)]
ok <- L[!grepl("ERROR|REF FAILED", L)]
g <- function(p, l = ok) sub(p, "\\1", l)
tb <- data.frame(
  tag = g("^WN (.*?) \\| m .*"),
  conv = as.numeric(g(".* conv (\\d) warned.*")),
  warned = grepl("warned TRUE", ok),
  r = g(".*r90/thr (\\S+) \\|.*"),
  gap = as.numeric(g(".* gap (\\S+) dzmax.*")),
  dz = as.numeric(g(".*dzmax (\\S+) moved.*")),
  moved = grepl("moved TRUE", ok),
  refconv = as.numeric(g(".*refconv (\\d+) refgr.*")),
  refgr = as.numeric(g(".*refgr (\\S+)$")),
  stringsAsFactors = FALSE)
cat(sprintf("designs %d, errors or no reference %d\n", length(L), length(bad)))
if (length(bad)) writeLines(substr(bad, 1, 150))
fam <- sub(" .*", "", tb$tag)
cat(sprintf("moved %d: fire %d, miss %d | unmoved %d: fire %d\n", sum(tb$moved),
            sum(tb$moved & tb$warned), sum(tb$moved & !tb$warned), sum(!tb$moved),
            sum(!tb$moved & tb$warned)))
cat("by family (moved fire/moved, unmoved fire/unmoved):\n")
for (f in unique(fam)) {
  z <- tb[fam == f, ]
  cat(sprintf("  %-12s %d/%d, %d/%d\n", f, sum(z$moved & z$warned), sum(z$moved),
              sum(!z$moved & z$warned), sum(!z$moved)))
}
options(width = 170)
cat("MISSES:\n")
print(tb[tb$moved & !tb$warned, c("tag", "conv", "r", "gap", "dz", "refconv", "refgr")],
      row.names = FALSE)
cat("FALSE ALARMS:\n")
print(tb[!tb$moved & tb$warned, c("tag", "conv", "r", "gap", "dz", "refconv", "refgr")],
      row.names = FALSE)
cat("unread r (NA):", sum(tb$r == "NA"), "\n")
cat("frm() above the reference by more than 0.5:\n")
print(tb[tb$gap < -0.5, c("tag", "conv", "r", "gap", "dz", "refconv", "refgr", "warned")],
      row.names = FALSE)
