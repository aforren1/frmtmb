# Reviewer of lane optrap, item 8: dev/optrap-rev-prec.R's four runs side
# by side (lane and rellib-r7, pinned to P and to E cores).
X <- lapply(c("laneP", "r7P", "laneE", "r7E"), function(b) {
  read.delim(sprintf("dev/optrap-rev-log/prec-%s.tsv", b))
})
m <- X[[1]][, c("case", "what")]
for (x in X) m[[paste0(x$build[1], "_abs")]] <- signif(x$max_abs, 3)
for (x in X) m[[paste0(x$build[1], "_se")]] <- signif(x$max_se_units, 3)
print(m, right = FALSE)
