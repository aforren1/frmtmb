# Reviewer recount of the 0.70.0 headline from the reviewer's own rerun
# (dev/remeasure-rev-wt/dev/remeasure-out/r8), straight from the
# per-vignette records and not through the lane's reducers.
d <- "dev/remeasure-rev-wt/dev/remeasure-out/r8"
fs <- list.files(file.path(d, "bv"), "[.]csv$", full.names = TRUE)
fs <- fs[basename(fs) != "drift.csv"]
H <- do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
cat("hand rows by kind:\n"); print(table(H$kind, H$ok))
nd <- H[H$kind != "data", ]
cat("non-data hand rows", nrow(nd), "ok", sum(nd$ok), "\n")
M <- readRDS(file.path(d, "port", "results-merged.rds"))
cat("merged columns:", paste(names(M), collapse = " "), "\n")
for (k in c("model", "post", "other")) {
  m <- M[M$kind == k, ]
  cat(sprintf("%-6s n %3d raw OK %3d spell OK %3d\n", k, nrow(m),
              sum(m$status_raw == "OK"), sum(m$status_spell == "OK")))
  print(table(m$class))
}
