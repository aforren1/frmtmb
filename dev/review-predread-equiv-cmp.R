# Item 9: all 172 equivalence pairs on the lane's build against the baseline of 0.78.0,
# with every number of each pair, not only the verdict. Run from the worktree root.
base <- "C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out"
out <- "dev/review-predread-equiv-out"
a <- read.delim(file.path(base, "results.tsv"), stringsAsFactors = FALSE)
b <- read.delim(file.path(out, "results.tsv"), stringsAsFactors = FALSE)
cat("pairs", nrow(a), nrow(b), " identical ids", identical(a$id, b$id), "\n")
reached <- readLines("dev/predread-log/equiv-reached.txt")
cat("lane's reached set:", length(reached), "\n")
cat("columns:", paste(names(a), collapse = " "), "\n")
vc <- a$verdict != b$verdict
cat("verdict changes:", sum(vc), "\n")
print(data.frame(a[vc, c("id", "feature_a", "feature_b")], was = a$verdict[vc], now = b$verdict[vc],
                 in_lane_set = a$id[vc] %in% reached), row.names = FALSE)
vm <- a$verdict_measured != b$verdict_measured
cat("measured verdict changes:", sum(vm), "\n")
print(data.frame(a[vm, c("id", "feature_a", "feature_b")], was = a$verdict_measured[vm], now = b$verdict_measured[vm],
                 in_lane_set = a$id[vm] %in% reached), row.names = FALSE)
num <- names(a)[vapply(a, is.numeric, NA)]
moved <- list()
for (cn in num) {
  x <- a[[cn]]
  y <- b[[cn]]
  d <- (is.na(x) != is.na(y)) | (!is.na(x) & !is.na(y) & abs(x - y) > 1e-8 * pmax(1, abs(x), abs(y)))
  d[is.infinite(x) & is.infinite(y)] <- FALSE
  if (any(d)) moved[[cn]] <- data.frame(col = cn, id = a$id[d], was = x[d], now = y[d], in_lane_set = a$id[d] %in% reached)
}
mv <- do.call(rbind, moved)
if (!is.null(mv)) {
  cat("numeric cells that moved:", nrow(mv), " pairs ", length(unique(mv$id)), " outside the lane's set: ",
      length(unique(mv$id[!mv$in_lane_set])), "\n")
  print(mv[!mv$in_lane_set, ], row.names = FALSE)
  print(head(mv[mv$in_lane_set, ], 30), row.names = FALSE)
}
out_set <- setdiff(a$id, reached)
pick <- function(d) d$level
cat("levels outside the lane's set, changes:", sum(a$level[a$id %in% out_set] %in% NA != b$level[a$id %in% out_set] %in% NA), "\n")
for (d in c(base, out)) {
  s <- readLines(file.path(d, "summary.md"))
  cat(basename(d), ":", s[grepl("agree with brms [(]the headline", s)], "\n")
}
failed <- grep("FAILED", list.files(file.path(out, "logs"), full.names = TRUE), value = TRUE)
lg <- list.files(file.path(out, "logs"), full.names = TRUE)
cat("logs with FAILED:", sum(vapply(lg, function(f) any(grepl("^FAILED", readLines(f, warn = FALSE))), NA)), "\n")
