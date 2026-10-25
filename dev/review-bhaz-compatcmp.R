# Reviewer: compat statuses, rellib-r12 against the lane build.
L <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-bhaz/dev/review-bhaz-log/compat-lane.rds")
B <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-bhaz/dev/review-bhaz-log/compat-base.rds")
cat("rules base", nrow(B), "lane", nrow(L), "\n")
cat("keys only in lane:", paste(setdiff(L$key, B$key), collapse = "; "), "\n")
cat("keys only in base:", paste(setdiff(B$key, L$key), collapse = "; "), "\n")
m <- merge(B, L, by = "key", suffixes = c(".base", ".lane"))
ch <- m[m$status.base != m$status.lane, ]
cat("shared keys", nrow(m), " status changed", nrow(ch), "\n")
print(ch)
