f <- readRDS("dev/review-spatial-log/fc-after.rds")
for (m in names(f$fits)) cat(m, vapply(f$fits[[m]], function(x) sprintf("%d/%d", sum(is.finite(unlist(x))), length(unlist(x))), ""), "\n")
cat(names(f$fits[[1]]), "\n")
