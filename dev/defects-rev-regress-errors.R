# Reviewer of lane defects: which items errored in each regression arm.
r <- readRDS("dev/defects-rev-log/regress-lane.rds")
b <- readRDS("dev/defects-rev-log/regress-base.rds")
iserr <- function(v) is.character(v) && length(v) == 1 && grepl("^ERROR", v)
for (m in names(r)) {
  e <- names(Filter(iserr, r[[m]])); eb <- names(Filter(iserr, b[[m]]))
  if (length(e) || length(eb)) {
    cat(m, "| lane:", e, "| base:", eb, "\n")
    for (k in e) cat("    lane", k, ":", substr(r[[m]][[k]], 1, 150), "\n")
  }
}
