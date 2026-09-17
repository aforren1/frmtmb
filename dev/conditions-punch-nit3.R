# Lane wt-conditions, punch round 1, review NIT 3: why the lane's sweep
# (dev/conditions-sweep-onefile.R, dev/conditions-log/lane-*) counted
# 2,032 expect_error() conditions and the review's
# (dev/conditions-rev-sweep-onefile.R, dev/conditions-rev-log/sweep2-lane)
# counted 2,053 on the same build. Per file, the two counts side by side.
#   Rscript dev/conditions-punch-nit3.R > dev/conditions-punch-log/nit3.txt
lane <- do.call(rbind, lapply(
  list.files(c("dev/conditions-log/lane-core", "dev/conditions-log/lane-ext"),
             "[.]tsv$", full.names = TRUE),
  function(f) utils::read.table(f, sep = "\t", header = TRUE, quote = "\"",
                                comment.char = "", stringsAsFactors = FALSE)))
rev <- do.call(rbind, lapply(
  list.files("dev/conditions-rev-log/sweep2-lane", "[.]rds$",
             full.names = TRUE), readRDS))
rev <- rev[grepl("^expect_error ", rev$key), ]
cat("lane rows", nrow(lane), " review expect_error rows", nrow(rev), "\n")
n_lane <- table(lane$file)
n_rev <- table(rev$file)
files <- union(names(n_lane), names(n_rev))
d <- data.frame(file = files,
                lane = as.integer(n_lane[files]),
                review = as.integer(n_rev[files]))
d[is.na(d)] <- 0L
d <- d[d$lane != d$review, ]
print(d, row.names = FALSE)
cat("files that differ", nrow(d), " review minus lane",
    sum(d$review - d$lane), "\n")
# Which review rows have no lane counterpart, by file and message
for (f in d$file) {
  a <- table(substr(lane$message[lane$file == f], 1, 60))
  b <- table(substr(gsub("\\\\n", " ", rev$message[rev$file == f]), 1, 60))
  cat("\n==", f, "\n")
  k <- union(names(a), names(b))
  diff <- as.integer(b[k]) - as.integer(a[k])
  diff[is.na(diff)] <- as.integer(b[k])[is.na(diff)]
  print(data.frame(message = k, review_minus_lane = diff)[diff != 0, ],
        row.names = FALSE)
}
