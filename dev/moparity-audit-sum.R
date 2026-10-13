# dev/moparity-audit.sh's records: per package, files, objective calls
# with a theta, and the smallest theta evaluated; every file below the
# floor's threshold (log sd -308.45).
#   Rscript dev/moparity-audit-sum.R
x <- readLines("dev/moparity-audit-out/all.txt")
f <- function(k) sub(paste0(".* ", k, " ([^ ]+).*"), "\\1", x)
d <- data.frame(pkg = sub("^AUDIT ([^ ]+) .*", "\\1", x),
                file = sub("^AUDIT [^ ]+ ([^ ]+) .*", "\\1", x),
                calls = as.numeric(f("calls")),
                min = as.numeric(f("min_theta")),
                fail = f("fail"))
cat("files with a record:", nrow(d), "| distinct", length(unique(paste(
  d$pkg, d$file))), "\n")
a <- do.call(rbind, lapply(split(d, d$pkg), function(s) {
  data.frame(pkg = s$pkg[1], files = nrow(s), with_theta = sum(s$calls > 0),
             calls = sum(s$calls), min_theta = min(s$min),
             files_failing = sum(s$fail != "0"))
}))
print(a, row.names = FALSE)
cat("smallest theta over every file:", min(d$min), "in",
    d$pkg[which.min(d$min)], d$file[which.min(d$min)], "\n")
cat("files below -308.45:", sum(d$min < -308.45), "\n")
o <- d[order(d$min), ][1:5, ]
print(o, row.names = FALSE)
