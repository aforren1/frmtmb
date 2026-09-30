# Per-package totals of the lane's full run (dev/ordinal-suite-all) with
# the files rerun after their fix (dev/ordinal-suite-rerun) taking the
# place of their first result. Output: dev/ordinal-p1-log-suite-sum.txt
rd <- function(f) {
  x <- grep(" RESULT ", readLines(f), value = TRUE)
  m <- regmatches(x, regexec(paste0("^(\\S+) RESULT (\\S+) pass=(\\d+) ",
                                    "fail=(\\d+) err=(\\d+) skip=(\\d+) ",
                                    "warn=(\\d+)"), x))
  do.call(rbind, lapply(m, function(v) {
    data.frame(pkg = v[2], file = v[3], pass = as.integer(v[4]),
               fail = as.integer(v[5]), err = as.integer(v[6]),
               skip = as.integer(v[7]), warn = as.integer(v[8]))
  }))
}
a <- rd("dev/ordinal-suite-p1all/p1all.log")
r <- a[0, ]
cat("first run:", nrow(a), "files; rerun:", nrow(r), "files\n")
key <- function(d) paste(d$pkg, d$file)
a <- a[!key(a) %in% key(r), ]
f <- rbind(a, r)
tot <- aggregate(cbind(files = 1L, pass, fail, err, skip, warn) ~ pkg,
                 data = transform(f, files = 1L), FUN = sum)
print(tot, row.names = FALSE)
cat("\nfiles with a failure, an error, a skip or a warning:\n")
print(f[f$fail + f$err + f$skip + f$warn > 0, ], row.names = FALSE)
