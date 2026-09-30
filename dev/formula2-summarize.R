# Count the suite logs in dev/formula2-suite/: one RESULT line per file
# is expected; a log without one is reported as MISSING, never counted.
d <- "C:/Users/adf44/source/r/frmtmb-wt-formula2/dev/formula2-suite"
logs <- list.files(d, pattern = "[.]log$", full.names = TRUE)
logs <- logs[basename(logs) != "driver.log"]
rows <- lapply(logs, function(f) {
  l <- grep("^RESULT ", readLines(f, warn = FALSE), value = TRUE)
  if (!length(l)) return(data.frame(log = basename(f), pkg = NA, file = NA,
                                    pass = NA, fail = NA, err = NA,
                                    skip = NA, warn = NA))
  m <- regmatches(l, regexec(paste0("RESULT (\\S+) (\\S+) pass=(\\d+) ",
                                    "fail=(\\d+) err=(\\d+) skip=(\\d+) ",
                                    "warn=(\\d+)"), l))[[1]]
  data.frame(log = basename(f), pkg = m[2], file = m[3],
             pass = as.integer(m[4]), fail = as.integer(m[5]),
             err = as.integer(m[6]), skip = as.integer(m[7]),
             warn = as.integer(m[8]))
})
tab <- do.call(rbind, rows)
for (p in unique(na.omit(tab$pkg))) {
  t <- tab[tab$pkg %in% p, ]
  cat(sprintf("SUITE %s files=%d pass=%d fail=%d err=%d skip=%d warn=%d\n",
              p, nrow(t), sum(t$pass), sum(t$fail), sum(t$err),
              sum(t$skip), sum(t$warn)))
}
bad <- tab[is.na(tab$pass) | tab$fail > 0 | tab$err > 0 | tab$warn > 0, ]
if (nrow(bad)) {
  cat("NOT CLEAN:\n")
  print(bad, row.names = FALSE)
}
sk <- tab[!is.na(tab$skip) & tab$skip > 0, c("pkg", "file", "skip")]
if (nrow(sk)) {
  cat("WITH SKIPS:\n")
  print(sk, row.names = FALSE)
}
