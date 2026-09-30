# Reviewer: pair the before/after RESULT lines per file and list differences
dir <- "C:/Users/adf44/source/r/frmtmb-wt-formula2/dev/formula2-rev-p1"
pat <- commandArgs(TRUE)[1]   # "" for ungated, "g-" for gated
rd <- function(arm) {
  fs <- list.files(dir, pattern = paste0("^", pat, arm, "-.*[.]log$"),
                   full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    l <- grep("^RESULT", readLines(f, warn = FALSE), value = TRUE)
    nm <- sub(paste0("^", pat, arm, "-"), "", sub("[.]log$", "", basename(f)))
    if (!length(l)) return(data.frame(file = nm, pass = NA, fail = NA,
                                      err = NA, skip = NA, warn = NA))
    g <- function(k) as.integer(sub(paste0(".*", k, "=([0-9]+).*"), "\\1", l[1]))
    data.frame(file = nm, pass = g("pass"), fail = g("fail"), err = g("err"),
               skip = g("skip"), warn = g("warn"))
  }))
}
b <- rd("before"); a <- rd("after")
m <- merge(b, a, by = "file", suffixes = c(".b", ".a"), all = TRUE)
cat(sprintf("files=%d  before: pass=%d fail=%d err=%d skip=%d warn=%d\n",
            nrow(b), sum(b$pass, na.rm = TRUE), sum(b$fail, na.rm = TRUE),
            sum(b$err, na.rm = TRUE), sum(b$skip, na.rm = TRUE),
            sum(b$warn, na.rm = TRUE)))
cat(sprintf("files=%d  after:  pass=%d fail=%d err=%d skip=%d warn=%d\n",
            nrow(a), sum(a$pass, na.rm = TRUE), sum(a$fail, na.rm = TRUE),
            sum(a$err, na.rm = TRUE), sum(a$skip, na.rm = TRUE),
            sum(a$warn, na.rm = TRUE)))
cols <- c("pass", "fail", "err", "skip", "warn")
diff <- apply(m[paste0(cols, ".b")] != m[paste0(cols, ".a")], 1,
              function(r) any(r | is.na(r)))
cat("DIFFERING FILES:\n"); print(m[diff, ], row.names = FALSE)
bad <- m$fail.a > 0 | m$err.a > 0 | m$warn.a > 0 | is.na(m$pass.a)
cat("AFTER NOT CLEAN:\n"); print(m[bad, ], row.names = FALSE)
cat("ALL FILES:\n"); print(m, row.names = FALSE)
