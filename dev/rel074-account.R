# The plain-tier counts of the merged build against the base, file by
# file: the base is lane mvnlatent's "before" run (rellib-r11, frmtmb
# 0.73.0, the 19b162b9 test trees, NOT_CRAN only, one file per process;
# dev/mvnlatent-findings.md section 5), the merged build this release's
# r074-plain scan. Lists every file whose pass count changed, every file
# only in one run, and the totals per package.
#   Rscript dev/rel074-account.R > dev/rel074-log/account.txt
setwd("C:/Users/adf44/source/r/frmtmb-wt-release")
rd <- function(dir) {
  fs <- list.files(dir, "[.]txt$", full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    x <- readLines(f, warn = FALSE)
    r <- tail(grep("^RESULT ", x, value = TRUE), 1)
    g <- function(k) {
      if (!length(r)) return(NA_integer_)
      as.integer(sub(paste0(".*", k, "=([0-9]+).*"), "\\1", r))
    }
    b <- sub("[.]txt$", "", basename(f))
    data.frame(pkg = sub("--.*$", "", b), file = sub("^.*--", "", b),
               pass = g("pass"), fail = g("fail"), skip = g("skip"))
  }))
}
B <- rd("C:/Users/adf44/source/r/frmtmb-wt-mvnlatent/dev/mvnlatent-out/before")
M <- rd("dev/ciharden-log/r074-plain-suite")
cat("base files", nrow(B), "pass", sum(B$pass), "fail", sum(B$fail),
    "| merged files", nrow(M), "pass", sum(M$pass), "fail", sum(M$fail), "\n")
tb <- function(x) aggregate(cbind(files = 1, pass, fail, skip) ~ pkg, x, sum)
cat("\nbase by package\n"); print(tb(B), row.names = FALSE)
cat("\nmerged by package\n"); print(tb(M), row.names = FALSE)
m <- merge(B, M, by = c("pkg", "file"), all = TRUE,
           suffixes = c(".base", ".merged"))
ch <- m[is.na(m$pass.base) | is.na(m$pass.merged) |
          m$pass.base != m$pass.merged, ]
ch$delta <- ifelse(is.na(ch$pass.merged), 0, ch$pass.merged) -
  ifelse(is.na(ch$pass.base), 0, ch$pass.base)
cat("\nfiles only in the base, only merged, or with another pass count:",
    nrow(ch), "\n")
print(ch[order(ch$pkg, ch$file), c("pkg", "file", "pass.base",
                                   "pass.merged", "delta")],
      row.names = FALSE)
cat("\nsum of deltas:", sum(ch$delta), "\n")
