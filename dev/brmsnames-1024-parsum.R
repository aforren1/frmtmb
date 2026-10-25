## Per-package totals of the lane's run-par.sh tiers, for the findings.
##
##   Rscript dev/brmsnames-1024-parsum.R > dev/brmsnames-1024-log/brmsnames-par-summary.txt
##
## Reads dev/release/brmsnames-<tier>.log and the per-file outputs (for
## the library each file loaded) of the plain and gated tiers, and of
## gated2: the two files the gated tier skipped for the drmTMB and fuzz
## gates, rerun with FRMTMB_DRMTMB_FIT_TESTS and FRMTMB_FUZZ as
## dev/release/run-gated.ps1 sets them.
for (tier in c("plain", "gated")) {
  log <- file.path("dev/release", paste0("brmsnames-", tier, ".log"))
  if (!file.exists(log)) next
  x <- readLines(log)
  cat("== tier", tier, ":", tail(x, 1L), "\n")
  r <- grep("^RESULT", x, value = TRUE)
  pk <- NULL
  rows <- list()
  for (ln in x) {
    if (grepl("^== .* files ==$", ln)) pk <- sub("^== (\\S+) .*$", "\\1", ln)
    if (grepl("^RESULT", ln)) {
      num <- function(k) as.integer(sub(paste0("^.* ", k, "=([0-9]+).*$"),
                                        "\\1", ln))
      rows[[length(rows) + 1L]] <- data.frame(
        pkg = pk, file = sub("^RESULT (\\S+) .*$", "\\1", ln),
        pass = num("pass"), fail = num("fail"), err = num("err"),
        skip = num("skip"), warn = num("warn"))
    }
  }
  d <- do.call(rbind, rows)
  a <- aggregate(cbind(files = 1L, pass, fail, err, skip, warn) ~ pkg,
                 data = transform(d, files = 1L), FUN = sum)
  print(a, row.names = FALSE)
  cat("files with fail, err or warn:", sum(d$fail + d$err + d$warn > 0),
      "\n")
  fs <- list.files(file.path("dev/release",
                             paste0("brmsnames-", tier, "-files")),
                   full.names = TRUE)
  libs <- vapply(fs, function(f) {
    l <- grep("^lib:", readLines(f, warn = FALSE), value = TRUE)[1L]
    if (is.na(l)) "NONE" else sub("^.*; frmtmb \\S+ ", "", l)
  }, "")
  cat("core library of each file:\n")
  print(table(libs, dnn = NULL))
  sk <- d[d$skip > 0, ]
  cat("files with skips:", nrow(sk), "\n\n")
}
