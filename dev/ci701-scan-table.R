# Lane ci701: per-package totals of one scan run (dev/ci701-scan.sh),
# with the frmtmb library every file loaded and the trace counts.
#
#   Rscript dev/ci701-scan-table.R <run> [<run> ...]
runs <- commandArgs(TRUE)
logd <- "C:/Users/adf44/source/r/frmtmb-wt-ci701/dev/ci701-log"
for (r in runs) {
  fs <- list.files(file.path(logd, paste0(r, "-suite")), pattern = "[.]txt$",
                   full.names = TRUE)
  rows <- lapply(fs, function(f) {
    l <- readLines(f, warn = FALSE)
    res <- grep("^RESULT ", l, value = TRUE)
    lib <- grep("^lib: ", l, value = TRUE)
    num <- function(k) {
      if (!length(res)) return(NA_integer_)
      as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", res[1]))
    }
    data.frame(pkg = sub("--.*", "", basename(f)),
               pass = num("pass"), fail = num("fail"), err = num("err"),
               skip = num("skip"), warn = num("warn"),
               corelib = sub(".*; frmtmb ", "", lib[1]),
               cor = sum(startsWith(l, "CI701COR")),
               edge = sum(startsWith(l, "CI701EDGE")))
  })
  d <- do.call(rbind, rows)
  cat("RUN", r, ": files", nrow(d), ", without a RESULT line",
      sum(is.na(d$pass)), "\n")
  cat("  core library on every file:", paste(unique(d$corelib),
                                              collapse = " | "), "\n")
  a <- aggregate(cbind(files = 1, pass, fail, err, skip, warn, cor, edge) ~
                   pkg, data = d, FUN = sum)
  print(a, row.names = FALSE)
}
