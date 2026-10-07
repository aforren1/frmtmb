# Reviewer of the 0.69.0 consolidation: per-package totals of a
# dev/relrev069-scan.sh run, and every file with a non-pass.
#   Rscript dev/relrev069-scan-sum.R <run>
run <- commandArgs(trailingOnly = TRUE)[1]
d <- file.path("C:/Users/adf44/source/r/frmtmb-wt-release/dev/relrev069-log",
               paste0(run, "-suite"))
fs <- list.files(d, pattern = "[.]txt$", full.names = TRUE)
rows <- lapply(fs, function(f) {
  l <- readLines(f, warn = FALSE)
  r <- grep("^RESULT", l, value = TRUE)
  lib <- grep("^lib:", l, value = TRUE)
  p <- sub("--.*", "", basename(f))
  if (!length(r)) return(data.frame(pkg = p, file = basename(f), pass = NA,
                                    fail = NA, err = NA, skip = NA, warn = NA,
                                    libok = FALSE))
  r <- r[length(r)]
  g <- function(k) as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", r))
  data.frame(pkg = p, file = basename(f), pass = g("pass"), fail = g("fail"),
             err = g("err"), skip = g("skip"), warn = g("warn"),
             libok = length(lib) == 1 &&
               lengths(regmatches(lib, gregexpr("rellib-r7", lib))) == 2)
})
x <- do.call(rbind, rows)
cat("files:", nrow(x), "with RESULT:", sum(!is.na(x$pass)),
    "lib lines naming rellib-r7 twice:", sum(x$libok), "\n\n")
cat("| package | files | pass | fail | error | skip | warn |\n|---|---|---|---|---|---|---|\n")
for (p in unique(x$pkg)) {
  y <- x[x$pkg == p, ]
  cat(sprintf("| %s | %d | %d | %d | %d | %d | %d |\n", p, nrow(y),
              sum(y$pass, na.rm = TRUE), sum(y$fail, na.rm = TRUE),
              sum(y$err, na.rm = TRUE), sum(y$skip, na.rm = TRUE),
              sum(y$warn, na.rm = TRUE)))
}
cat(sprintf("| **total** | **%d** | **%d** | **%d** | **%d** | **%d** | **%d** |\n",
            nrow(x), sum(x$pass, na.rm = TRUE), sum(x$fail, na.rm = TRUE),
            sum(x$err, na.rm = TRUE), sum(x$skip, na.rm = TRUE),
            sum(x$warn, na.rm = TRUE)))
bad <- x[is.na(x$pass) | x$fail > 0 | x$err > 0 | x$skip > 0 | x$warn > 0, ]
cat("\nfiles with a failure, error, skip, escaped warning or no RESULT:",
    nrow(bad), "\n")
for (i in seq_len(nrow(bad))) cat("-", bad$file[i], "\n")
