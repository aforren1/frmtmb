# Tier counts for dev/predfix-findings.md, generated from the logs.
#   Rscript dev/predfix-tiers.R > dev/predfix-log/tiers.txt
d <- "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/"
row <- function(lbl, f) {
  x <- readLines(paste0(d, f), warn = FALSE)
  r <- grep("^RESULT ", x, value = TRUE)
  num <- function(k) sum(as.integer(sub(paste0(".* ", k, "=([0-9]+).*"),
                                        "\\1", grep(paste0(" ", k, "="), r,
                                                    value = TRUE))))
  bad <- grep(" (fail|err)=[1-9]", r, value = TRUE)
  cat(sprintf("| %s | %d | %d | %d | %d | %d | %d | %d | %s |\n", lbl,
              length(r), num("pass"), num("fail"), num("err"), num("skip"),
              length(grep("^NO RESULT LINE", x)),
              length(grep("LOADERROR", r)),
              format(file.mtime(paste0(d, f)), "%Y-%m-%d %H:%M:%S")))
  if (length(bad)) cat("  files with a failure or error:", bad, sep = "\n  ")
  tail(grep("^(SUITE|GATED) ran", x, value = TRUE), 1)
}
cat("| tier | files | pass | fail | error | skip | no RESULT line |",
    "load error | log mtime |\n|---|---|---|---|---|---|---|---|---|\n")
a <- row("ungated, 1 process per file", "suite.log")
b <- row("gated, all env vars set", "gated.log")
cat("\n", a, "\n", b, "\n")
