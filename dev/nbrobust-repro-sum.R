# Lane nbrobust: one line per log of dev/nbrobust-repro.R, counted from
# its RESULT lines.
#   Rscript dev/nbrobust-repro-sum.R [dir]
dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else
  "dev/nbrobust-out"
fs <- list.files(dir, pattern = "^(base|lib)-.*-(warning|ignore)[.]txt$",
                 full.names = TRUE)
cat("| log | lib | fits | code 0 | code 1 | limit msg | grad warn |",
    "any warn | log shapes (min, max) |\n|---|---|---|---|---|---|---|---|---|\n")
for (f in fs) {
  x <- readLines(f)
  lib <- sub("^lib: (\\S+).*", "\\1", grep("^lib:", x, value = TRUE))
  r <- grep("^RESULT", x, value = TRUE)
  stopifnot(any(grepl("^END", x)))
  num <- function(key) {
    as.numeric(sub(paste0(".* ", key, " (\\S+).*"), "\\1", r))
  }
  flag <- function(key) sub(paste0(".* ", key, " (\\S+).*"), "\\1", r)
  code <- num("code")
  sh <- abs(num("shape"))
  cat(sprintf("| %s | %s | %d | %d | %d | %d | %d | %d | %.4f, %.4f |\n",
              basename(f), basename(lib), length(r), sum(code == 0),
              sum(code == 1), sum(flag("lim") == "TRUE"),
              sum(flag("gradwarn") == "TRUE"), sum(num("nwarn") > 0),
              min(sh), max(sh)))
}
