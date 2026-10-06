# Reviewer: house style of the lane's prose. Counts, per record, lines
# over 80 columns OUTSIDE generated blocks, tables and code fences;
# em and en dashes; non-ASCII symbols outside the generated blocks;
# spaced hyphens; British spellings; control characters anywhere.
#
#   Rscript dev/vigport-rev-style.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
# only the lane's added lines count for the two older records
added <- function(f) {
  d <- system2("git", c("-C", root, "diff", "-U0", "--", f), stdout = TRUE)
  sub("^[+]", "", grep("^[+][^+]", d, value = TRUE))
}
files <- list(
  "vigport-findings.md" = readLines(file.path(root, "vigport-findings.md"),
                                    warn = FALSE),
  "brms-vignette-port.md" = added("brms-vignette-port.md"),
  "brms-vignette-audit.md" = added("brms-vignette-audit.md"))
brit <- paste0("\\b(behaviour|colour|modell|labell|favour|analys(e|ed|ing)",
               "\\b|organis|recognis|summaris|optimis|minimis|realis|",
               "parameteris|centre|licence|programme)")
for (nm in names(files)) {
  x <- files[[nm]]
  gen <- FALSE; fence <- FALSE; prose <- logical(length(x))
  for (i in seq_along(x)) {
    if (grepl("^<!-- BEGIN generated", x[i])) gen <- TRUE
    if (grepl("^```", x[i])) fence <- !fence
    prose[i] <- !gen && !fence && !grepl("^[|]", x[i])
    if (grepl("^<!-- END generated", x[i])) gen <- FALSE
  }
  long <- which(prose & nchar(x, type = "width") > 80)
  cat(sprintf(paste0("%-24s lines %d; prose over 80: %d; em/en dash: %d;",
                     " spaced hyphen in prose: %d; British: %d;",
                     " control chars: %d; non-ASCII in prose: %d\n"),
              nm, length(x), length(long),
              sum(grepl("\u2014|\u2013", x)),
              sum(prose & grepl("[^ ] - [^ ]", x) & !grepl("^ *- ", x)),
              sum(grepl(brit, x, ignore.case = TRUE)),
              sum(grepl("[\001-\010\013\014\016-\037]", x, useBytes = TRUE)),
              sum(prose & grepl("[^\001-\177]", x, useBytes = TRUE))))
  for (i in long) cat("    long:", substr(x[i], 1, 70), "\n")
  for (i in which(grepl(brit, x, ignore.case = TRUE)))
    cat("    British?:", substr(x[i], 1, 90), "\n")
  for (i in which(prose & grepl("[^ ] - [^ ]", x) & !grepl("^ *- ", x)))
    cat("    spaced hyphen:", substr(x[i], 1, 90), "\n")
}
# R sources the lane wrote or changed: lines over 80
rs <- c(list.files(root, "^vigport-[a-z-]*[.]R$", full.names = TRUE),
        file.path(root, "brms-port", c("brms-fit.R", "env.R", "plausibility.R",
                                       "shim.R", "run-vignette.R",
                                       "summarize.R", "run-all.R",
                                       "port-lib.R")),
        file.path(root, "brms-vignettes", c("_drift.R", "_workarounds.R",
                                            "_harness.R")))
rs <- rs[!grepl("vigport-rev-", rs)]
for (f in rs) {
  x <- readLines(f, warn = FALSE)
  n <- sum(nchar(x, type = "width") > 80)
  if (n) cat(sprintf("%-40s R lines over 80: %d\n", basename(f), n))
}
