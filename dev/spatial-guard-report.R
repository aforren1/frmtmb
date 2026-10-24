# Guardrail-1 classification of every brms 2.23.0 name, on the library
# in SPATIAL_LIBS (core alone in the child), with the probes and the
# classifier of tests/testthat/test-guardrail1.R. On a build without
# frmtmb:::brms_vocabulary (rellib-r11), the vocabulary is read from
# this worktree's R/brms-vocabulary.R.
libs <- strsplit(Sys.getenv("SPATIAL_LIBS"), ";", fixed = TRUE)[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
if (!exists("brms_vocabulary", ns)) {
  e <- new.env()
  src <- readLines("R/brms-vocabulary.R")
  a <- grep("^brms_vocabulary <- list", src)
  b <- grep("^# Formula specials an extension supplies", src) - 1L
  eval(parse(text = src[a:b]), envir = e)
  assignInNamespace <- NULL
  voc <- e$brms_vocabulary
} else {
  voc <- ns$brms_vocabulary
}
src <- readLines("tests/testthat/test-guardrail1.R")
first <- grep("^test_that", src)[1]
code <- paste(src[seq_len(first - 1)], collapse = "\n")
code <- gsub("frmtmb:::brms_vocabulary", "voc", code, fixed = TRUE)
eval(parse(text = code))
pr <- guard_probes()
res <- guard_run(pr, strsplit(Sys.getenv("SPATIAL_PKGS"), ",")[[1]])
cls <- mapply(guard_classify, res, pr$name)
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
for (i in seq_len(nrow(pr))) {
  cat(sprintf("%-14s %-28s %-16s %s\n", pr$kind[i], pr$name[i], cls[i],
              substr(res[i], 1, 110)))
}
print(table(pr$kind, cls))
