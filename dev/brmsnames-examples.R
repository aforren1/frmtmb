## Every example of one installed package, donttest included, one topic
## at a time, reporting which topics error.
##   Rscript dev/brmsnames-examples.R frmtmb.sample
pkg <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/brmsnames-libs.R"); brmsnames_libs("lane")
suppressMessages(library(pkg, character.only = TRUE))
db <- tools::Rd_db(pkg)
bad <- character(0)
for (nm in names(db)) {
  f <- tempfile(fileext = ".R")
  tools::Rd2ex(db[[nm]], f, commentDonttest = FALSE)
  if (!file.exists(f)) next
  r <- tryCatch({
    suppressWarnings(suppressMessages(utils::capture.output(
      source(f, local = new.env(), echo = FALSE))))
    "ok"
  }, error = function(e) conditionMessage(e))
  if (!identical(r, "ok")) bad <- c(bad, paste(nm, ":", substr(r, 1, 120)))
}
cat("topics:", length(db), " failing:", length(bad), "\n")
cat(bad, sep = "\n")
cat("DONE\n")
