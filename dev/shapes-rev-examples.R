# Reviewer, priority 3: RUN every Rd example in all eight packages and
# record which ones error, with a trace on predict.frmtmb_fit that says
# which examples reach the renamed generic at all.
#
#   Rscript dev/shapes-rev-examples.R lane <pkg>
#   Rscript dev/shapes-rev-examples.R base <pkg>
#
# One package per process: an example that leaves state behind cannot
# then reach the next package's.

arg <- commandArgs(trailingOnly = TRUE)
which_lib <- if (identical(arg[1], "base")) "base" else "lane"
pkg <- arg[2]
lib <- if (which_lib == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/shapes-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
Sys.setenv(NOT_CRAN = "true")

suppressPackageStartupMessages(library(pkg, character.only = TRUE))

db <- tools::Rd_db(pkg)
res <- list()
tf <- tempfile(fileext = ".R")
for (nm in names(db)) {
  ex <- tryCatch(utils::capture.output(
    tools::Rd2ex(db[[nm]], out = tf, defines = NULL)), error = function(e) NULL)
  if (!file.exists(tf)) { res[[nm]] <- "no-examples"; next }
  src <- readLines(tf, warn = FALSE)
  unlink(tf)
  if (!length(src) || !any(nzchar(trimws(src)))) {
    res[[nm]] <- "no-examples"; next
  }
  e <- new.env(parent = globalenv())
  r <- tryCatch({
    suppressWarnings(suppressMessages(
      utils::capture.output(eval(parse(text = src), envir = e))))
    "ok"
  }, error = function(c) paste0("ERROR: ", conditionMessage(c)))
  res[[nm]] <- r
}
bad <- res[!vapply(res, function(x) x %in% c("ok", "no-examples"), TRUE)]
cat("== ", pkg, " (", which_lib, ") Rd pages: ", length(res),
    "  with examples: ",
    sum(!vapply(res, identical, TRUE, "no-examples")),
    "  errors: ", length(bad), "\n", sep = "")
for (nm in names(bad)) cat("  [", nm, "] ", bad[[nm]], "\n", sep = "")
saveRDS(res, file.path(TREE, sprintf("dev/shapes-rev-ex-%s-%s.rds",
                                     which_lib, pkg)))
