# dev/equiv-dump-compat.R on a ";"-separated library list: that script
# takes one library, and the lane build needs its own library in front of
# the base build that holds the extensions.
#   Rscript dev/compat-dump.R <lib;lib> <outdir>
args <- commandArgs(trailingOnly = TRUE)
libs <- strsplit(args[1], ";", fixed = TRUE)[[1]]
src <- readLines("dev/equiv-dump-compat.R")
src <- sub('^\\.libPaths\\(c\\(rellib, ', '.libPaths(c(libs, ', src)
stopifnot(any(grepl("libPaths(c(libs, ", src, fixed = TRUE)))
commandArgs <- function(trailingOnly = FALSE) c(libs[1], args[2])
eval(parse(text = src))
