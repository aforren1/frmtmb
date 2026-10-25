# Compare the rds files dev/bhaz-bitwise.R wrote: every lane run against
# every base run, element by element, with identical().
#   Rscript dev/bhaz-bitwise-cmp.R <dir>
d <- commandArgs(trailingOnly = TRUE)[[1L]]
lane <- list.files(d, "^lane-.*[.]rds$", full.names = TRUE)
base <- list.files(d, "^base-.*[.]rds$", full.names = TRUE)
cat("lane runs", length(lane), " base runs", length(base), "\n")
for (l in lane) for (b in base) {
  L <- readRDS(l)
  B <- readRDS(b)
  stopifnot(identical(names(L), names(B)))
  for (nm in names(L)) {
    same <- vapply(names(B[[nm]]), function(k) identical(L[[nm]][[k]],
                                                         B[[nm]][[k]]), NA)
    cat(sprintf("%s vs %s  %-9s %s\n", basename(l), basename(b), nm,
                paste(names(same), ifelse(same, "identical", "DIFFER"),
                      collapse = " ")))
  }
}
