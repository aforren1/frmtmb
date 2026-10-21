# Regenerate the bin-2 port into dev/cdfs-out/bin2/ instead of the
# tracked dev/brmsport-bin2/ (the consolidation regenerates that once),
# so the distributions and log_lik blocks can be rerun with this lane's
# shims (dev/brmsport-bin2-adapters.R, -spec.R).
#
#   Rscript dev/cdfs-bin2-gen.R
src <- readLines("dev/brmsport-bin2-gen.R")
hit <- grep('^out_root <- "dev/brmsport-bin2"$', src)
stopifnot(length(hit) == 1L)
src[hit] <- 'out_root <- "dev/cdfs-out/bin2"'
eval(parse(text = src), envir = globalenv())
cat("generated under dev/cdfs-out/bin2:\n")
print(list.files("dev/cdfs-out/bin2", recursive = TRUE))
