# After the 09:24 crash: read back every recovery RDS; delete the ones
# that do not read (written part-way), and list, per arm, the seeds that
# have a fit without an error in recov3, recov4 or recov5.
# Output: dev/phase3b-log/rds-check.txt
dirs <- c("dev/phase3b-log/recov3", "dev/phase3b-log/recov4",
          "dev/phase3b-log/recov5")
out <- character(0)
bad <- character(0)
have <- list()
for (dd in dirs) {
  for (f in list.files(dd, "[.]rds$", full.names = TRUE)) {
    r <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(r) || is.null(r$seed) || is.null(r$arm)) {
      bad <- c(bad, f)
      next
    }
    if (is.null(r$error)) have[[r$arm]] <- sort(unique(c(have[[r$arm]], r$seed)))
  }
}
out <- c(out, sprintf("unreadable: %d", length(bad)), bad)
if (length(bad)) file.remove(bad)
target <- list(cleft = 1:80, contfix = 1:80, collapse = 1:60, contdl = 1:80,
               cont = 1:60, cens = 1:60)
for (arm in names(target)) {
  miss <- setdiff(target[[arm]], have[[arm]])
  out <- c(out, sprintf("%-8s fitted %3d of %3d; missing: %s", arm,
                        length(intersect(target[[arm]], have[[arm]])),
                        length(target[[arm]]), paste(miss, collapse = " ")))
}
writeLines(out, "dev/phase3b-log/rds-check.txt")
cat(out, sep = "\n")
