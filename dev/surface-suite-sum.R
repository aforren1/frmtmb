# Lane surface: totals of one suite run, per package, from the RESULT
# lines of dev/surface-suite-<tag>/, with every file that failed,
# errored, warned, skipped or wrote no RESULT line, and the lib line
# check.
#
#   Rscript dev/surface-suite-sum.R <tag> [rerun tag]
#
# With a rerun tag, a file present in the rerun replaces its first
# result, and both results are listed.
a <- commandArgs(TRUE)
dir1 <- file.path("dev", paste0("surface-suite-", a[1]))
dir2 <- if (length(a) > 1) file.path("dev", paste0("surface-suite-", a[2]))
read_dir <- function(d) {
  fs <- list.files(d, pattern = "[.]txt$", full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    x <- readLines(f, warn = FALSE)
    r <- grep("^RESULT ", x, value = TRUE)
    lib <- grep("^lib: ", x, value = TRUE)
    lane_lib <- length(lib) &&
      grepl("wt-surface-lib/frmtmb ", paste0(lib, " ")) &&
      grepl("wt-surface-lib/frmtmb$|wt-surface-lib/frmtmb ",
            sub(".*[|] frmtmb ", "", lib))
    nm <- sub("[.]txt$", "", basename(f))
    pk <- sub("--.*$", "", nm)
    if (!length(r)) {
      return(data.frame(pkg = pk, file = nm, pass = NA, fail = NA,
                        error = NA, skip = NA, warn = NA,
                        lanecore = lane_lib))
    }
    v <- as.integer(sub(".*=", "", strsplit(r[1], " ")[[1]][4:8]))
    data.frame(pkg = pk, file = nm, pass = v[1], fail = v[2], error = v[3],
               skip = v[4], warn = v[5], lanecore = lane_lib)
  }))
}
t1 <- read_dir(dir1)
if (!is.null(dir2)) {
  t2 <- read_dir(dir2)
  cat("Files run again in", dir2, ":\n")
  for (i in seq_len(nrow(t2))) {
    j <- match(t2$file[i], t1$file)
    cat(sprintf("- %s: first pass %s fail %s err %s warn %s, rerun pass %s fail %s err %s warn %s\n",
                t2$file[i], t1$pass[j], t1$fail[j], t1$error[j], t1$warn[j],
                t2$pass[i], t2$fail[i], t2$error[i], t2$warn[i]))
    t1[j, ] <- t2[i, ]
  }
  cat("\n")
}
cat("files:", nrow(t1), "| with a RESULT line:", sum(!is.na(t1$pass)),
    "| whose frmtmb is the lane's:", sum(t1$lanecore), "\n\n")
agg <- aggregate(cbind(files = 1, pass, fail, error, skip, warn) ~ pkg,
                 data = transform(t1, files = 1), FUN = sum,
                 na.action = na.pass)
ord <- c("frmtmb", "frmtmb.sample", "frmtmb.coupling", "frmtmb.eam",
         "frmtmb.latent", "frmtmb.learn", "frmtmb.ode", "frmtmb.spline")
agg <- agg[match(ord, agg$pkg), ]
cat("| package | files | pass | fail | error | skip | warn |\n")
cat("|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(agg))) {
  cat(sprintf("| %s | %d | %d | %d | %d | %d | %d |\n", agg$pkg[i],
              agg$files[i], agg$pass[i], agg$fail[i], agg$error[i],
              agg$skip[i], agg$warn[i]))
}
cat(sprintf("| **total** | **%d** | **%d** | **%d** | **%d** | **%d** | **%d** |\n",
            sum(agg$files), sum(agg$pass), sum(agg$fail), sum(agg$error),
            sum(agg$skip), sum(agg$warn)))
bad <- t1[is.na(t1$pass) | t1$fail > 0 | t1$error > 0 | t1$warn > 0, ]
cat("\nfiles with a failure, an error, a warning or no RESULT line:",
    nrow(bad), "\n")
for (i in seq_len(nrow(bad))) {
  cat(sprintf("- %s: pass %s fail %s error %s warn %s\n", bad$file[i],
              bad$pass[i], bad$fail[i], bad$error[i], bad$warn[i]))
}
sk <- t1[!is.na(t1$skip) & t1$skip > 0, ]
cat("\nfiles with a skip:", nrow(sk), "\n")
for (i in seq_len(nrow(sk))) {
  cat(sprintf("- %s: %d\n", sk$file[i], sk$skip[i]))
}
