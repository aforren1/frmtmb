# Lane predraw2: compare dev/predraw2-measure.R's output of the lane and
# the base build. Per group: fits, fits where the check fired, fits whose
# draws are bitwise those of the base build (both fingerprints identical),
# fits that changed without firing (must be 0), conditions said, and the
# time of the first and second predict() (200 draws) on each build.
# Usage: Rscript dev/predraw2-meas-cmp.R <group>... (reads, with
#        PREDRAW2_MEAS=<prefix> the lane arm from dev/predraw2-<prefix>-*)
#        dev/predraw2-meas-<group>-<lane|base>.rds)
for (g in commandArgs(TRUE)) {
  a <- readRDS(sprintf("dev/rel083-log/kt/predraw2/predraw2-%s-%s-lane.rds",
                       Sys.getenv("PREDRAW2_MEAS", "meas"), g))
  b <- readRDS(sprintf("dev/rel083-log/kt/predraw2/predraw2-meas-%s-base.rds", g))
  ids <- intersect(names(a), names(b))
  fire <- vapply(ids, function(i) isTRUE(a[[i]]$fire), TRUE)
  same <- vapply(ids, function(i) identical(a[[i]]$fp, b[[i]]$fp), TRUE)
  t1a <- vapply(ids, function(i) a[[i]]$t_first, 0)
  t1b <- vapply(ids, function(i) b[[i]]$t_first, 0)
  t2a <- vapply(ids, function(i) a[[i]]$t_second, 0)
  ex <- vapply(ids, function(i) a[[i]]$excess %||% NA_real_, 0)
  cat(sprintf(paste("%-7s fits %3d fired %3d identical %3d changed-unfired",
                    "%d | excess max unfired %.3f | first predict lane %.2f",
                    "s base %.2f s (sum), second lane %.2f s\n"), g,
              length(ids), sum(fire), sum(same), sum(!same & !fire),
              max(ex[!fire], na.rm = TRUE), sum(t1a), sum(t1b), sum(t2a)))
  for (i in ids[fire | !same]) {
    cat(sprintf("   %-36s fire %-5s same %-5s excess %8.3f %s\n", i, fire[i],
                same[i], ex[i], paste(substr(a[[i]]$conds, 1, 70),
                                      collapse = " | ")))
  }
}
