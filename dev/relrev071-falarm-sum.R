# Release reviewer of 0.71.0: per arm, over the false-alarm battery runs
# of dev/relrev071-falarm.sh, the silent fits (code 0, no warning) whose
# coefficient SEs are off their reference by more than 1 percent, and
# the max_rel moves between the merged build and each lane alone.
#   Rscript dev/relrev071-falarm-sum.R <dir> [<dir> ...]
dirs <- commandArgs(trailingOnly = TRUE)
rd <- function(f) utils::read.delim(f, colClasses = "character")
arms <- c("merged", "base", "optstart", "illse")
for (cfg in c("ref", "ob")) {
  X <- list()
  for (arm in arms) {
    X[[arm]] <- do.call(rbind, lapply(dirs, function(d) do.call(rbind,
      lapply(c("nl", "lin", "poly", "mixed", "bounded"), function(g)
        rd(file.path(d, "falarm", sprintf("%s-%s-%s.tsv", g, cfg, arm)))))))
  }
  cat("\n==== BLAS", cfg, "fits per arm", nrow(X$merged), "\n")
  for (arm in arms) {
    x <- X[[arm]]
    r <- suppressWarnings(as.numeric(x$max_rel))
    silent <- x$code == "0" & x$se_warn == "0" & x$other_warn == "0"
    bad <- silent & is.finite(r) & r > 0.01
    cat(sprintf("%-9s silent %4d; silent with max_rel > 1%%: %3d; > 10%%: %3d;",
                arm, sum(silent), sum(bad),
                sum(silent & is.finite(r) & r > 0.1)),
        " models:", paste(names(table(x$model[bad])), table(x$model[bad]),
                          collapse = ", "), "\n")
  }
  m <- suppressWarnings(as.numeric(X$merged$max_rel))
  for (arm in c("optstart", "illse")) {
    o <- suppressWarnings(as.numeric(X[[arm]]$max_rel))
    ok <- is.finite(m) & is.finite(o) & m != o
    cat(sprintf(paste0("merged vs %-8s max_rel differs %d: merged smaller",
                       " %d, larger %d; largest ratio merged/%s %.3g (%s",
                       " seed %s: %s vs %s)\n"),
                arm, sum(ok), sum(m[ok] < o[ok]), sum(m[ok] > o[ok]), arm,
                max(m[ok] / o[ok]), X$merged$model[ok][which.max(m[ok] / o[ok])],
                X$merged$seed[ok][which.max(m[ok] / o[ok])],
                m[ok][which.max(m[ok] / o[ok])],
                o[ok][which.max(m[ok] / o[ok])]))
    big <- ok & m > 0.01 & m > 2 * o
    if (any(big)) {
      print(data.frame(model = X$merged$model[big], seed = X$merged$seed[big],
                       merged = m[big], other = o[big]))
    }
  }
}
