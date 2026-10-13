# Release 0.70.0: lane optrap's studies on rellib-r8 after its
# consolidation fixes (dev/rel070-optrap-jobs-c3.txt) against the same
# studies on the first merged build (dev/rel070-optrap-jobs.txt), every
# column but the timing ones.
#   Rscript dev/rel070-c3-compare.R > dev/rel070-log/studies/optrap/c3/compare.txt
S <- "dev/rel070-log/studies/optrap"
rd <- function(f) utils::read.delim(f, colClasses = "character")
cmp <- function(lab, a, b, drop = c("sec", "secs", "time", "elapsed")) {
  x <- do.call(rbind, lapply(a, rd)); y <- do.call(rbind, lapply(b, rd))
  keep <- setdiff(intersect(names(x), names(y)), drop)
  same <- vapply(seq_len(nrow(x)), function(i)
    identical(unname(unlist(x[i, keep])), unname(unlist(y[i, keep]))), NA)
  cat(sprintf("%-8s rows %d / %d | identical in every non-timing column: %d\n",
              lab, nrow(x), nrow(y), sum(same)))
  if (any(!same)) {
    for (k in keep) {
      d <- which(x[[k]] != y[[k]])
      if (length(d)) cat("   column", k, "differs on", length(d), "rows\n")
    }
  }
}
mo1 <- sort(list.files(file.path(S, "mo"), "^r8-.*[.]tsv$", full.names = TRUE))
mo1 <- mo1[order(as.integer(vapply(strsplit(basename(mo1), "-"), `[`, "", 2)))]
mo2 <- file.path(S, "c3", paste0("mo-", c("1-100", "101-200", "201-300",
                                          "301-400"), ".tsv"))
cmp("mo()", mo1, mo2)
cmp("scale", file.path(S, "optrap-log/scale-r8.tsv"),
    file.path(S, "optrap-log/scale-r8c3.tsv"))
cmp("wide", file.path(S, "rev3-log/b1-wide-r8.tsv"),
    file.path(S, "rev3-log/b1-wide-r8c3.tsv"))
h1 <- readLines(file.path(S, "healthy.txt")); h2 <- readLines(file.path(S, "c3/healthy.txt"))
cat("healthy  lines", length(h1), "/", length(h2), "| identical:",
    identical(h1[-1], h2[-1]), "(after the first line)\n")
