# Reviewer (lane ceplot): counts generated from the RESULT lines of a
# dev/ceplot-rev-par.sh log, per arm and package, with every file that
# failed, errored, warned or aborted.
#   Rscript dev/ceplot-rev-counts.R <tag> [<tag> ...]
lg <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-log/"
for (tag in commandArgs(TRUE)) {
  x <- readLines(paste0(lg, tag, ".log"), warn = FALSE)
  res <- grep("^RESULT ", x, value = TRUE)
  noresult <- grep("^NO RESULT", x, value = TRUE)
  loaderr <- grep("LOADERROR", res, value = TRUE)
  res <- setdiff(res, loaderr)
  num <- function(k) as.integer(sub(paste0(".*", k, "=([0-9]+).*"), "\\1", res))
  arm <- sub("^RESULT (\\S+) .*", "\\1", res)
  pkg <- sub("^RESULT \\S+ (\\S+) .*", "\\1", res)
  df <- data.frame(arm, pkg, tests = num("tests"), failed = num("failed"),
                   error = num("error"), skipped = num("skipped"),
                   warning = num("warning"), passed = num("passed"))
  cat("## ", tag, ": ", length(res) + length(loaderr), " RESULT lines, ",
      length(noresult), " without one; ", tail(x, 1), "\n", sep = "")
  for (k in unique(paste(df$arm, df$pkg))) {
    s <- df[paste(df$arm, df$pkg) == k, ]
    cat(sprintf("%-22s files %3d expectations %5d failed %d error %d skipped %d warning %d passed %d\n",
                k, nrow(s), sum(s$tests), sum(s$failed), sum(s$error),
                sum(s$skipped), sum(s$warning), sum(s$passed)))
  }
  bad <- res[df$failed > 0 | df$error > 0 | df$warning > 0]
  if (length(bad) || length(noresult) || length(loaderr)) cat("not clean:\n")
  for (b in c(bad, loaderr, noresult)) cat("  ", b, "\n")
  fails <- grep("^  (fails|warns):", x, value = TRUE)
  for (f in fails) cat("  ", f, "\n")
  sk <- res[df$skipped > 0]
  cat("files with skips:", length(sk), "\n")
}
