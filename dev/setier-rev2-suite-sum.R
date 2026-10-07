# Lane setier: compare suite runs of dev/setier-suite.sh.
#   Rscript dev/setier-suite-sum.R <run> [<run> ...]
# Per run: totals per package, then (against the first run) every file
# whose counts differ, and the firings of each condition the lane cares
# about, counted from the .cond logs (raised, whether or not a test
# muffled it).
runs <- commandArgs(TRUE)
root <- "dev/setier-rev2-log"
pat <- c(boundary = "^M\tBoundary [(]singular[)] fit",
         se_lost = "^W\tStandard errors are not available",
         separation = "^W\tThe data separate the outcomes",
         pred_lost = "^W\t[0-9]+ of [0-9]+ .* move along a direction",
         nonfinite_cov = "^W\tSome standard errors are not finite",
         not_conv = "^W\tOptimizer did not report convergence",
         uninformative = "^W\tUninformative interval")
read_run <- function(run) {
  s <- readLines(file.path(root, paste0(run, ".sum")))
  s <- grep(" RESULT ", s, value = TRUE)
  m <- regmatches(s, regexec(paste0("^(\\S+) RESULT (\\S+) pass=(\\d+) ",
                                    "fail=(\\d+) err=(\\d+) skip=(\\d+) ",
                                    "warn=(\\d+)"), s))
  df <- do.call(rbind, lapply(m, function(x) {
    data.frame(pkg = x[2], file = x[3], pass = as.integer(x[4]),
               fail = as.integer(x[5]), err = as.integer(x[6]),
               skip = as.integer(x[7]), warn = as.integer(x[8]))
  }))
  cf <- file.path(root, run, paste0(df$pkg, "--", sub("[.]R$", "", df$file),
                                    ".cond"))
  for (k in names(pat)) {
    df[[k]] <- vapply(cf, function(f) {
      if (!file.exists(f)) return(0L)
      sum(grepl(pat[[k]], readLines(f, warn = FALSE)))
    }, 0L)
  }
  df
}
all <- lapply(runs, read_run)
names(all) <- runs
for (r in runs) {
  df <- all[[r]]
  cat("==", r, ": files", nrow(df), "\n")
  agg <- aggregate(df[, c("pass", "fail", "err", "skip", "warn")],
                   list(pkg = df$pkg), sum)
  agg$files <- as.integer(table(df$pkg)[agg$pkg])
  print(agg, row.names = FALSE)
  cat("totals:", paste(names(colSums(df[, -(1:2)])),
                       colSums(df[, -(1:2)]), collapse = " "), "\n")
  cat("files with a firing:",
      paste(names(pat), vapply(names(pat), function(k) sum(df[[k]] > 0), 0L),
            collapse = " "), "\n")
}
if (length(runs) > 1) {
  a <- all[[1]]
  for (r in runs[-1]) {
    b <- all[[r]]
    k <- merge(a, b, by = c("pkg", "file"), suffixes = c(".a", ".b"))
    cols <- c("pass", "fail", "err", "skip", "warn", names(pat))
    diff <- Reduce(`|`, lapply(cols, function(c) {
      k[[paste0(c, ".a")]] != k[[paste0(c, ".b")]]
    }))
    cat("\n== files differing,", runs[1], "->", r, ":", sum(diff), "\n")
    if (any(diff)) {
      show <- k[diff, ]
      for (i in seq_len(nrow(show))) {
        x <- show[i, ]
        d <- vapply(cols, function(c) {
          va <- x[[paste0(c, ".a")]]; vb <- x[[paste0(c, ".b")]]
          if (va != vb) paste0(c, " ", va, "->", vb) else ""
        }, "")
        cat(sprintf("%-16s %-40s %s\n", x$pkg, x$file,
                    paste(d[nzchar(d)], collapse = ", ")))
      }
    }
  }
}
