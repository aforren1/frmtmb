# Lane illse, item c: summarize the convergence-warning probe of a suite
# scan (dev/illse-run1scan.R, dev/illse-run1scan2.R): per warning kind,
# how many fits are at the maximum by a restart of nlminb from the stop
# (gain at most 1e-4 in log-likelihood, a tenth of grad_tol) and how
# many are short, and, where the build has conv_headroom(), how its
# verdict (at most grad_tol) agrees with the restart.
#   Rscript dev/illse-convsum.R <probe dir>
d <- commandArgs(TRUE)[1]
fs <- list.files(d, "[.]tsv$", full.names = TRUE)
rows <- do.call(rbind, lapply(fs, function(f) {
  l <- readLines(f, warn = FALSE)
  l <- l[!grepl("^HOOKERROR", l)]
  if (!length(l)) return(NULL)
  x <- strsplit(l, "\t")
  data.frame(pkg = sub("--.*", "", basename(f)),
             file = vapply(x, `[`, "", 1), kinds = vapply(x, `[`, "", 2),
             code = as.integer(vapply(x, `[`, "", 3)),
             np = as.integer(vapply(x, `[`, "", 4)),
             nre = as.integer(vapply(x, `[`, "", 5)),
             head = suppressWarnings(as.numeric(vapply(x, `[`, "", 8))),
             gain = suppressWarnings(as.numeric(vapply(x, `[`, "", 9))),
             ch = suppressWarnings(as.numeric(vapply(x, function(v) {
               if (length(v) >= 13) v[10] else NA_character_
             }, ""))),
             stringsAsFactors = FALSE)
}))
herr <- sum(vapply(fs, function(f) sum(grepl("^HOOKERROR",
                                              readLines(f, warn = FALSE))),
                   0L))
cat("probe rows", nrow(rows), "from", length(fs), "files; probe errors",
    herr, "\n")
rows$kind <- ifelse(grepl("code", rows$kinds), "code != 0",
                    ifelse(grepl("grad", rows$kinds), "large gradient",
                           ifelse(rows$kinds == "", "none of these",
                                  rows$kinds)))
at <- is.finite(rows$gain) & rows$gain <= 1e-4
short <- is.finite(rows$gain) & rows$gain > 1e-4
cat(sprintf("%-16s %6s %8s %8s %8s\n", "kind", "fits", "at max", "short",
            "no restart"))
for (k in unique(rows$kind)) {
  i <- rows$kind == k
  cat(sprintf("%-16s %6d %8d %8d %8d\n", k, sum(i), sum(at & i),
              sum(short & i), sum(!is.finite(rows$gain) & i)))
}
if (any(is.finite(rows$ch))) {
  i <- rows$kind == "code != 0"
  okc <- is.finite(rows$ch) & rows$ch <= 1e-3
  cat("\ncode != 0 fits: conv_headroom() within grad_tol", sum(okc & i),
      "| above it", sum(is.finite(rows$ch) & rows$ch > 1e-3 & i),
      "| not measured", sum(!is.finite(rows$ch) & i), "\n")
  cat("against the restart (at max by restart / short by restart):\n")
  print(table(headroom = ifelse(!is.finite(rows$ch[i]), "NA",
                                ifelse(okc[i], "<= tol", "> tol")),
              restart = ifelse(at[i], "at max", ifelse(short[i], "short",
                                                       "none"))))
}
cat("\nby file (fits, at max by restart):\n")
tb <- aggregate(cbind(fits = 1, at_max = as.integer(at)) ~ pkg + file,
                data = rows, FUN = sum)
print(tb[order(-tb$fits), ], row.names = FALSE)
