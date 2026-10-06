# Reviewer: summarize the tier-1 noise-row log of an instrumented suite
# run. Usage: Rscript dev/cifixrev-revsum.R <run name>
args <- commandArgs(TRUE)
dir <- file.path("C:/Users/adf44/source/r/frmtmb-wt-cifix/dev/cifixrev-log",
                 paste0(args[1], "-suite"))
fs <- list.files(dir, pattern = "[.]rev$", full.names = TRUE)
rows <- do.call(rbind, lapply(fs, function(f) {
  x <- readLines(f, warn = FALSE)
  x <- x[nzchar(x)]
  do.call(rbind, lapply(strsplit(x, "\t", fixed = TRUE), function(v) {
    length(v) <- 8
    as.data.frame(t(v))
  }))
}))
names(rows) <- c("tag", "site", "np", "pars", "est", "rowmax", "se", "call")
rows$key <- paste(rows$tag, rows$call, rows$pars, rows$est)
u <- rows[!duplicated(rows$key), ]
cat("log lines:", nrow(rows), " distinct fits (file, call, params, estimate):",
    nrow(u), " test files:", length(unique(u$tag)), "\n")
cat("by package:\n")
print(table(sub("/.*", "", u$tag)))
cat("sites seen per distinct fit:\n")
print(table(tapply(rows$site, rows$key, function(s) {
  paste(sort(unique(s)), collapse = "+")
})))
se <- suppressWarnings(as.numeric(unlist(strsplit(u$se, ","))))
est <- suppressWarnings(as.numeric(unlist(strsplit(u$est, ","))))
rm <- suppressWarnings(as.numeric(unlist(strsplit(u$rowmax, ","))))
cat("noise-row parameters:", length(se), "\n")
cat("their SE (optimizer scale) quantiles:\n")
print(signif(stats::quantile(se, c(0, 0.1, 0.5, 0.9, 1), na.rm = TRUE), 3))
cat("their estimates quantiles:\n")
print(signif(stats::quantile(est, c(0, 0.1, 0.5, 0.9, 1), na.rm = TRUE), 3))
cat("their row max quantiles:\n")
print(signif(stats::quantile(rm, c(0, 0.1, 0.5, 0.9, 1), na.rm = TRUE), 3))
cat("parameter names:\n")
print(table(unlist(strsplit(u$pars, ","))))
cat("\nper test file (distinct fits):\n")
print(sort(table(u$tag), decreasing = TRUE))
cat("\nsmallest SEs reported from a noise row:\n")
o <- order(se)
print(head(data.frame(se = se, est = est, rowmax = rm)[o, ], 10))
utils::write.table(u[, c("tag", "site", "np", "pars", "est", "rowmax", "se",
                         "call")],
                   file.path(dirname(dir), paste0(args[1], "-rev.tsv")),
                   sep = "\t", row.names = FALSE, quote = FALSE)
