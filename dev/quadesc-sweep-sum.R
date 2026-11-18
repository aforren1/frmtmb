# Lane quadesc: summarize the sweep logs of dev/quadesc-sweep.sh.
#   Rscript dev/quadesc-sweep-sum.R <tag>
# Per log: fits read, fits stopped on alpha's point (every alpha
# coefficient below 1 in absolute value), the k of those, the range of the
# fit's exact marginal log-likelihood, the smallest margin over the Laplace
# fit's exact value, and the warning classes seen.
tag <- commandArgs(TRUE)[1]
fs <- Sys.glob(sprintf("dev/quadesc-log/sw-%s-*.txt", tag))
for (f in fs) {
  x <- readLines(f)
  fit <- grep("^FIT ", x, value = TRUE)
  lapx <- grep("^LAPX ", x, value = TRUE)
  fld <- function(l, key) {
    vapply(strsplit(l, " +"), function(v) v[which(v == key)[1] + 1], "")
  }
  k <- as.integer(fld(fit, "k"))
  al <- strsplit(fld(fit, "alpha"), ",")
  # the alpha field can be split by format() padding: read the alphas from
  # the text between "alpha" and "code"
  al <- lapply(sub(".* alpha +(.*) code .*", "\\1", fit), function(s) {
    as.numeric(strsplit(trimws(s), "[, ]+")[[1]])
  })
  on <- vapply(al, function(a) all(abs(a) < 1), NA)
  ex <- as.numeric(fld(fit, "exact"))
  mg <- as.numeric(fld(lapx, "margin"))
  wcls <- sub("^ *WARN ([a-z_]+):.*", "\\1",
              grep("^ *WARN ", x, value = TRUE))
  cat(sprintf(paste0("%-40s fits %2d on_point %2d k=[%s] exact %.4f..%.4f",
                     " min_margin %s exit %s\n  warns: %s\n"),
              basename(f), length(fit), sum(on), paste(k[on], collapse = ","),
              min(ex), max(ex),
              if (length(mg)) format(min(mg), digits = 4) else "NA",
              sub("exit ", "", tail(grep("^exit", x, value = TRUE), 1)),
              paste(names(table(wcls)), table(wcls), sep = "=",
                    collapse = " ")))
}
