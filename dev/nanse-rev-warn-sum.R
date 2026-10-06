# Reviewer: per fit (top-level frm() index within a test file), which
# kinds of warning frm_warning() raised, from the .warn logs of
# dev/nanse-rev-suite.sh. Depth 0 means raised outside frm() (a later
# vcov()/summary() on the most recent fit, attributed to it).
#   Rscript dev/nanse-rev-warn-sum.R <suite dir>
d <- commandArgs(trailingOnly = TRUE)[1]
fs <- list.files(d, "[.]warn$", full.names = TRUE)
kind <- function(m) {
  if (grepl("^Standard errors are not available", m)) return("SE")
  if (grepl("^Some standard errors are not finite", m)) return("VCOV")
  if (grepl("^Optimizer did not report|^The gradient at the reported|gradient", m) &&
      grepl("^Optimizer|^The gradient|^Convergence|^Maximum absolute gradient|grad", m)) {
    return("CONV")
  }
  if (grepl("are not identified: at the optimum", m)) return("NLFLAT")
  if (grepl("^Hessian is not positive definite", m)) return("PDHESS")
  if (grepl("has a single level|gives every observation its own", m)) {
    return("STRUCT")
  }
  if (grepl("^mixture[(]|^thres[(]|zero_one_inflated|xbeta|disc", m)) {
    return("FAMILY")
  }
  "OTHER"
}
rows <- list()
for (f in fs) {
  L <- readLines(f, warn = FALSE)
  if (!length(L)) next
  p <- strsplit(L, "\t", fixed = TRUE)
  rows[[f]] <- data.frame(file = sub("[.]warn$", "", basename(f)),
                          fit = vapply(p, `[`, "", 1),
                          depth = vapply(p, `[`, "", 2),
                          kind = vapply(p, function(x) kind(x[3]), ""),
                          msg = vapply(p, function(x) substr(x[3], 1, 110), ""))
}
w <- do.call(rbind, rows)
cat("warning lines:", nrow(w), " files:", length(unique(w$file)), "\n")
print(table(w$kind, w$depth))
key <- paste(w$file, w$fit)
by_fit <- split(w, key)
se_fits <- Filter(function(x) any(x$kind == "SE"), by_fit)
cat("\nfits with the SE warning:", length(se_fits), "\n")
two <- Filter(function(x) sum(x$kind %in% c("SE", "VCOV")) > 1 ||
                (any(x$kind %in% c("SE", "VCOV")) &&
                   any(x$kind %in% c("CONV", "NLFLAT", "STRUCT", "FAMILY",
                                     "PDHESS"))), by_fit)
cat("fits with an SE/VCOV warning AND another SE-explaining warning:",
    length(two), "\n")
for (x in two) {
  cat("--", x$file[1], "fit", x$fit[1], "\n")
  for (i in seq_len(nrow(x))) {
    cat("   ", x$kind[i], "d", x$depth[i], ":", x$msg[i], "\n")
  }
}
cat("\nSE warnings per file:\n")
print(sort(table(w$file[w$kind == "SE"]), decreasing = TRUE))
cat("\nVCOV warnings per file:\n")
print(sort(table(w$file[w$kind == "VCOV"]), decreasing = TRUE))
