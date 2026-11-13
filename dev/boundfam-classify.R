# Lane boundfam: every firing of the limit rule in a traced lane suite,
# beside the same analysis on the base build.
#   Rscript dev/boundfam-classify.R <lane dir> <base dir>
# One block per LIM line: the file, the record number, what the rule
# found (name:kind, estimates, directions, gains on the ray), and the
# standard errors that changed (base -> lane) in the COV record that
# follows it. A base record that is missing prints "no base".
a <- commandArgs(trailingOnly = TRUE)
dl <- a[1]
db <- a[2]
num <- function(s) suppressWarnings(as.numeric(strsplit(s, ",", fixed = TRUE)[[1]]))
covs <- function(x) {
  x <- x[startsWith(x, "COV ")]
  k <- as.integer(sub("^COV ([0-9]+) .*", "\\1", x))
  stats::setNames(lapply(x, function(l) {
    list(se = num(sub(".*se=", "", l)),
         lost = sub(" se=.*", "", sub(".*lost=", "", l)))
  }), k)
}
tot <- c(fired = 0L, already_nan = 0L, newly_nan = 0L, no_base = 0L)
for (f in Sys.glob(file.path(dl, "*.trace"))) {
  xl <- readLines(f, warn = FALSE)
  li <- which(startsWith(xl, "LIM "))
  if (!length(li)) next
  bf <- file.path(db, basename(f))
  cb <- if (file.exists(bf)) covs(readLines(bf, warn = FALSE))
  cl <- covs(xl)
  for (i in li) {
    tot["fired"] <- tot["fired"] + 1L
    k <- as.integer(sub("^LIM ([0-9]+) .*", "\\1", xl[i])) + 1L
    cat(basename(f), " #", k, "\n  ", substr(xl[i], 1, 400), "\n", sep = "")
    lane <- cl[[as.character(k)]]
    base <- if (!is.null(cb)) cb[[as.character(k)]]
    if (is.null(lane)) next
    if (is.null(base) || length(base$se) != length(lane$se)) {
      tot["no_base"] <- tot["no_base"] + 1L
      cat("  no base; lane lost", lane$lost, "\n")
      next
    }
    both <- !is.na(lane$se) & !is.na(base$se)
    if (any(both)) {
      r <- lane$se[both] / base$se[both]
      cat("  kept SEs, lane / base: min", signif(min(r), 6), "max",
          signif(max(r), 6), "\n")
    }
    new <- which(is.na(lane$se) & !is.na(base$se))
    if (length(new)) {
      tot["newly_nan"] <- tot["newly_nan"] + 1L
      cat("  base SE -> NaN at positions", new, ":",
          signif(base$se[new], 4), "| lane lost", lane$lost, "\n")
    } else {
      tot["already_nan"] <- tot["already_nan"] + 1L
      cat("  no SE changed; base lost", base$lost, "\n")
    }
  }
}
cat("TOTAL", paste(names(tot), tot, sep = "=", collapse = " "), "\n")
