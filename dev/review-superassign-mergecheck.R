# Review of lane superassign: read the pairwise trial merges written by
# dev/review-superassign-trialmerge.sh. For every merged file without a
# conflict: parse it, count `<<-` tokens, and for every function compare
# its free variables with those of the two parents. A free variable in
# the merge that neither parent has is a reference one side's change
# left dangling in the other's code.
root <- "C:/Users/adf44/source/r/rev-superassign-src/merge"
defs <- function(path) {
  ex <- parse(path, keep.source = FALSE)
  out <- list()
  for (e in ex) {
    if (is.call(e) && as.character(e[[1L]])[1L] %in% c("<-", "=") &&
          is.name(e[[2L]])) {
      rhs <- e[[3L]]
      if (is.call(rhs) && identical(rhs[[1L]], as.name("function"))) {
        out[[as.character(e[[2L]])]] <- eval(rhs, baseenv())
      }
    }
  }
  out
}
glob <- function(f) unique(codetools::findGlobals(f, merge = FALSE)$variables)
cat("lib: none (source analysis of trial merges)\n")
n_bad <- 0L
for (lane in list.dirs(root, recursive = FALSE, full.names = FALSE)) {
  fl <- list.files(file.path(root, lane), recursive = TRUE, pattern = "\\.R$",
                   full.names = TRUE)
  for (m in fl) {
    txt <- readLines(m, warn = FALSE)
    if (any(grepl("^<<<<<<< ", txt))) {
      cat(sprintf("%-9s %-40s CONFLICTED, not parsed\n", lane, basename(m)))
      next
    }
    pd <- utils::getParseData(parse(m, keep.source = TRUE))
    n_sup <- sum(pd$token == "LEFT_ASSIGN" & pd$text == "<<-")
    dm <- defs(m)
    ds <- defs(paste0(m, ".sa"))
    dt <- defs(paste0(m, ".", lane))
    bad <- character(0)
    for (nm in names(dm)) {
      g <- glob(dm[[nm]])
      ref <- union(if (!is.null(ds[[nm]])) glob(ds[[nm]]),
                   if (!is.null(dt[[nm]])) glob(dt[[nm]]))
      new <- setdiff(g, ref)
      if (length(new)) bad <- c(bad, paste0(nm, ": ", paste(new, collapse = ",")))
    }
    n_bad <- n_bad + length(bad)
    cat(sprintf("%-9s %-40s <<- %d; functions %d; dangling %d\n", lane,
                sub(".*/merge/[^/]+/", "", m), n_sup, length(dm), length(bad)))
    for (b in bad) cat("    DANGLING", b, "\n")
  }
}
cat("RESULT dangling free variables in merged functions:", n_bad, "\n")
