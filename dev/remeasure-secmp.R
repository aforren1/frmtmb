# The standard errors of every vignette fit, build against build. A call
# that still runs can still regress: brms_missings fit_imp2 kept its
# optimum at 0.70.0 and lost every standard error. This lists every fit
# whose estimates, standard errors or optimizer code moved.
#
#   Rscript dev/remeasure-secmp.R <tag> [<tag> ...]
#
# Reads the spell pass of dev/remeasure-out/<tag>/port/results-spell,
# where run-vignette.R keeps fixef() estimates and standard errors.
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
tags <- commandArgs(trailingOnly = TRUE)
rd <- function(t) {
  fs <- list.files(file.path(root, "remeasure-out", t, "port",
                             "results-spell"), "[.]rds$", full.names = TRUE)
  out <- list()
  for (f in fs) out <- c(out, readRDS(f))
  out[vapply(out, function(r) identical(r$kind, "model") &&
               identical(r$status, "OK"), NA)]
}
P <- lapply(tags, rd)
names(P) <- tags
ids <- sort(unique(unlist(lapply(P, names))))
rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-8)
cat("fits compared:", length(ids), "; builds:", paste(tags, collapse = " "),
    "\n")
n_moved <- 0L
for (id in ids) {
  base <- P[[1]][[id]]
  lines <- character()
  for (t in tags[-1]) {
    r <- P[[t]][[id]]
    if (is.null(base) || is.null(r)) {
      lines <- c(lines, sprintf("  %-6s %s", t, if (is.null(r))
        "does not run" else "runs (not in the first build)"))
      next
    }
    k <- intersect(names(base$fit$fixef_se), names(r$fit$fixef_se))
    se0 <- base$fit$fixef_se[k]; se1 <- r$fit$fixef_se[k]
    b0 <- base$fit$fixef[k]; b1 <- r$fit$fixef[k]
    nan_new <- k[!is.finite(se1) & is.finite(se0)]
    nan_gone <- k[is.finite(se1) & !is.finite(se0)]
    ok <- is.finite(se0) & is.finite(se1)
    dse <- if (any(ok)) max(rel(se1[ok], se0[ok])) else 0
    db <- if (length(k)) max(rel(b1, b0), na.rm = TRUE) else 0
    cv <- c(base$fit$conv %||% NA, r$fit$conv %||% NA)
    if (length(nan_new) || length(nan_gone) || dse > 0.01 || db > 0.01 ||
        !identical(cv[1], cv[2])) {
      lines <- c(lines, sprintf(paste(
        "  %-6s code %s -> %s; max rel change: estimate %.3g, SE %.3g;",
        "SE lost: %s; SE gained: %s"), t, cv[1], cv[2], db, dse,
        if (length(nan_new)) paste(nan_new, collapse = " ") else "-",
        if (length(nan_gone)) paste(nan_gone, collapse = " ") else "-"))
    }
  }
  if (length(lines)) {
    n_moved <- n_moved + 1L
    cat("\n", id, " (against ", tags[1], ")\n", sep = "")
    cat(lines, sep = "\n")
    cat("\n")
  }
}
cat("\nfits with a change:", n_moved, "\n")
