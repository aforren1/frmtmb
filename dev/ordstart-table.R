# Lane ordstart: summarise dev/ordstart-starts.R output. Per seed: the
# default (run 0), the best of the jittered runs, and the best of those
# the degeneracy rule calls proper; each classified by which part of
# mixture_ord_degeneracy()'s rule fires (step, runoff, collapsed), plus
# weight (a component's mean mixing weight below 1e-3) and big (an outer
# parameter above 30 in absolute value, the 0.68.0 review's label).
#   Rscript dev/ordstart-table.R file.tsv [file.tsv ...]
args <- commandArgs(trailingOnly = TRUE)
X <- do.call(rbind, lapply(args, function(f) {
  utils::read.delim(f, stringsAsFactors = FALSE)
}))
kind <- function(r) {
  k <- c(
    step = any(r[c("sharpen1", "sharpen2")] > -0.1, na.rm = TRUE),
    runoff = any(r[c("near1", "near2")] > 50, na.rm = TRUE),
    collapsed = any(as.logical(r[c("collapsed1", "collapsed2")]),
                    na.rm = TRUE),
    weight = any(pmin(r$weight1, r$weight2) < 1e-3, na.rm = TRUE),
    big = isTRUE(r$maxabs > 30))
  if (!any(k)) "proper" else paste(names(k)[k], collapse = "+")
}
X$kind <- vapply(seq_len(nrow(X)), function(i) kind(X[i, ]), "")
X$rule <- X$kind != "proper" &
  grepl("step|runoff|collapsed", X$kind)
out <- NULL
for (key in unique(paste(X$family, X$seed))) {
  S <- X[paste(X$family, X$seed) == key & !is.na(X$ll), ]
  d0 <- S[S$run == 0, ]
  J <- S[S$run > 0, ]
  b <- J[which.max(J$ll), ]
  P <- J[!J$rule & !J$degen_warn, ]
  # the best proper over the default and the jittered runs
  A <- rbind(if (!d0$rule && !d0$degen_warn) d0, P)
  bp <- if (nrow(A)) A[which.max(A$ll), ] else NULL
  out <- rbind(out, data.frame(
    family = S$family[1], seed = S$seed[1],
    ll0 = d0$ll, code0 = d0$code, warn0 = d0$degen_warn, kind0 = d0$kind,
    ll_best = b$ll, gain = b$ll - d0$ll, warn_best = b$degen_warn,
    kind_best = b$kind,
    ll_proper = if (is.null(bp)) NA else bp$ll,
    gain_proper = if (is.null(bp)) NA else bp$ll - d0$ll,
    n_proper = nrow(P), n_runs = nrow(J)))
}
op <- options(width = 200)
print(out, digits = 7, row.names = FALSE)
cat("\nrule vs warning, every run: ",
    "agree", sum(X$rule == X$degen_warn, na.rm = TRUE),
    "of", sum(!is.na(X$degen_warn)), "\n")
print(table(kind = X$kind, warn = X$degen_warn, useNA = "ifany"))
for (fm in unique(out$family)) {
  o <- out[out$family == fm, ]
  g <- o$gain > 1e-3
  cat(fm, ": seeds", nrow(o), "; best of", o$n_runs[1],
      "above the default by > 1e-3:", sum(g),
      "; of those degenerate (warning):", sum(g & o$warn_best),
      "; non-degenerate:", sum(g & !o$warn_best),
      "; default degenerate (warning):", sum(o$warn0),
      "; best proper above the default by > 1e-3:",
      sum(o$gain_proper > 1e-3, na.rm = TRUE),
      "; default degenerate and a proper point found:",
      sum(o$warn0 & !is.na(o$ll_proper)), "\n")
}
g <- out$gain > 1e-3
cat("ALL: seeds", nrow(out), "; best above default:", sum(g),
    "; of those degenerate:", sum(g & out$warn_best),
    "; kinds of the degenerate better points:\n")
print(table(out$kind_best[g & out$warn_best]))
