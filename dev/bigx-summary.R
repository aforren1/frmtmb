# Lane bigx: summarize a grid TSV (dev/bigx-grid.R) per family, group
# term and offset. A cell is "right" when the slope is there, within
# 0.01 SE of the centered fit, the log-likelihood within 1e-4 of it, and
# the SE ratio within 1e-3 of 1. "silent wrong" is a cell that is not
# right and raised no warning.
#   Rscript dev/bigx-summary.R <grid.tsv> [<grid.tsv> ...]
for (f in commandArgs(TRUE)) {
  t <- read.delim(f, stringsAsFactors = FALSE)
  t$dropped <- !t$w_err & !t$w_haslope %in% TRUE
  t$right <- !t$w_err & t$w_haslope %in% TRUE & abs(t$slope_z) < 0.01 &
    abs(t$dll) < 1e-4 & abs(t$se_ratio - 1) < 1e-3
  t$right[is.na(t$right)] <- FALSE
  t$silent_wrong <- !t$right & t$w_warn == 0 & !t$w_err
  cat("==", f, "cells", nrow(t), "\n")
  key <- interaction(t$family, t$re, t$m, drop = TRUE, lex.order = TRUE)
  s <- do.call(rbind, lapply(split(t, key), function(d) {
    data.frame(family = d$family[1], re = d$re[1], m = d$m[1],
               n = nrow(d), right = sum(d$right), dropped = sum(d$dropped),
               warned = sum(d$w_warn > 0), silent_wrong = sum(d$silent_wrong),
               max_abs_slope_z = suppressWarnings(max(abs(d$slope_z),
                                                     na.rm = TRUE)),
               min_dll = suppressWarnings(min(d$dll, na.rm = TRUE)),
               max_se_dev = suppressWarnings(max(abs(d$se_ratio - 1),
                                                 na.rm = TRUE)))
  }))
  s <- s[order(s$family, -s$re, s$m), ]
  print(s, row.names = FALSE, digits = 3)
  cat("TOTAL cells", nrow(t), "right", sum(t$right), "dropped",
      sum(t$dropped), "warned", sum(t$w_warn > 0), "silent wrong",
      sum(t$silent_wrong), "errors", sum(t$w_err), "\n")
  big <- t$m >= 1e7
  cat("m >= 1e7: cells", sum(big), "dropped", sum(t$dropped[big]),
      "dll range of dropped", format(range(t$dll[big & t$dropped]),
                                     digits = 3), "\n")
}
