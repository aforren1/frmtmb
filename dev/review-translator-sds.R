# Reviewer scratch (review of lane translator, 2026-11-09): the smooth
# SDs at frm()'s estimate for each pair the lane lifted to L1 through a
# smooth rule. A block whose SD sits near 0 adds almost nothing to the
# log density at the three points, so L1 cannot see its map or its model.
# Records nothing.
.libPaths(c("C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
e <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = e)
}
smp <- read.delim("dev/equiv-data/sample.tsv")
for (id in c("p042", "p047", "p143", "p048", "p044")) {
  row <- smp[smp$id == id, ]
  key <- paste(row$feature_a, row$feature_b, sep = " x ")
  fit <- suppressWarnings(e$eq_frm_fit(e$eq_specs[[key]]()))
  th <- fit$estimates[["theta"]]
  bl <- vapply(fit$frame$re_blocks, function(b) {
    paste(b$label %||% "", "dim", b$dim %||% NA)
  }, "")
  cat(sprintf("%s %-22s theta (log sd) %s -> sd %s | blocks: %s\n", id, key,
              paste(signif(th, 3), collapse = ","),
              paste(signif(exp(th), 3), collapse = ","),
              paste(bl, collapse = "; ")))
}
