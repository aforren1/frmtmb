## Reviewer, claims 4 and 6: does each r_ column hold what its name says,
## on the multi-component models? Reads the saved draws of
## dev/brmsnames-rev-frm.R and checks draws 1, 50 and 100 against the ML
## ranef() of a fit with its estimates replaced by that draw.
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
for (nm in c("C1", "C3", "C4", "C5", "C7")) {
  fo <- readRDS(sprintf("dev/stan-cache/brmsnames-rev-frm-%s.rds", nm))
  ds <- fo$ds; fit <- fo$fit
  idx <- frmtmb.sample:::draws_par_index(fit)
  cat("==", nm, "== blocks:", vapply(fit$frame$re_blocks, `[[`, "", "term_label"), "\n")
  for (bk in fit$frame$re_blocks) {
    cat("  block", bk$term_label, "cnms:", bk$cnms, " components:",
        vapply(bk$components, function(c) paste0(c$lp_key, "[", paste(c$cnms, collapse = ","), "]"), ""), "\n")
  }
  checked <- 0; wrong <- 0; ex <- character(0)
  for (i in c(1, 50, 100)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, i, idx)
    re <- unclass(ranef(sh))
    for (b in seq_along(re)) {
      bk <- fit$frame$re_blocks[[b]]
      M <- re[[b]]
      for (k in seq_len(ncol(M))) {
        # which component does column k come from: its dpar/resp prefix
        pre <- ""
        for (cp in bk$components) {
          if (k > cp$offset && k <= cp$offset + cp$dim) {
            lp <- fit$frame$linpreds[[cp$lp_key]]
            parts <- c(if (lp$dpar != "mu") lp$dpar,
                       if (length(fit$spec$responses) > 1) lp$resp)
            pre <- paste(parts, collapse = "_")
            cf <- gsub("[()]", "", cp$cnms[k - cp$offset])
          }
        }
        if (!length(bk$components)) cf <- gsub("[()]", "", colnames(M)[k])
        for (lv in rownames(M)) {
          col <- paste0("r_", bk$group_name, if (nzchar(pre)) paste0("__", pre),
                        "[", lv, ",", cf, "]")
          if (!col %in% colnames(ds$draws)) { ex <- c(ex, col); next }
          checked <- checked + 1
          if (!identical(unname(ds$draws[i, col]), unname(M[lv, k]))) wrong <- wrong + 1
        }
      }
    }
  }
  cat("  r_ cells checked", checked, " wrong", wrong, " names not found", length(unique(ex)),
      if (length(ex)) paste("e.g.", ex[1]), "\n")
}
