# Markdown table of dev/drmtmb-log/agree.rds for dev/drmtmb-findings.md.
res <- readRDS(
  "C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-log/agree.rds")
res2 <- readRDS(
  "C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-log/agree2.rds")
s <- do.call(rbind, lapply(c(res, res2), `[[`, "summary"))
f <- function(x) formatC(x, format = "e", digits = 1)
cat("| Model | logLik frmtmb | logLik drmTMB | diff | gap at drm opt |",
    "gap at frm opt | gap off opt | max diff/SE | max SE ratio - 1 |\n")
cat("|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(s))) {
  r <- s[i, ]
  cat(sprintf("| %s | %.6f | %.6f | %s | %s | %s | %s | %s | %s |\n",
              gsub("|p|", "labelled", r$model, fixed = TRUE),
              r$ll_frm, r$ll_drm, f(r$ll_diff), f(r$obj_gap_at_drm_opt),
              f(r$obj_gap_at_frm_opt), f(r$obj_gap_off_opt),
              f(r$max_abs_diff_over_se),
              f(expm1(r$max_abs_log_se_ratio))))
}
