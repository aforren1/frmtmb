# Dump extracted + transformed code for eyeballing. Rscript dump.R <vignette>
args <- commandArgs(trailingOnly = TRUE)
root <- "C:/Users/adf44/source/r/frmtmb-wt-audit"
source(file.path(root, "dev/brms-port/port-lib.R"))
for (v in args) {
  ch <- extract_vignette(v)
  cat("########## ", v, " : ", length(ch), " chunks\n", sep = "")
  for (k in ch) {
    tr <- transform_code(k$code)
    for (j in seq_along(tr)) {
      t <- tr[[j]]
      cat(sprintf("--- [%s.%d.%d] %s%s\n", v, k$idx, j, t$status,
                  if (length(t$dropped)) paste0(" drop:", paste(t$dropped, collapse = ",")) else ""))
      cat(t$src, "\n")
    }
  }
}
