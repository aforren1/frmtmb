# Lane wcheck2: frm_compat()'s rows for asym_laplace and
# zero_inflated_asym_laplace, one line each.
#   Rscript dev/wcheck2-compat-ald.R <base|lane|lib>
lib_arg <- commandArgs(TRUE)[1]
source("dev/wcheck2-common.R")
for (f in c("asym_laplace", "zero_inflated_asym_laplace")) {
  cm <- tryCatch(frm_compat(f), error = function(e) conditionMessage(e))
  if (is.character(cm)) { cat(f, ": ", cm, "\n"); next }
  for (i in seq_len(nrow(cm))) {
    cat(sprintf("%s | %s (%s) | %s | %s\n", cm$feature_a[i], cm$feature_b[i],
                cm$kind_b[i], cm$status[i], substr(cm$note[i], 1, 70)))
  }
}
cat("DONE\n")
