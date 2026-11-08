# Each remeasured nl pair (dev/nlstart-equiv.sh) against the record
# (dev/equiv-data/results.tsv): the verdict columns and the numbers they
# rest on. Usage: Rscript dev/nlstart-equiv-cmp.R dev/nlstart-equiv-out
out <- commandArgs(TRUE)[1]
rd <- function(f) utils::read.delim(f, colClasses = "character",
                                    na.strings = "NA")
res <- rd("dev/equiv-data/results.tsv")
cols <- c("verdict", "level", "agree", "frm_ok", "brms_ok", "l1_ok",
          "l2_ok", "l4_ok")
nums <- c("l1_const", "l1_grad", "l2_metric", "l4_z")
fs <- list.files(out, "^p[0-9]+[.]tsv$", full.names = TRUE)
for (f in fs) {
  now <- rd(f)
  id <- now$id[1]
  rec <- res[res$id == id, ]
  same <- vapply(cols, function(k) identical(now[[k]][1], rec[[k]][1]), NA)
  cat(sprintf("%s %s x %s: verdict %s -> %s; %s\n", id, rec$feature_a,
              rec$feature_b, rec$verdict, now$verdict,
              if (all(same)) "verdict columns unchanged" else
                paste("CHANGED:", paste(cols[!same], collapse = ","))))
  for (k in nums) {
    a <- suppressWarnings(as.numeric(rec[[k]][1]))
    b <- suppressWarnings(as.numeric(now[[k]][1]))
    if (!identical(is.na(a), is.na(b)) || (!is.na(a) && a != b)) {
      cat(sprintf("    %s: %s -> %s\n", k, rec[[k]][1], now[[k]][1]))
    }
  }
  if (!is.na(now$frm_error[1]) && nzchar(now$frm_error[1])) {
    cat("    frm_error:", substr(now$frm_error[1], 1, 160), "\n")
  }
}
cat("pairs compared:", length(fs), "\n")
