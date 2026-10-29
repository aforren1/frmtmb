# Release fix R4: each remeasured pair (dev/rel078-equiv-rerun.sh) against
# the equiv lane's own measurement on rellib-r15
# (../frmtmb-wt-equiv/dev/equiv-out/<id>.tsv, read only) and against the
# adjudicated record (dev/equiv-data/results.tsv): the verdict columns,
# and every number the verdict rests on (l1_const, l1_grad, l2_metric,
# l4_z), as the review's cmp-rerun compared them.
#   Rscript dev/rel078-equiv-cmp.R <outdir>
out <- commandArgs(TRUE)[1]
lane <- "C:/Users/adf44/source/r/frmtmb-wt-equiv/dev/equiv-out"
rd <- function(f) utils::read.delim(f, colClasses = "character",
                                    na.strings = "NA")
res <- rd("dev/equiv-data/results.tsv")
cols <- c("verdict", "level", "agree", "frm_ok", "brms_ok", "l1_ok",
          "l2_ok", "l4_ok", "method_reached")
nums <- c("l1_const", "l1_grad", "l2_metric", "l4_z")
fs <- sort(list.files(out, "^p[0-9]+[.]tsv$"))
moved <- 0L
for (f in fs) {
  id <- sub("[.]tsv$", "", f)
  a <- rd(file.path(out, f))
  b <- rd(file.path(lane, f))
  r <- res[res$id == id, ]
  cat("==", id, a$feature_a, "x", a$feature_b, "| frmtmb", a$frmtmb, "\n")
  cat("   rerun:", a$verdict, a$level, "| lane:", b$verdict, b$level,
      "| record:", r$verdict, r$level, "\n")
  dv <- cols[vapply(cols, function(k) !identical(a[[k]], b[[k]]), NA)]
  num <- vapply(nums, function(k) {
    x <- suppressWarnings(as.numeric(a[[k]]))
    y <- suppressWarnings(as.numeric(b[[k]]))
    if (is.na(x) && is.na(y)) return(0)
    if (is.na(x) || is.na(y)) return(Inf)
    abs(x - y) / max(abs(y), .Machine$double.eps)
  }, 0)
  if (length(dv)) {
    moved <- moved + 1L
    for (k in dv) cat("   MOVED", k, ":", b[[k]], "->", a[[k]], "\n")
  }
  cat("   numbers, relative difference to the lane's:",
      paste(sprintf("%s %.3g", nums, num), collapse = ", "), "\n")
  if (!identical(a$verdict, r$verdict)) {
    cat("   raw verdict differs from the record's adjudicated one:",
        r$verdict, "(adjudication.tsv)\n")
  }
}
cat("pairs:", length(fs), "| verdict columns moved against the lane:",
    moved, "\n")
