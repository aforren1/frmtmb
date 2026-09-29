# REVIEW script 04b: compare the two builds' saved results bitwise.
#
#   Rscript dev/arcovsample-rev-04b-cmp.R

a <- readRDS("dev/arcovsample-rev-log/pp-lane.rds")
b <- readRDS("dev/arcovsample-rev-log/pp-ref.rds")
cat("lane newexports ", a$newexports, "  ref newexports ", b$newexports,
    "\n\n", sep = "")
stopifnot(isTRUE(a$newexports), isFALSE(b$newexports))

for (nm in grep("^pp_", names(a), value = TRUE)) {
  for (k in c("predict", "epred")) {
    x <- a[[nm]][[k]]; y <- b[[nm]][[k]]
    if (is.character(x) || is.character(y)) {
      cat(nm, " ", k, ": lane=", if (is.character(x)) "ERR" else "ok",
          " ref=", if (is.character(y)) "ERR" else "ok",
          " same message: ", identical(x, y), "\n", sep = "")
      if (is.character(x)) cat("    lane: ", x, "\n", sep = "")
      if (is.character(y)) cat("    ref : ", y, "\n", sep = "")
      next
    }
    cat(nm, " ", k, ": identical=", identical(x, y),
        "  max|diff|=", format(max(abs(x - y)), digits = 6),
        "  dim ", paste(dim(x), collapse = "x"), "\n", sep = "")
  }
}

cat("\n---- refusal messages, lane against reference ----\n")
for (nm in grep("^msg_", names(a), value = TRUE)) {
  same <- identical(a[[nm]], b[[nm]])
  cat(nm, ": identical=", same, "\n", sep = "")
  if (!same) {
    cat("  LANE: ", a[[nm]], "\n", sep = "")
    cat("  REF : ", b[[nm]], "\n", sep = "")
  }
}
cat("\nDONE\n")
