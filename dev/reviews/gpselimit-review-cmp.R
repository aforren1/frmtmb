# Review of lane gpselimit, check 3: base against lane, field by field,
# for the 30-fit sweep and the review's cases (separate processes).
# Usage: Rscript gpselimit-review-cmp.R
out <- "C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/reviews/gpselimit-review-out"
cmp <- function(tag, b, l) {
  if (is.null(b) || is.null(l)) { cat("CMP", tag, "missing\n"); return() }
  se_same <- identical(unname(b$se), unname(l$se))
  rel <- suppressWarnings(max(abs(l$se - b$se) / abs(b$se), na.rm = TRUE))
  nanchg <- which(is.nan(b$se) != is.nan(l$se))
  fin <- is.finite(b$se) & is.finite(l$se)
  cat(sprintf(paste0("CMP %-26s est %d par %d ll %d se %d lost %d | ",
                     "se max rel diff (both finite) %.3g | NaN changed at ",
                     "[%s] base %s lane %s\n"),
              tag, identical(b$est, l$est), identical(b$par, l$par),
              identical(b$ll, l$ll), se_same,
              identical(unname(b$lost), unname(l$lost)),
              if (any(fin)) max(abs(l$se[fin] - b$se[fin]) / b$se[fin]) else NA,
              paste(nanchg, collapse = ","),
              paste(format(b$se[nanchg], digits = 6), collapse = ","),
              paste(format(l$se[nanchg], digits = 6), collapse = ",")))
  if (any(fin) && !se_same) {
    cat("   base", paste(names(b$se), format(b$se, digits = 17), collapse = " "), "\n")
    cat("   lane", paste(names(l$se), format(l$se, digits = 17), collapse = " "), "\n")
  }
}
sb <- readRDS(file.path(out, "sweep-base-j0.rds"))
sl <- readRDS(file.path(out, "sweep-lane-j0.rds"))
for (i in 1:30) cmp(sprintf("sweep%02d", i), sb[[i]], sl[[i]])
for (f in c("ABCDEFGHI")) {
  cb <- readRDS(file.path(out, sprintf("cases-base-%s.rds", f)))
  cl <- readRDS(file.path(out, sprintf("cases-lane-%s.rds", f)))
  for (n in union(names(cb), names(cl))) cmp(n, cb[[n]], cl[[n]])
}
