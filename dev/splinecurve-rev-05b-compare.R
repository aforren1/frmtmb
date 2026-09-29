# Reviewer, lane splinecurve: compare the two arms of
# splinecurve-rev-05-bitwise.R.
#   Rscript splinecurve-rev-05b-compare.R <before.rds> <after.rds>
a <- commandArgs(trailingOnly = TRUE)
b <- readRDS(a[1])
f <- readRDS(a[2])
stopifnot(identical(names(b), names(f)))
drop_flag <- function(v) {
  s <- attr(v, "spec")
  if (!is.null(s)) {
    s$allow_new_levels <- NULL
    attr(v, "spec") <- s
  }
  v
}
n_raw <- 0L
n_ok <- 0L
for (nm in names(b)) {
  raw <- identical(b[[nm]], f[[nm]])
  vb <- b[[nm]]$value
  vf <- f[[nm]]$value
  ok <- identical(drop_flag(vb), drop_flag(vf)) &&
    identical(b[[nm]]$warnings, f[[nm]]$warnings)
  n_raw <- n_raw + raw
  n_ok <- n_ok + ok
  cat(sprintf(paste0("%-24s identical raw: %-5s  without ",
                     "spec$allow_new_levels: %s%s\n"),
              nm, raw, ok,
              if (inherits(vb, "err")) paste0("  (both ERROR: ",
                                             substr(unclass(vb), 1, 60), ")")
              else ""))
  if (!ok) {
    print(all.equal(drop_flag(vb), drop_flag(vf)))
    if (inherits(vf, "err")) cat("  after arm ERROR:", unclass(vf), "\n")
  }
}
cat("\n", length(b), "calls;", n_raw, "identical raw;", n_ok,
    "identical once spec$allow_new_levels is removed\n")
flags <- vapply(f, function(x) {
  s <- attr(x$value, "spec")
  if (is.null(s)) NA else isTRUE(s$allow_new_levels)
}, logical(1))
cat("after arm spec$allow_new_levels values:",
    paste(names(table(flags, useNA = "ifany")), table(flags, useNA = "ifany"),
          sep = "=", collapse = " "), "\n")
