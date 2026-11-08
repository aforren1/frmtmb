# Review of lane superassign: compare the probe results of the base and
# lane builds (dev/review-superassign-probe.R), element by element, with
# doubles as 17-significant-digit strings, and print what each probe saw.
# Usage: Rscript dev/review-superassign-probe-cmp.R <a.rds> <b.rds>
a <- commandArgs(trailingOnly = TRUE)
x <- readRDS(a[1])
y <- readRDS(a[2])
cat("lib: none (compares", basename(a[1]), "with", basename(a[2]), ")\n")
nd <- 0L
for (p in union(names(x), names(y))) {
  same <- identical(x[[p]], y[[p]])
  if (!same) nd <- nd + 1L
  cat(sprintf("PROBE %-3s identical: %s\n", p, same))
  if (!same) {
    fx <- unlist(x[[p]])
    fy <- unlist(y[[p]])
    nm <- union(names(fx), names(fy))
    for (k in nm) {
      if (!identical(fx[[k]], fy[[k]])) {
        cat("   DIFF", k, ":", fx[[k]], "|", fy[[k]], "\n")
      }
    }
  }
}
cat("\n-- what the probes saw (first file)\n")
cat("p1 outA:", x$p1$outA, " count:", x$p1$A[[1]], " tally:", x$p1$tally,
    " best:", x$p1$A[[3]]$val, " B count:", x$p1$B[[1]],
    " A2 out:", x$p1$outA2, " B2 count:", x$p1$B2[[1]], "\n")
for (i in seq_along(x$p2)) {
  cat("p2 fit", i, ": nonfinite", x$p2[[i]]$nonfinite %||% NA, " obj",
      x$p2[[i]]$obj %||% x$p2[[i]]$error, " evals",
      paste(x$p2[[i]]$evals, collapse = "/"), "\n")
}
cat("p3 rep1 obj", x$p3[[1]]$obj, " rep2 obj", x$p3[[2]]$obj, "\n")
for (i in seq_along(x$p4)) {
  cat("p4 call", i, ": value", x$p4[[i]]$value, "; warnings",
      length(x$p4[[i]]$w), "; failures",
      if (is.null(x$p4[[i]]$failures)) "NULL" else
        paste(x$p4[[i]]$failures$groups, collapse = ","), "\n")
}
cat("p5:", paste(names(unlist(x$p5)), unlist(x$p5), sep = "=",
                  collapse = " "), "\n")
cat("p6 rows: b1", nrow(x$p6$b1), " b2", nrow(x$p6$b2), "\n")
cat("p7 calls:", x$p7$calls, "\n")
for (nm in names(x$p8)) {
  e <- x$p8[[nm]]
  if (is.character(e)) e <- list(error = e)
  cat("p8", nm, ":", if (!is.null(e$error)) paste("error:", e$error) else
    paste("entries", length(e$entries), "lower", length(e$lower),
          "upper", length(e$upper)), "\n")
}
for (nm in names(x$p9)) {
  cat("p9", nm, ": obj", x$p9[[nm]]$obj %||% x$p9[[nm]]$error,
      " nonfinite", x$p9[[nm]]$nonfinite %||% NA, "\n")
}
cat("RESULT probes differing:", nd, "\n")
