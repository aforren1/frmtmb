# The generated blocks of dev/latres-findings.md, from the lane's logs: L1 over
# term x family (dev/latres-l1.R), L1 beside other features (dev/latres-l1x.R),
# and the grid's counts before and after (dev/latres-grid-summary.R).
rd <- function(pat) {
  fs <- list.files("dev/latres-out", pat, full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    x <- readLines(f, warn = FALSE)
    h <- strsplit(x[1L], "\t", fixed = TRUE)[[1L]]
    m <- lapply(strsplit(x[-1L], "\t", fixed = TRUE), function(v) {
      c(v, rep("", length(h) - length(v)))[seq_along(h)]
    })
    d <- as.data.frame(do.call(rbind, m), stringsAsFactors = FALSE)
    names(d) <- h
    d
  }))
}
num <- function(x) suppressWarnings(as.numeric(x))
l1 <- rd("^l1-w[0-9]+[.]tsv$")
cat("### L1 over term x family\n\n")
cat(nrow(l1), "cells. Each cell: frm() fit, then brms's flat-prior program.\n\n")
worst <- function(d) {
  cc <- abs(c(num(d$const), num(d$const_up), num(d$const_dn)))
  cc[is.finite(cc)]
}
cat(sprintf("Largest |constant| over all cells and points: %.3g.\n",
            max(worst(l1))))
cat(sprintf("Largest inner gradient at the estimate: %.3g.\n",
            max(num(l1$grad), na.rm = TRUE)))
cat(sprintf(paste0("Largest |logLik(fit) - Laplace value on brms's program ",
                   "(by differences of Stan's gradient in zerr)|: %.3g.\n"),
            max(abs(num(l1$laplace_gap)), na.rm = TRUE)))
na_pts <- sum(is.na(num(l1$const_up))) + sum(is.na(num(l1$const_dn)))
cat("Stepped points brms's Stan program rejects (no constant):", na_pts, "\n")
cat("Cells with an error instead of a number:",
    sum(!nzchar(l1$const) | is.na(num(l1$const))), "\n")
cat("frm() convergence codes, cells by code:\n")
print(table(l1$conv))
cat("\n| term | ", paste(unique(l1$family), collapse = " | "), " |\n", sep = "")
cat("|---|", paste(rep("---", length(unique(l1$family))), collapse = "|"),
    "|\n", sep = "")
for (tm in unique(l1$term)) {
  v <- vapply(unique(l1$family), function(f) {
    r <- l1[l1$term == tm & l1$family == f, ]
    if (!nrow(r)) return("")
    sprintf("%.1e", max(worst(r)))
  }, "")
  cat("| ", tm, " | ", paste(v, collapse = " | "), " |\n", sep = "")
}
lx <- rd("^l1x-w[0-9]+[.]tsv$")
cat("\n### L1 beside other features\n\n")
cat("| spec | const | up | down | grad | logLik | conv | note |\n")
cat("|---|---|---|---|---|---|---|---|\n")
lx <- lx[order(lx$spec), ]
for (i in seq_len(nrow(lx))) {
  r <- lx[i, ]
  cat("| ", paste(c(r$spec, r$const, r$const_up, r$const_dn, r$grad, r$loglik,
                    r$conv, substr(gsub("[;|]", " ", r$note), 1L, 90L)),
                  collapse = " | "), " |\n", sep = "")
}
