# Paste the logs a section of dev/famuni-findings.md rests on into it,
# verbatim, so that no number there is typed. A placeholder line
# PASTE-<NAME>, or an earlier generated block of that name, is replaced
# by a fenced block of the named logs.
# Usage: Rscript dev/famuni-fill.R
setwd("C:/Users/adf44/source/r/frmtmb-wt-famuni")
o <- "dev/famuni-out/"
rd <- function(f, drop = "^g[+][+] ") {
  p <- if (file.exists(file.path(o, f))) file.path(o, f) else f
  x <- readLines(p, warn = FALSE)
  x <- x[!grepl(drop, x)]
  c(paste0("[", f, "]"), x)
}
blk <- function(nm, ...) {
  c(paste0("<!-- BEGIN generated: ", nm, " -->"), "```text", ..., "```",
    paste0("<!-- END generated: ", nm, " -->"))
}
blocks <- list(
  DENSITY = function() rd("density-ref.txt"),
  MOMENTS = function() c(rd("moments.txt"), "", rd("gevmean-mp.txt")),
  BRMSFIT = function() rd("brmsfit.txt"),
  RECOVERY = function() c(rd("recovery-frechet-200-final.txt"), "",
                          rd("recovery-gev-200-final.txt"), "",
                          rd("recovery-dw-200.txt"), "",
                          rd("recovery-ctl-gauss-200.txt")),
  SUITE = function() {
    s <- "dev/famuni-suite-full/SUMMARY.txt"
    if (file.exists(s)) rd(s) else "(suite not run yet)"
  },
  GATED = function() {
    s <- "dev/famuni-suite-gated/SUMMARY.txt"
    if (file.exists(s)) rd(s) else "(gated tier not run yet)"
  }
)
f <- "dev/famuni-findings.md"
md <- readLines(f)
for (nm in names(blocks)) {
  b0 <- which(md == paste0("<!-- BEGIN generated: ", nm, " -->"))
  b1 <- which(md == paste0("<!-- END generated: ", nm, " -->"))
  i <- which(md == paste0("PASTE-", nm))
  new <- blk(nm, blocks[[nm]]())
  if (length(b0) == 1L && length(b1) == 1L) {
    md <- c(md[seq_len(b0 - 1L)], new, md[-seq_len(b1)])
  } else if (length(i) == 1L) {
    md <- c(md[seq_len(i - 1L)], new, md[-seq_len(i)])
  }
}
writeLines(md, f)
cat("filled\n")
