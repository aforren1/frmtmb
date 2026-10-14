# Paste the fix-round-1 logs into dev/famuni-findings.md section 12.
setwd("C:/Users/adf44/source/r/frmtmb-wt-famuni")
o <- "dev/famuni-out/"
rd <- function(f, keep = NULL, drop = "^g[+][+] |NaNs produced") {
  x <- readLines(file.path(o, f), warn = FALSE)
  x <- x[!grepl(drop, x)]
  if (!is.null(keep)) x <- x[grepl(keep, x)]
  c(paste0("[", f, "]"), x)
}
blk <- function(nm, ...) c(paste0("<!-- BEGIN generated: ", nm, " -->"),
                           "```text", ..., "```",
                           paste0("<!-- END generated: ", nm, " -->"))
final <- function() {
  fs <- c("final-installs.txt", "final-install-sample2.txt",
          "final-suite.txt", "final-gated.txt", "final-drmtmb.txt",
          "final-reach.txt", "final-openblas.txt")
  out <- character(0)
  for (f in fs) {
    if (!file.exists(file.path(o, f))) {
      out <- c(out, paste0("[", f, "] (not run yet)"))
      next
    }
    keep <- if (f == "final-drmtmb.txt") "^RESULT|^lib:" else NULL
    out <- c(out, rd(f, keep), "")
  }
  ck <- function(p) {
    fl <- file.path("dev/famuni-check", p, "check.log")
    if (!file.exists(fl)) return(paste(p, ": not run"))
    con <- file(fl, encoding = "UTF-16LE")
    x <- gsub("[^ -~]", "", readLines(con, warn = FALSE))
    close(con)
    paste0(p, " R CMD check --as-cran: ",
           paste(grep("^Status|NOTE$|WARNING$|ERROR$", x, value = TRUE),
                 collapse = " | "))
  }
  c(out, ck("frmtmb"), ck("frmtmb.eam"), ck("frmtmb.sample"))
}
blocks <- list(
  FOLD = function() rd("fold-probe.txt"),
  WGUARD = function() rd("fix1-wguard2.txt", keep = "^seed"),
  TMPXI = function() rd("tmpxi.txt"),
  TMPXISAMPLE = function() rd("tmpxi-sample.txt",
    keep = "===|divergences|outside the support|xi median|^xi "),
  FINAL = final)
f <- "dev/famuni-findings.md"
md <- readLines(f)
for (nm in names(blocks)) {
  b0 <- which(md == paste0("<!-- BEGIN generated: ", nm, " -->"))
  b1 <- which(md == paste0("<!-- END generated: ", nm, " -->"))
  i <- which(md == paste0("PASTE-", nm))
  new <- blk(nm, blocks[[nm]]())
  if (length(b0) == 1L) md <- c(md[seq_len(b0 - 1L)], new, md[-seq_len(b1)])
  else if (length(i) == 1L) md <- c(md[seq_len(i - 1L)], new, md[-seq_len(i)])
}
writeLines(md, f)
cat("filled\n")
