# Fills the GEN-P3-* blocks of dev/postfit2-findings.md (section 7e,
# punch round 3) from the logs.
#   Rscript dev/postfit2-findings-gen5.R
wt <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2"
lg <- file.path(wt, "dev/postfit2-log")
rd <- function(f) readLines(file.path(lg, f), warn = FALSE)
res <- function(f) {
  r <- grep("^RESULT", rd(f), value = TRUE)
  if (!length(r)) paste("NO RESULT", f) else paste0(r, "   [", f, "]")
}
fails <- function(f) {
  unique(sub("^---- [a-z_]+ (at line [0-9]+ )?in ", "  fails: ",
             grep("^---- ", rd(f), value = TRUE)))
}
x <- readLines(file.path(wt, "dev/postfit2-findings.md"), warn = FALSE)
put <- function(x, tag, val) {
  i <- which(x == tag)
  if (length(i) != 1L) stop("marker ", tag, " found ", length(i), " times")
  c(x[seq_len(i - 1L)], val, x[-seq_len(i)])
}
x <- put(x, "GEN-P3-BRMS", c("R = 30:", rd("p3-brms.txt"), "R = 120:",
                             rd("p3-brms-120.txt")))
fs <- c("p3-ce-levels-r2.txt", "p3-postfit-draws-r2.txt",
        "p3-lane-test-ce-levels.R.txt", "p3-lane-test-ce-options.R.txt",
        "p3-lane-test-postfit-draws.R.txt",
        "p3-mmsplit-test-ce-levels.R.txt",
        "p3-mmsplit-test-ce-options.R.txt",
        "p3-mmsplit-test-postfit-draws.R.txt")
x <- put(x, "GEN-P3-TESTS", unlist(lapply(fs, function(f) {
  c(res(f), fails(f))
})))
x <- put(x, "GEN-P3-CHECK", vapply(c("frmtmb", "frmtmb.sample"), function(p) {
  f <- file.path(wt, "dev/postfit2-check", p, paste0(p, ".Rcheck"),
                 "00check.log")
  s <- grep("^Status:", readLines(f, warn = FALSE), value = TRUE)
  paste(p, format(file.mtime(f), "%H:%M"),
        if (length(s)) s else "NO STATUS LINE")
}, ""))
writeLines(x, file.path(wt, "dev/postfit2-findings.md"))
cat("filled\n")
