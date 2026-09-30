# Fills the GEN-P2-* blocks of dev/postfit2-findings.md (section 7d,
# punch round 2) from the logs, so no count in it is typed by hand.
#   Rscript dev/postfit2-findings-gen4.R
wt <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2"
lg <- file.path(wt, "dev/postfit2-log")
rd <- function(f) readLines(file.path(lg, f), warn = FALSE)
res <- function(f) {
  p <- file.path(lg, f)
  if (!file.exists(p)) return(paste("MISSING", f))
  r <- grep("^RESULT", readLines(p, warn = FALSE), value = TRUE)
  if (!length(r)) paste("NO RESULT", f) else paste0(r, "   [", f, "]")
}
fails <- function(f) {
  x <- readLines(file.path(lg, f), warn = FALSE)
  unique(sub("^---- [a-z_]+ (at line [0-9]+ )?in ", "  fails: ",
             grep("^---- ", x, value = TRUE)))
}
x <- readLines(file.path(wt, "dev/postfit2-findings.md"), warn = FALSE)
put <- function(x, tag, val) {
  i <- which(x == tag)
  if (length(i) != 1L) stop("marker ", tag, " found ", length(i), " times")
  c(x[seq_len(i - 1L)], val, x[-seq_len(i)])
}
cut80 <- function(s) ifelse(nchar(s) > 200, paste0(substr(s, 1, 197), "..."),
                            s)
x <- put(x, "GEN-P2-CHECKS", cut80(c(rd("p2-checks-r1.txt"),
                                     rd("p2-checks-lane.txt"))))
x <- put(x, "GEN-P2-BRMS", cut80(rd("p2-brms.txt")))
arm2 <- function(pkg, file) {
  c(paste0("## ", pkg, " ", file),
    paste0("r1: ", res(sprintf("p2-%s-r1.txt", file))),
    fails(sprintf("p2-%s-r1.txt", file)),
    paste0("lane: ", res(sprintf("p2-%s-lane.txt", file))))
}
x <- put(x, "GEN-P2-TESTS", c(arm2("frmtmb", "ce-levels"),
                              arm2("frmtmb.sample", "postfit-draws")))
mfiles <- c(sprintf("p2-mut2-core-%s.txt",
                    c("byall", "mmskip", "mmgv", "rekey", "nosmooth")),
            sprintf("p2-mut2-sample-%s.txt", c("byall", "mmskip", "mmgv")),
            "p2-mut-p1_nolock.txt", "p2-mut-p1_noshare.txt",
            "p2-mut-sample-p1_noshare.txt")
x <- put(x, "GEN-P2-MUTANT", unlist(lapply(mfiles, function(f) {
  c(res(f), fails(f),
    sub("^ *[0-9]+ ", "  fails: ",
        grep("^[0-9]+ ", rd(f), value = TRUE)))
})))
tiers <- c("target", "core", "gated", "sample", "ext")
x <- put(x, "GEN-P2-TIERS", c(
  unlist(lapply(tiers, function(t) rd(sprintf("p2-tier-%s.txt", t)))),
  "",
  paste("the core tier's three failing files with frmtmb attached, as",
        "tests/testthat.R attaches it (dev/postfit2-rev-runtest.R):"),
  unlist(lapply(c("conditions", "data2", "id-kron"), function(f) {
    c(res(sprintf("p2-testenv-lane-%s.txt", f)),
      res(sprintf("p2-testenv-base-%s.txt", f)))
  })),
  "",
  "R CMD check --as-cran (dev/postfit2-check.ps1):",
  vapply(c("frmtmb", "frmtmb.sample"), function(p) {
    f <- file.path(wt, "dev/postfit2-check", p, paste0(p, ".Rcheck"),
                   "00check.log")
    if (!file.exists(f)) return(paste(p, "NO CHECK LOG"))
    s <- grep("^Status:", readLines(f, warn = FALSE), value = TRUE)
    paste(p, if (length(s)) s else "NO STATUS LINE")
  }, "")))
writeLines(x, file.path(wt, "dev/postfit2-findings.md"))
cat("filled\n")
