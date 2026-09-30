# Fills the GEN-P1-* blocks of dev/postfit2-findings.md (section 7c,
# punch round 1) from the logs, so no count in it is typed by hand.
#   Rscript dev/postfit2-findings-gen3.R
wt <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2"
lg <- file.path(wt, "dev/postfit2-log")
rd <- function(f) readLines(file.path(lg, f), warn = FALSE)
res <- function(files) {
  unlist(lapply(files, function(f) {
    p <- file.path(lg, f)
    if (!file.exists(p)) return(paste("MISSING", f))
    r <- grep("^RESULT", readLines(p, warn = FALSE), value = TRUE)
    if (!length(r)) paste("NO RESULT", f) else paste0(r, "   [", f, "]")
  }))
}
fails <- function(f) {
  x <- readLines(file.path(lg, f), warn = FALSE)
  unique(sub("^---- [a-z_]+ in ", "  fails: ", grep("^---- ", x, value = TRUE)))
}
x <- readLines(file.path(wt, "dev/postfit2-findings.md"), warn = FALSE)
put <- function(x, tag, val) {
  i <- which(x == tag)
  if (length(i) != 1L) stop("marker ", tag, " found ", length(i), " times")
  c(x[seq_len(i - 1L)], val, x[-seq_len(i)])
}
x <- put(x, "GEN-P1-SIMBOOT", c(rd("p1-simboot-base.txt"),
                                rd("p1-simboot-lane.txt")))
arms <- function(file) {
  c(paste0("## ", file),
    unlist(lapply(c("r0", "base", "lane"), function(a) {
      f <- sprintf("p1-%s-%s", a, file)
      c(paste0(a, ": ", res(f)), if (a != "lane") fails(f))
    })))
}
x <- put(x, "GEN-P1-TESTS", c(
  arms("core-test-ce-levels.R.txt"),
  arms("core-test-simulate-newdata.R.txt"),
  arms("sample-test-postfit-draws.R.txt")))
x <- put(x, "GEN-P1-MUTANT", c(tail(rd("p1-mutant-none.txt"), 3),
                               tail(rd("p1-mutant-keydrop.txt"), 3)))
pick <- function(f, pat) grep(pat, rd(f), value = TRUE)
x <- put(x, "GEN-P1-REV", substr(c(
  pick("p1-rev-adv-fit-lane.txt",
       "^A1 list[(]g='99'[)] anl=FALSE|^A0|^A2|^B |^D |^E "),
  pick("p1-rev-nested-lane.txt", "^C-nested"),
  pick("p1-rev-guards-lane.txt", "^nl"),
  pick("p1-rev-adv-draws-lane.txt", "unseen|mixed|crossed|nested"),
  pick("p1-rev-brms-lane.txt", "^1")), 1, 240))
x <- put(x, "GEN-P1-BRMS", rd("p1-brms.txt"))
x <- put(x, "GEN-P1-SMOOTH", c(rd("p1-smooth-base.txt"),
                               rd("p1-smooth-lane.txt")))
chk <- vapply(c("frmtmb", "frmtmb.sample"), function(p) {
  f <- file.path(wt, "dev/postfit2-check", p, paste0(p, ".Rcheck"),
                 "00check.log")
  if (!file.exists(f)) return(paste(p, "NO CHECK LOG"))
  s <- grep("^Status:", readLines(f, warn = FALSE), value = TRUE)
  paste(p, if (length(s)) s else "NO STATUS LINE")
}, "")
gated_core <- list.files(lg, pattern = "^p1-gated-core-test-brms-suite")
x <- put(x, "GEN-P1-TIERS", c(
  "core and frmtmb.sample, whole ungated suites (dev/postfit2-suite.sh):",
  rd("p1-suite-summary.txt"),
  "", "gated tier (dev/postfit2-gated.sh):", rd("p1-gated-tier.txt"),
  "", "frmtmb.sample, whole suite gated (dev/postfit2-sample-suite.sh):",
  rd("p1-sample-suite.txt"),
  "", "core brms-suite tier, gated:", res(gated_core),
  "", "extension files calling frm_bootstrap(), simulate() or refit(),",
  "this lane's core (dev/postfit2-ext-suite.sh lane):",
  rd("p1-ext-suite-lane.txt"),
  "", "the same files on the base build (dev/postfit2-ext-suite.sh base):",
  rd("p1-ext-suite-base.txt"),
  "", "R CMD check --as-cran (dev/postfit2-check.ps1):", chk))
writeLines(x, file.path(wt, "dev/postfit2-findings.md"))
cat("filled\n")
