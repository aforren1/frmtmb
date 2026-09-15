# Splice the generated blocks into dev/samplegen-findings.md. Each
# <<BLOCK:name>> marker is replaced by the file it names; a missing file
# leaves a visible PENDING line rather than a silent gap.
#   Rscript dev/samplegen-assemble.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-samplegen"
out <- file.path(root, "dev/samplegen-out")
fence <- function(x) c("```", x, "```")
rd <- function(f) {
  p <- file.path(out, f)
  if (file.exists(p)) readLines(p, warn = FALSE) else
    paste("PENDING:", f, "has not been produced")
}
strip_fence <- function(x) x[!x %in% c("```")]
blocks <- list(
  audit = fence(c("== ownership audit, dev/samplegen-audit.R ==",
                  grep("^<environment", rd("audit.txt"), value = TRUE,
                       invert = TRUE))),
  collision = rd("block-collision.txt"),
  brmsout = rd("block-brmsout.txt"),
  tests = rd("block-tests.txt"),
  loadcost = fence(c("== load cost, dev/samplegen-loadcost.R ==",
                     rd("loadcost.txt"))),
  shadowload = fence(c(
    "== partial owners, dev/samplegen-shadowload.R ==",
    "rellib-r3 is BASE, samplegen-lib is FIX", rd("shadowload.txt"))),
  run = c(
    "The test files this change can reach, one process each, two at a",
    "time, `NOT_CRAN` and `FRMTMB_BRMS_FIT_TESTS` true. That is every",
    "frmtmb.sample file, because the 28 names are used across them, and",
    "the core files that touch the sampling API or the bindings.",
    "",
    rd("block-suite-sample.txt"), "",
    rd("block-suite-core.txt"), "",
    "The same frmtmb.sample suite in ONE process, the way `R CMD check`",
    "runs it, which is the run that caught core's `hypothesis()` defect:",
    "",
    fence(c("== dev/samplegen-oneproc.R ==",
            grep("^(FILES|BLOCKS|ASSERT|FAIL|ERROR|SKIP|brms loaded|BAD)",
                 rd("oneproc-sample.txt"), value = TRUE))),
    "",
    "`R CMD check --as-cran`, once each, built with vignettes and without",
    "`--no-manual`, on a quiet machine (`dev/samplegen-cran.sh`):",
    "",
    fence(c("== dev/samplegen-cran.sh ==", unlist(lapply(
      c(frmtmb = "cran-core.txt", frmtmb.sample = "cran-sample.txt"),
      function(f) {
        x <- rd(f)
        i <- grep("^Status:", x)
        j <- grep("NOTE$|WARNING$|ERROR$", x)
        c(paste0("-- ", f), x[unique(sort(c(j, j + 1L, i)))])
      })))), "",
    "Before the check was allowed, the sections a NAMESPACE, signature or",
    "Rd change can break were run on the installed packages directly:",
    "",
    fence(rd("toolscheck.txt")), "",
    "Every changed Rd page, rendered with `Rd2txt` and grepped:",
    "",
    fence(rd("rdcheck.txt"))))
f <- file.path(root, "dev/samplegen-findings.md")
txt <- readLines(f, warn = FALSE)
for (nm in names(blocks)) {
  i <- which(txt == paste0("<<BLOCK:", nm, ">>"))
  if (length(i) != 1L) stop("marker for ", nm, " found ", length(i), " times")
  txt <- c(txt[seq_len(i - 1L)], blocks[[nm]], txt[-seq_len(i)])
}
con <- file(f, "wb")
writeBin(charToRaw(paste0(paste(txt, collapse = "\n"), "\n")), con)
close(con)
cat("assembled", length(blocks), "blocks\n")
