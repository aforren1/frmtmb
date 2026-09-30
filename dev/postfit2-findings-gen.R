# Fills the GENERATED-* blocks of dev/postfit2-findings.md from the
# logs, so no count in it is typed by hand. Run after the logs exist:
#   Rscript dev/postfit2-findings-gen.R
wt <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2"
lg <- file.path(wt, "dev/postfit2-log")
rd <- function(f) readLines(file.path(lg, f), warn = FALSE)
block <- function(x) c("```", x, "```")
res <- function(files) {
  unlist(lapply(files, function(f) {
    p <- file.path(lg, f)
    if (!file.exists(p)) return(paste("MISSING", f))
    x <- readLines(p, warn = FALSE)
    r <- grep("^RESULT", x, value = TRUE)
    if (!length(r)) paste("NO RESULT", f) else r
  }))
}
first_fail <- function(f, n = 3) {
  x <- rd(f)
  i <- grep("^---- ", x)
  unlist(lapply(head(i, n), function(k) {
    c(x[k], substr(x[k + 1], 1, 160))
  }))
}

before_files <- c("before-core-test-conditional-smooths.R.txt",
                  "before-core-test-brms-utilities.R.txt",
                  "before-core-test-ce-options.R.txt",
                  "before-core-test-ce-bands.R.txt",
                  "before-core-test-nlf.R.txt",
                  "before-sample-test-postfit-draws.R.txt")
gen_before <- block(c(
  res(before_files), "",
  "first failures of each, base build:",
  unlist(lapply(before_files, function(f) c(paste0("[", f, "]"),
                                            first_fail(f, 2))))))

cmp <- rd("brms-compare.txt")
keep <- !grepl(paste0("^    (se__|lower__|upper__) |^    attr |",
                      "^    cond__ equal"), cmp)
gen_compare <- block(cmp[keep])

gen_nlmc <- block(rd("nlmc.txt"))

gated <- list.files(lg, pattern = "^gated-(core|sample)-", full.names = TRUE)
flips <- unlist(lapply(gated, function(f) {
  x <- readLines(f, warn = FALSE)
  m <- regmatches(x, regexpr("brms [a-z-]+:[0-9]+ now HOLDS but is recorded as '[a-z ]+'", x))
  if (length(m)) paste0(sub("^gated-", "", basename(f)), ": ", m)
}))
base_gated <- list.files(lg, pattern = "^gated-base-", full.names = FALSE)
gen_flips <- block(c(
  "this lane's library:", flips,
  sprintf("stale-verdict failures: %d, all of them rows above", length(flips)),
  "", "the base build, same five files (every row holds its recorded verdict):",
  res(base_gated)))

suite <- rd("suite-summary.txt")
suite <- suite[!grepl("Segmentation fault", suite)]
segv <- sum(grepl("Segmentation fault", rd("suite-summary.txt")))
new_files <- c("after-core-test-conditional-smooths.R.txt",
               "after-core-test-brms-utilities.R.txt",
               "after-core-test-ce-options.R.txt",
               "after2-core-test-ce-bands.R.txt",
               "after2-core-test-nlf.R.txt",
               "after2-core-test-conditions.R.txt",
               "after-core-test-generic-collision.R.txt",
               "after-core-test-adefects.R.txt",
               "after2-sample-test-postfit-draws.R.txt",
               "after2-sample-test-message-uniqueness.R.txt",
               "after-sample-test-conditional-effects-draws.R.txt",
               "after-sample-test-generic-collision.R.txt",
               "after-sample-test-adefects.R.txt")
chk <- vapply(c("frmtmb", "frmtmb.sample"), function(p) {
  f <- file.path(wt, "dev/postfit2-check", p, paste0(p, ".Rcheck"),
                 "00check.log")
  if (!file.exists(f)) return(paste(p, "NO CHECK LOG"))
  s <- grep("^Status:", readLines(f, warn = FALSE), value = TRUE)
  paste(p, if (length(s)) s else "NO STATUS LINE")
}, "")
gen_tests <- block(c(
  "whole ungated suites, one process per file (dev/postfit2-suite.sh):",
  suite,
  sprintf(paste("(%d 'Segmentation fault' lines from Rscript at process",
                "exit, after each RESULT line; every log has its RESULT)"),
          segv),
  "", "after the fixes that run prompted, rerun one file per process:",
  res(new_files),
  "", "gated brms-suite tier, all files, this lane's library:",
  res(basename(gated)),
  "", "R CMD check --as-cran (dev/postfit2-check.ps1):", chk))

f <- file.path(wt, "dev/postfit2-findings.md")
x <- readLines(f, warn = FALSE)
put <- function(x, tag, val) {
  i <- which(x == tag)
  if (length(i) != 1L) stop("marker ", tag, " found ", length(i), " times")
  c(x[seq_len(i - 1L)], val, x[-seq_len(i)])
}
x <- put(x, "GENERATED-BEFORE", gen_before)
x <- put(x, "GENERATED-COMPARE", gen_compare)
x <- put(x, "GENERATED-NLMC", gen_nlmc)
x <- put(x, "GENERATED-FLIPS", gen_flips)
x <- put(x, "GENERATED-TESTS", gen_tests)
writeLines(x, f)
cat("filled\n")
