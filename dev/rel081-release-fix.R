# The brief's merge steps for frmtmb 0.81.0 that edit text, applied once
# to the merged release tree (after dev/rel081-merge.sh):
#
#   Rscript dev/rel081-release-fix.R
#
# Each edit must match exactly once, or the script stops before writing
# anything more. A file keeps its own line endings.
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel081"
R <- "C:/Users/adf44/source/r"
rd <- function(p) {
  raw <- readBin(p, "raw", file.info(p)$size)
  s <- rawToChar(raw)
  list(text = gsub("\r\n", "\n", s, fixed = TRUE),
       crlf = grepl("\r\n", s, fixed = TRUE))
}
wr <- function(p, f) {
  s <- if (f$crlf) gsub("\n", "\r\n", f$text, fixed = TRUE) else f$text
  writeBin(charToRaw(s), p)
}
fix <- function(rel, old, new) {
  p <- file.path(REL, rel)
  f <- rd(p)
  n <- lengths(regmatches(f$text, gregexpr(old, f$text, fixed = TRUE)))
  if (n != 1L) stop(rel, ": the edit matched ", n, " times: ", substr(old, 1, 60))
  f$text <- sub(old, new, f$text, fixed = TRUE)
  wr(p, f)
  cat("edited", rel, "\n")
}

## step 4: the new importance and nested-quadrature test takes about 20 s
fix("tests/testthat/test-families-extreme.R",
'test_that("importance and nested quadrature fit a GEV off its support", {
',
'test_that("importance and nested quadrature fit a GEV off its support", {
  # an importance = 500 fit and a nested quadrature fit, about 20 s, so
  # off CRAN as test-importance.R keeps its long fits
  skip_on_cran()
')

## step 4: the support record names a von_mises() bound past pi too
## (the short check's note 1, dev/review-support-fix2.R)
fix("dev/support-findings.md",
"    von_mises trunc(lb = -4)                 dlogLik 5.68e-14  max|dbeta| 4.56e-12
    von_mises trunc(lb = -pi)                dlogLik 5.68e-14  max|dbeta| 4.56e-12
",
"    von_mises trunc(lb = -4)                 dlogLik 5.68e-14  max|dbeta| 4.56e-12
    von_mises trunc(lb = -pi)                dlogLik 5.68e-14  max|dbeta| 4.56e-12

A bound past pi was as wrong on the base and is fixed by the same clamp
(added at consolidation from the review's short check,
`dev/review-support-fix2.R`, seed 406, n 300; dlogLik and max|dpar|):

    von_mises trunc(ub = 4)                  base -38.6 (0.262)   fixed 5.7e-14 (4.6e-12)
    von_mises trunc(lb = -4, ub = 4)         base -85.2 (0.738)   fixed 5.7e-14 (4.6e-12)
    von_mises cens + trunc(lb = -4, ub = 4)  base -188 (69.1)     fixed 0 (1.6e-15)
")
fix("dev/support-findings.md",
"| `R/families-cdf.R`, `von_mises_cdf_log()` at q outside [-pi, pi] | S | von_mises | an arc of negative length; logLik -6.42 (measured) | fit | FIXED, tape and doubles |",
"| `R/families-cdf.R`, `von_mises_cdf_log()` at q outside [-pi, pi] | S | von_mises | an arc of negative length; logLik -6.42 at `lb = -4`, -38.6 at `ub = 4`, -188 with `cens()` and both (measured) | fit | FIXED, tape and doubles |")

## step 3 and survey step 4: versions and floors
fix("DESCRIPTION", "Version: 0.80.0", "Version: 0.81.0")
fix("extensions/frmtmb.sample/DESCRIPTION", "Version: 0.27.0", "Version: 0.28.0")
fix("extensions/frmtmb.sample/DESCRIPTION", "frmtmb (>= 0.80.0)", "frmtmb (>= 0.81.0)")
fix("extensions/frmtmb.survival/DESCRIPTION", "Version: 0.2.2", "Version: 0.2.3")

## step 8: the release harness reads the round's library
fix("dev/release/run-tests.R", 'LIB <- "C:/Users/adf44/source/r/rellib-r19"',
    'LIB <- "C:/Users/adf44/source/r/rellib-r20"')
fix("dev/release/run-session.R", 'LIB <- "C:/Users/adf44/source/r/rellib-r19"',
    'LIB <- "C:/Users/adf44/source/r/rellib-r20"')

## step 1: the translator's adjudication edits (a DROP row removes the
## id's row; any other row replaces it)
adj_p <- file.path(REL, "dev/equiv-data/adjudication.tsv")
ed_p <- file.path(R, "frmtmb-wt-translator/dev/translator-adjudication-edits.tsv")
adj <- rd(adj_p)
lines <- strsplit(adj$text, "\n", fixed = TRUE)[[1L]]
ed <- strsplit(gsub("\r", "", rawToChar(readBin(ed_p, "raw", file.info(ed_p)$size))),
               "\n", fixed = TRUE)[[1L]]
ed <- ed[nzchar(ed)][-1L]
ids <- sub("\t.*", "", lines)
changed <- 0L
for (e in ed) {
  id <- sub("\t.*", "", e)
  k <- which(ids == id)
  if (length(k) != 1L) stop("adjudication: ", id, " matched ", length(k), " rows")
  if (grepl("^[^\t]+\tDROP", e)) {
    lines <- lines[-k]; ids <- ids[-k]
  } else {
    lines[k] <- e
  }
  changed <- changed + 1L
}
adj$text <- paste0(paste(lines, collapse = "\n"), "\n")
wr(adj_p, adj)
cat("adjudication.tsv: ", changed, " rows changed; ids now: ",
    paste(ids[-1L], collapse = " "), "\n", sep = "")
cat("RELEASE FIX DONE\n")
