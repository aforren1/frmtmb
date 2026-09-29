# Lane wt-arcovsample, punch round 1: the mechanical nits, in one pass.
# Written as a FILE because `Rscript -e` with backticks segfaulted twice
# on this box in round 0.
#
#   Rscript dev/arcovsample-punch1.R

wr <- function(p, t) {
  con <- file(p, open = "wb")
  writeBin(charToRaw(paste0(paste(t, collapse = "\r\n"), "\r\n")), con)
  close(con)
}

# nit 1 and 2: an 86-column test_that() title, and no blank line before
# the test_that() that follows the new block
p <- "tests/testthat/test-mv-gaps.R"
t <- readLines(p)
i <- grep('^test_that\\("the pointwise rescor density holds at a large nu, as the objective does", \\{$',
          t)
stopifnot(length(i) == 1L)
t[i] <- 'test_that("the pointwise rescor density holds at a large nu", {'
j <- grep('^test_that\\("student\\(\\) rescor refuses what brms refuses", \\{$',
          t)
stopifnot(length(j) == 1L)
if (nzchar(t[j - 1L])) t <- append(t, "", after = j - 1L)
wr(p, t)
cat(p, ": title now", nchar(t[i]), "columns; blank line before the next",
    "test_that()\n")
stopifnot(sum(nchar(readLines(p)) > 80) == 0L)

# nit 3: skip_sampler() in the new brms block, which builds its draws
# with stanfit = NULL and needs no sampler
p <- "extensions/frmtmb.sample/tests/testthat/test-loo.R"
t <- readLines(p)
i <- grep("^test_that\\(\"log_lik\\(\\) matches brms row by row on a cov = FALSE ARMA\", \\{$",
          t)
stopifnot(length(i) == 1L,
          identical(t[i + 1L], "  skip_unless_brms_fit()"),
          identical(t[i + 2L], "  skip_sampler()"))
t[i + 2L] <- paste0("  # no skip_sampler(): the draws are a fixed ",
                    "matrix with stanfit = NULL,")
t <- append(t, "  # so tmbstan is not on this block's path at all",
            after = i + 2L)
wr(p, t)
cat(p, ": skip_sampler() removed from the brms block\n")
# NOT a whole-file width check here: lines 86 and 264 of this file are
# 85 columns and are PRE-EXISTING (`git show HEAD:` has them at the same
# line numbers). Only the lines this script writes are checked.
stopifnot(nchar(t[i + 2L]) <= 80L, nchar(t[i + 3L]) <= 80L)
cat("DONE\n")
