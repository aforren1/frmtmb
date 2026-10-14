# Fix round 1: can the ported bin-2P rows see an error in FRMTMB's
# density, as opposed to one in the adapter? (The weights construction
# of dev/remeasure-bin2-absent.R changes the adapter's own arithmetic,
# so it shows only that log_lik:32 tests the adapter.)
#
#   Rscript dev/remeasure-bin2-sayno.R \
#     > dev/remeasure-out/probe/bin2-sayno.txt
#
# A scratch copy of dev/brmsport-bin2-adapters.R stands in for a frmtmb
# whose row density is wrong by a relative 1e-6: frmtmb::row_lpdf()'s
# value, and each ordinal family's ord_lpdf_make() row, is multiplied by
# 1 + 1e-6 before the adapter returns it. Nothing else changes. Every
# log_lik row that holds on the true build is recorded again. A row
# that compares a value against brms's or base R's density must stop
# holding; a row that checks only a length or a dim cannot see it, and
# the count says how many of each.
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = "frmtmb", FRMTMB_PORT_ROOT = getwd())
stage <- "dev/remeasure-out/sayno-stage"
unlink(stage, recursive = TRUE)
dir.create(stage, recursive = TRUE)
file.copy(list.files("tests/testthat", "^helper.*[.]R$", full.names = TRUE),
          stage)
ad <- readLines("dev/brmsport-bin2-adapters.R")
h1 <- grep("ll <- frmtmb::row_lpdf(fam, y, y, dp$dpars, av, NULL)", ad,
           fixed = TRUE)
h2 <- grep("lp(y[s], dps, list(), list(tau_raw = th))", ad, fixed = TRUE)
stopifnot(length(h1) == 1L, length(h2) == 1L)
ad[h1] <- sub("NULL)", "NULL) * (1 + 1e-6)", ad[h1], fixed = TRUE)
ad[h2] <- sub("list(tau_raw = th))", "list(tau_raw = th)) * (1 + 1e-6)",
              ad[h2], fixed = TRUE)
cat("perturbed adapter lines:\n", ad[h1], "\n", ad[h2], "\n")
writeLines(ad, file.path(stage, "helper-brms-bin2-adapters.R"))
f <- "dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R"
file.copy(f, stage)
rec <- normalizePath(file.path(stage, "rec.tsv"), winslash = "/",
                     mustWork = FALSE)
Sys.setenv(FRMTMB_BRMSPORT_RECORD = rec)
suppressMessages({library(testthat); library(frmtmb)})
invisible(test_file(file.path(stage, basename(f)), package = "frmtmb",
                    env = testthat::test_env("frmtmb"), reporter = "silent",
                    stop_on_failure = FALSE))
rd <- function(p) {
  x <- utils::read.delim(p, header = FALSE, quote = "",
                         colClasses = "character")
  x[x[[1]] == "assert", c(3, 5)]
}
true <- rd("dev/brmsport-bin2-log/rec-frmtmb-log_lik.tsv")
mut <- rd(rec)
led <- utils::read.delim("dev/brmsport-bin2-ledger.tsv", quote = "",
                         stringsAsFactors = FALSE)
m <- merge(true, mut, by = "V3", suffixes = c("_true", "_mut"))
m$kind <- led$kind[match(m$V3, led$id)]
m <- m[m$V5_true == "TRUE" & led$outcome[match(m$V3, led$id)] == "pass", ]
cat("\nlog_lik passes on the true build:", nrow(m), "\n")
print(table(kind = m$kind, holds_with_wrong_density = m$V5_mut))
cat("\nvalue rows still holding with a density off by 1e-6:",
    paste(m$V3[m$kind == "value" & m$V5_mut == "TRUE"], collapse = " "),
    "\n")
