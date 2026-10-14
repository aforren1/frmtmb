# Reviewer round 3: can the G1 shim passes fail? Record the generated
# distributions test twice, with the true adapter and with brms_g1_d()'s
# row_lpdf() value multiplied by 1 + 1e-6, and compare per assertion.
# Staged under dev/remeasure-rev3-stage/ (deleted afterwards).
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = "frmtmb", FRMTMB_PORT_ROOT = getwd())
suppressMessages({library(testthat); library(frmtmb)})
f <- "dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-distributions.R"
run <- function(tag, perturb) {
  st <- file.path("dev/remeasure-rev3-stage", tag)
  unlink(st, recursive = TRUE); dir.create(st, recursive = TRUE)
  file.copy(list.files("tests/testthat", "^helper.*[.]R$",
                       full.names = TRUE), st)
  ad <- readLines("dev/brmsport-bin2-adapters.R")
  if (perturb) {
    h <- grep("ll <- as.numeric(frmtmb::row_lpdf(fam, x, x, dp, av, NULL))",
              ad, fixed = TRUE)
    stopifnot(length(h) == 1L)
    ad[h] <- sub("NULL))", "NULL)) * (1 + 1e-6)", ad[h], fixed = TRUE)
    cat("perturbed:", ad[h], "\n")
  }
  writeLines(ad, file.path(st, "helper-brms-bin2-adapters.R"))
  file.copy(f, st)
  rec <- normalizePath(file.path(st, "rec.tsv"), winslash = "/",
                       mustWork = FALSE)
  Sys.setenv(FRMTMB_BRMSPORT_RECORD = rec)
  invisible(test_file(file.path(st, basename(f)), package = "frmtmb",
                      env = testthat::test_env("frmtmb"),
                      reporter = "silent", stop_on_failure = FALSE))
  x <- utils::read.delim(rec, header = FALSE, quote = "",
                         colClasses = "character")
  x[x[[1]] == "assert", c(3, 5)]
}
a <- run("true", FALSE); b <- run("perturbed", TRUE)
m <- merge(a, b, by = "V3", suffixes = c("_true", "_pert"))
led <- utils::read.delim("dev/brmsport-bin2-ledger.tsv", quote = "",
                         stringsAsFactors = FALSE)
m$kind <- led$kind[match(m$V3, led$id)]
m <- m[m$V5_true == "TRUE", ]
g1 <- grepl("^distributions:([0-9]|[0-9][0-9]|1[0-9][0-9]|20[0-9])$", m$V3)
cat("G1 rows holding on the true adapter:", sum(g1), "\n")
print(table(kind = m$kind[g1], holds_perturbed = m$V5_pert[g1]))
print(m[g1 & m$kind == "value", ], row.names = FALSE)
