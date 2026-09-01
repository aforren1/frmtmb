# Standalone runner for the grammar fuzz tier.
#
#   Rscript dev/fuzz/run-fuzz.R [size] [out_prefix]
#
# Writes <prefix>.json (machine readable) and prints a triaged summary.
# The test tier (tests/testthat/test-fuzz.R) runs the same plan; this
# script exists so a full run can be done without testthat's reporter
# swallowing progress on a ten-minute job.

args <- commandArgs(trailingOnly = TRUE)
size <- if (length(args) >= 1) as.integer(args[1]) else 300L
prefix <- if (length(args) >= 2) args[2] else "dev/fuzz-findings"

if (!dir.exists("R") || !file.exists("DESCRIPTION")) {
  stop("run this from the package root: Rscript dev/fuzz/run-fuzz.R")
}

Sys.setenv(NOT_CRAN = "true", FRMTMB_FUZZ = "true")
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-fuzz.R")
invisible(requireNamespace("brms", quietly = TRUE))   # attach cost, once

plan <- fuzz_plan(seed = 20260901L, size = size)
cat("specs:", nrow(plan),
    " pair-cover rows:", attr(plan, "n_cover_rows"),
    " feasible pairs:", attr(plan, "n_pairs"),
    " covered:", attr(plan, "n_pairs_covered"), "\n\n")

res <- fuzz_run(plan, progress = TRUE)
sm <- fuzz_summary(res)

cat("\n=== ", round(res$elapsed, 1), "s, ", res$n_specs, " specs, ",
    res$n_fits, " fits ===\n", sep = "")
print(sm$counts)
cat("\n")
for (f in sm$triaged) cat(fuzz_format(f))

# Machine-readable dump. No jsonlite dependency: the records are flat
# enough to serialize by hand, and dev/ must not add an Import.
esc <- function(s) {
  s <- gsub("\\\\", "\\\\\\\\", as.character(s))
  s <- gsub("\"", "\\\\\"", s)
  s <- gsub("\n", "\\\\n", s)
  gsub("\t", "\\\\t", s)
}
as_json <- function(x) {
  if (is.null(x)) return("null")
  if (is.list(x)) {
    if (!is.null(names(x))) {
      return(paste0("{", paste(sprintf("\"%s\":%s", esc(names(x)),
                                       vapply(x, as_json, "")),
                               collapse = ","), "}"))
    }
    return(paste0("[", paste(vapply(x, as_json, ""), collapse = ","), "]"))
  }
  if (length(x) != 1) {
    return(paste0("[", paste(vapply(x, as_json, ""), collapse = ","), "]"))
  }
  if (is.numeric(x) && is.finite(x)) return(format(x, digits = 12))
  paste0("\"", esc(x), "\"")
}
writeLines(as_json(list(
  seed = 20260901L, n_specs = res$n_specs, n_fits = res$n_fits,
  elapsed_sec = round(res$elapsed, 2),
  feasible_pairs = attr(plan, "n_pairs"),
  covered_pairs = attr(plan, "n_pairs_covered"),
  cover_rows = attr(plan, "n_cover_rows"),
  findings = unname(sm$triaged)
)), paste0(prefix, ".json"))
cat("\nwrote ", prefix, ".json\n", sep = "")
