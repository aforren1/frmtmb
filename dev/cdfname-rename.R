# Rename a family's CDF slot from `lcdf` to `cdf` across the tree.
#
# Why a script: the slot is shared by every family file, and lanes of the
# same round add families that still fill the old name. The consolidation
# reruns this on the merged tree so those lines follow. It is idempotent:
# a second run on a renamed tree changes nothing.
#
#   Rscript dev/cdfname-rename.R [root] [--dry-run]
#
# root defaults to the working directory, which must be the package root
# (it holds DESCRIPTION and extensions/). --dry-run lists what would change
# and writes nothing. Run roxygen2::roxygenise() on core and on each
# changed extension afterwards; the script edits man/*.Rd too, so a
# roxygenise that changes an Rd shows a line the script and the roxygen
# source disagree on.
#
# What it leaves alone, and why:
# - NEWS.md files: history names the slot as it was.
# - tests/testthat/test-cdf-slot-name.R and any line that holds
#   `old_cdf_slot`: the refusal of the old name, and its test, must spell it.
# - Names with a prefix, `poisson_lcdf`, `ddm_rt_lcdf2`, `nlcdf.lower`:
#   Stan's functions, brms's, TMB's oneStepPredict() columns, and
#   frmtmb.eam's helpers that do return a log CDF. Their names are right.
# - dev/ scripts other than the three the ported-suite ledger runs: they
#   are records of past runs against the builds of their day.

args <- commandArgs(trailingOnly = TRUE)
dry <- "--dry-run" %in% args
args <- setdiff(args, "--dry-run")
root <- if (length(args)) args[[1L]] else "."
stopifnot(file.exists(file.path(root, "DESCRIPTION")),
          dir.exists(file.path(root, "extensions")))

pkgs <- c(".", file.path("extensions",
                         list.files(file.path(root, "extensions"))))
pkgs <- pkgs[file.exists(file.path(root, pkgs, "DESCRIPTION"))]
listed <- function(dir, pattern) {
  d <- file.path(root, dir)
  if (!dir.exists(d)) return(character(0))
  file.path(dir, list.files(d, pattern = pattern, recursive = TRUE))
}
files <- unlist(lapply(pkgs, function(p) c(
  listed(file.path(p, "R"), "[.]R$"),
  listed(file.path(p, "tests", "testthat"), "[.]R$"),
  listed(file.path(p, "man"), "[.]Rd$"),
  listed(file.path(p, "vignettes"), "[.]Rmd$"),
  listed(file.path(p, "inst"), "[.](R|Rmd)$"))))
files <- c(files, "SPEC.md", "dev/brmsport-bin2-adapters.R",
           "dev/brmsport-bin2-survey.R", "dev/brmsport-bin2-verdicts-2p.R")
files <- sub("^[.]/", "", unique(files))
files <- files[file.exists(file.path(root, files))]
files <- files[basename(files) != "test-cdf-slot-name.R"]

# a name with no prefix of letters, digits, `_` or `.`, and no suffix of
# letters or digits: `lcdf`, `lcdf_rows`, `$lcdf`, `"lcdf"`, `\code{lcdf}`
bare <- "(?<![A-Za-z0-9_.])lcdf(?![A-Za-z0-9])"

rename_line <- function(x, r_code) {
  x <- gsub(paste0("\\b([Aa])n ([`]|\\\\code\\{)?", bare), "\\1 \\2cdf", x,
            perl = TRUE)
  # `$` partial-matches: on a hand-built list with no `cdf` element
  # `fam$cdf` would return `cdf_log`, so package code reads it exactly
  if (r_code && !grepl("^\\s*#", x)) {
    x <- gsub("\\$lcdf(?![A-Za-z0-9_])", "[[\"cdf\"]]", x, perl = TRUE)
  }
  x <- gsub("\\bfam_lcdf\\b", "fam_cdf", x, perl = TRUE)
  x <- gsub("\\bddm_cont_lcdf\\b", "ddm_cont_cdf", x, perl = TRUE)
  gsub(bare, "cdf", x, perl = TRUE)
}

total <- 0L
for (f in files) {
  path <- file.path(root, f)
  txt <- readChar(path, file.size(path), useBytes = TRUE)
  # split on LF only, so a CRLF file keeps its CR on each line
  lines <- strsplit(txt, "\n", fixed = TRUE, perl = FALSE)[[1L]]
  r_code <- grepl("(^|/)R/[^/]+[.]R$", f)
  keep <- grepl("old_cdf_slot", lines, fixed = TRUE)
  new <- lines
  new[!keep] <- vapply(lines[!keep], rename_line, "", r_code = r_code,
                       USE.NAMES = FALSE)
  hit <- which(new != lines)
  if (!length(hit)) next
  total <- total + length(hit)
  cat(sprintf("%-60s %3d line(s)\n", f, length(hit)))
  if (dry) next
  out <- paste(new, collapse = "\n")
  # strsplit drops a trailing empty field, so restore the final newline
  if (endsWith(txt, "\n")) out <- paste0(out, "\n")
  con <- file(path, "wb")
  writeChar(out, con, eos = NULL, useBytes = TRUE)
  close(con)
}
cat(sprintf("RENAME %s: %d line(s) in %d candidate file(s)\n",
            if (dry) "dry run" else "done", total, length(files)))
