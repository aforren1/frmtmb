## The core suite, ONE TEST FILE PER PROCESS, joined against
## dev/suite-baseline.tsv so a count that FELL is reported as a question
## rather than hidden behind a green line.
##
## Usage: Rscript dev/argspell-suite.R [pkgdir]
## Writes dev/argspell-suite.tsv and prints the comparison.
args <- commandArgs(trailingOnly = TRUE)
pkgdir <- if (length(args)) args[1] else "."
pkg <- if (identical(pkgdir, ".")) "frmtmb" else basename(pkgdir)
root <- normalizePath(".")
tdir <- file.path(root, pkgdir, "tests", "testthat")
files <- sort(basename(list.files(tdir, pattern = "^test-.*[.]R$")))
cat("files:", length(files), "\n")
rows <- list()
rscript <- file.path(R.home("bin"), "Rscript")
t0 <- Sys.time()
for (f in files) {
  filt <- paste0("^", sub("^test-", "", sub("[.]R$", "", f)), "$")
  out <- suppressWarnings(system2(
    rscript, c(file.path(root, "dev", "argspell-run1.R"),
               shQuote(pkgdir), shQuote(filt)),
    stdout = TRUE, stderr = TRUE))
  line <- grep("^ARGSPELL ", out, value = TRUE)
  if (!length(line)) {
    rows[[length(rows) + 1L]] <- data.frame(
      file = f, pass = NA_integer_, fail = NA_integer_,
      error = NA_integer_, skip = NA_integer_, crashed = TRUE)
    cat(sprintf("%-46s CRASHED\n", f))
    writeLines(tail(out, 15))
    next
  }
  p <- as.integer(sub(".* pass ([0-9]+) .*", "\\1", line))
  fa <- as.integer(sub(".* fail ([0-9]+) .*", "\\1", line))
  er <- as.integer(sub(".* error ([0-9]+) .*", "\\1", line))
  sk <- as.integer(sub(".* skip ([0-9]+) .*", "\\1", line))
  rows[[length(rows) + 1L]] <- data.frame(file = f, pass = p, fail = fa,
                                          error = er, skip = sk,
                                          crashed = FALSE)
  cat(sprintf("%-46s pass %5d fail %3d error %3d skip %3d\n",
              f, p, fa, er, sk))
  if (fa > 0L || er > 0L) {
    writeLines(grep("Failure [(]|Error [(]|^Error", out, value = TRUE))
  }
}
tab <- do.call(rbind, rows)
utils::write.table(tab, file.path("dev", "argspell-suite.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat("\n---- GENERATED COUNTS (paste verbatim) ----\n")
cat(sprintf("package        : %s\n", pkg))
cat(sprintf("files run      : %d\n", nrow(tab)))
cat(sprintf("crashed        : %d\n", sum(tab$crashed)))
cat(sprintf("pass %d  fail %d  error %d  skip %d  in %.0f s\n",
            sum(tab$pass, na.rm = TRUE), sum(tab$fail, na.rm = TRUE),
            sum(tab$error, na.rm = TRUE), sum(tab$skip, na.rm = TRUE),
            secs))
bl <- utils::read.delim("dev/suite-baseline.tsv", header = FALSE,
                        col.names = c("pkg", "file", "pass", "skip"))
bl <- bl[bl$pkg == pkg, ]
m <- merge(tab, bl, by = "file", all.x = TRUE, suffixes = c("", "_base"))
drop <- m[!is.na(m$pass_base) & !is.na(m$pass) & m$pass < m$pass_base, ]
miss <- bl$file[!bl$file %in% tab$file]
cat(sprintf("files below their 0.57.0 baseline count: %d\n", nrow(drop)))
if (nrow(drop)) {
  print(drop[, c("file", "pass", "pass_base", "skip", "skip_base")],
        row.names = FALSE)
}
cat(sprintf("baseline files with no run: %d %s\n", length(miss),
            paste(miss, collapse = " ")))
cat("---- END GENERATED COUNTS ----\n")
