# Summarise every run of this lane into the marked block of
# dev/adefects-findings.md, and STOP if a file produced no RESULT line
# or produced one with zero expectations.
#
#   Rscript dev/adefects-runs.R
#
# dev/lane-rules.md: emit counts from the summariser into a marked block
# and paste that block verbatim, never type them. A runner that sums
# `failed` and not `error` prints a clean line for a file that aborted
# halfway, so `err` and the missing-line case are both counted here, and
# a file that reports zero expectations and zero skips is called out.
# Both happened on this lane: four files on one pass, because other R
# work was running beside the suite, and 43 on another, because the
# launcher did not export NOT_CRAN.
logs <- c(
  "core suite" = "dev/adefects-log/suite-core.txt",
  "frmtmb.coupling" = "dev/adefects-log/suite-coupling.txt",
  "frmtmb.eam" = "dev/adefects-log/suite-eam.txt",
  "frmtmb.latent" = "dev/adefects-log/suite-latent.txt",
  "frmtmb.learn" = "dev/adefects-log/suite-learn.txt",
  "frmtmb.ode" = "dev/adefects-log/suite-ode.txt",
  "frmtmb.sample" = "dev/adefects-log/suite-sample.txt",
  "frmtmb.spline" = "dev/adefects-log/suite-spline.txt",
  "gated tier" = "dev/adefects-log/gated.txt",
  "ported brms tier" = "dev/adefects-log/port-tier.txt")

num <- function(ln, key) {
  m <- regmatches(ln, regexpr(paste0(key, "=[0-9]+"), ln))
  if (!length(m)) return(NA_integer_)
  as.integer(sub(".*=", "", m))
}

rows <- list()
bad <- character(0)
for (nm in names(logs)) {
  f <- logs[[nm]]
  if (!file.exists(f)) {
    bad <- c(bad, paste("missing log:", f))
    next
  }
  ln <- readLines(f, warn = FALSE)
  res <- grep("^RESULT ", ln, value = TRUE)
  noout <- grep("NOOUTPUT|LOADERROR", res, value = TRUE)
  res <- setdiff(res, noout)
  p <- vapply(res, num, 1L, key = "pass")
  fl <- vapply(res, num, 1L, key = "fail")
  er <- vapply(res, num, 1L, key = "err")
  sk <- vapply(res, num, 1L, key = "skip")
  zero <- res[!is.na(p) & p == 0L & sk == 0L]
  # The ONE exemption, and it is only sound because another run covers
  # them: the generated brms-suite files skip WHOLE when
  # FRMTMB_BRMS_FIT_TESTS is unset, which is how the ungated suite sees
  # them, and dev/brmsport-tier.sh runs the same files with it set. The
  # block after the loop asserts that the ported-tier log really does
  # report expectations for them, so the exemption cannot hide a file
  # that runs nowhere.
  if (nm != "ported brms tier") {
    zero <- zero[!grepl("test-brms-suite-", zero, fixed = TRUE)]
  }
  if (length(noout)) bad <- c(bad, paste(nm, "no RESULT line:", noout))
  if (length(zero)) bad <- c(bad, paste(nm, "zero expectations:", zero))
  if (sum(fl, na.rm = TRUE) + sum(er, na.rm = TRUE) > 0) {
    bad <- c(bad, paste(nm, "failures or errors"))
  }
  rows[[length(rows) + 1L]] <- data.frame(
    what = nm, files = length(res), pass = sum(p, na.rm = TRUE),
    fail = sum(fl, na.rm = TRUE), err = sum(er, na.rm = TRUE),
    skip = sum(sk, na.rm = TRUE), stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows)

# The exemption's other half: every file the ungated suites skipped
# whole must have run, with expectations, in the ported tier.
pt <- readLines("dev/adefects-log/port-tier.txt", warn = FALSE)
pt <- grep("^RESULT ", pt, value = TRUE)
pt_files <- sub("^RESULT [^ ]+ ([^ ]+) .*$", "\\1", pt)
pt_pass <- vapply(pt, num, 1L, key = "pass")
core_gated <- list.files("tests/testthat", "^test-brms-suite-.*[.]R$")
miss <- setdiff(core_gated, pt_files[pt_pass > 0L])
if (length(miss)) {
  bad <- c(bad, paste("gated file with no expectations anywhere:", miss))
}

blk <- c("", "| run | files | pass | fail | error | skip |",
         "|---|---|---|---|---|---|")
for (i in seq_len(nrow(tab))) {
  blk <- c(blk, sprintf("| %s | %d | %d | %d | %d | %d |", tab$what[i],
                        tab$files[i], tab$pass[i], tab$fail[i],
                        tab$err[i], tab$skip[i]))
}
blk <- c(blk,
         sprintf(paste0("| **total** | **%d** | **%d** | **%d** | ",
                        "**%d** | **%d** |"),
                 sum(tab$files), sum(tab$pass), sum(tab$fail),
                 sum(tab$err), sum(tab$skip)), "")

f <- "dev/adefects-findings.md"
ln <- readLines(f, warn = FALSE)
b <- grep("^<!-- BEGIN GENERATED: dev/adefects-runs.R -->$", ln)
e <- grep("^<!-- END GENERATED -->$", ln)
stopifnot(length(b) == 1L)
e <- e[e > b][1L]
writeLines(c(ln[seq_len(b)], blk, ln[e:length(ln)]), f)
cat(blk, sep = "\n")
if (length(bad)) {
  cat("\nPROBLEMS\n")
  cat(bad, sep = "\n")
  stop(length(bad), " problem(s)")
}
cat("\nevery run clean\n")
