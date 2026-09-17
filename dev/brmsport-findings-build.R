# Paste the generated blocks into dev/brmsport-findings.md, so no count
# in it is typed: the ledger summary, the tier run, and the number of
# generated test lines past 80 columns.
#
#   Rscript dev/brmsport-findings-build.R
p <- "dev/brmsport-findings.md"
doc <- readLines(p)

splice <- function(doc, tag, body) {
  b <- grep(sprintf("^<!-- BEGIN GENERATED: %s -->$", tag), doc)
  e <- grep("^<!-- END GENERATED -->$", doc)
  e <- e[e > b][1]
  stopifnot(length(b) == 1L, !is.na(e))
  c(doc[seq_len(b)], body, doc[e:length(doc)])
}

ledger <- readLines("dev/brmsport-log/ledger-summary.md")
ledger <- ledger[!grepl("^<!-- generated", ledger)]
doc <- splice(doc, "dev/brmsport-ledger.R", ledger)

tier <- readLines("dev/brmsport-log/tier-final.txt")
res <- grep("^RESULT|^TIER|^NO RESULT", tier, value = TRUE)
stopifnot(length(grep("^TIER ran", res)) == 1L)
n <- as.integer(sub("^RESULT .* pass=([0-9]+) .*$", "\\1",
                    grep("^RESULT", res, value = TRUE)))
bad <- grep("fail=[1-9]|err=[1-9]|skip=[1-9]", res, value = TRUE)
gen <- c(Sys.glob("tests/testthat/test-brms-suite-*.R"),
         Sys.glob("extensions/frmtmb.sample/tests/testthat/test-brms-suite-*.R"))
long <- sum(vapply(gen, function(f) sum(nchar(readLines(f)) > 80), 1L))
doc <- splice(doc, "tier", c(
  "```", res, "```", "",
  sprintf(paste("%d expectations over %d files; %d files with a failure,",
                "error or skip. Generated lines past 80 columns: %d."),
          sum(n), length(n), length(bad), long)))
writeLines(doc, p)
cat("spliced; tier expectations", sum(n), "files", length(n), "bad",
    length(bad), "long lines", long, "\n")
