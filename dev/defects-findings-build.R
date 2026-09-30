# Lane wt-defects: paste the generated blocks into
# dev/defects-findings.md, so that no count there is typed.
#   Rscript dev/defects-ledger-diff.R > dev/defects-log/ledger-diff.md
#   Rscript dev/defects-findings-build.R
p <- "dev/defects-findings.md"
doc <- readLines(p)
splice <- function(doc, tag, body) {
  b <- grep(sprintf("^<!-- BEGIN GENERATED: %s -->$", tag), doc)
  e <- grep("^<!-- END GENERATED -->$", doc)
  e <- e[e > b][1]
  stopifnot(length(b) == 1L, !is.na(e))
  c(doc[seq_len(b)], body, doc[e:length(doc)])
}
doc <- splice(doc, "ledger-diff", readLines("dev/defects-log/ledger-diff.md"))

# one line per suite log: files, and the summed counts from its RESULT
# lines; a file with no RESULT line is counted and named
sumlog <- function(f, label) {
  x <- readLines(f)
  r <- grep("^RESULT", x, value = TRUE)
  num <- function(k) sum(as.integer(sub(paste0(".* ", k, "= *([0-9]+).*"),
                                        "\\1", r)))
  hdr <- grep("^(SUITE|GATED|TIER)", x, value = TRUE)
  nores <- grep("^NO RESULT", x, value = TRUE)
  if (grepl("tier", f)) {
    n <- function(k) sum(as.integer(sub(paste0(".* ", k, "=([0-9]+).*"),
                                         "\\1", r)))
    return(sprintf("| %s | %s | %d | %d | %d | %d | - |", label, hdr,
                   n("pass"), n("fail"), n("err"), n("skip")))
  }
  sprintf("| %s | %s%s | %d | %d | %d | %d | %d |", label, hdr,
          if (length(nores)) paste0("; ", paste(nores, collapse = ", "))
          else "", num("pass"), num("fail"), num("err"), num("skip"),
          num("warn"))
}
logs <- c(
  "frmtmb, every file" = "dev/defects-log/suite-core-final.txt",
  "frmtmb.sample, every file" = "dev/defects-log/suite-sample-final.txt",
  "frmtmb, gated files" = "dev/defects-log/gated-core-final.txt",
  "frmtmb.sample, gated files" = "dev/defects-log/gated-sample-final.txt",
  "frmtmb.coupling" = "dev/defects-log/suite-ext-frmtmb.coupling.txt",
  "frmtmb.eam" = "dev/defects-log/suite-ext-frmtmb.eam.txt",
  "frmtmb.latent" = "dev/defects-log/suite-ext-frmtmb.latent.txt",
  "frmtmb.learn" = "dev/defects-log/suite-ext-frmtmb.learn.txt",
  "frmtmb.ode" = "dev/defects-log/suite-ext-frmtmb.ode.txt",
  "frmtmb.spline" = "dev/defects-log/suite-ext-frmtmb.spline.txt",
  "ported brms tier, verdicts asserted" = "dev/defects-log/tier-final.txt")
body <- c("| suite | files | pass | fail | error | skip | warn |",
          "|---|---|---|---|---|---|---|")
for (k in names(logs)) {
  if (file.exists(logs[[k]])) body <- c(body, sumlog(logs[[k]], k))
}
doc <- splice(doc, "suites", body)
writeLines(doc, p)
