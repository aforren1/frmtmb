# Review of lane gpgap, item 8: the headline under each defensible verdict for
# the exact gp() pairs. A scratch tree holds the survey's tables and the
# lane's per-pair results; dev/equiv-summary.R runs unchanged in it, and only
# the adjudication differs. Nothing in dev/equiv-data/ is touched.
#
#   Rscript dev/review-gpgap-rates.R
root <- normalizePath("dev/review-gpgap-rates", winslash = "/", mustWork = FALSE)
unlink(root, recursive = TRUE)
dir.create(file.path(root, "dev"), recursive = TRUE)
file.copy("dev/equiv-data", file.path(root, "dev"), recursive = TRUE)
file.copy("dev/equiv-summary.R", file.path(root, "dev"))
base <- "C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out"
src_lane <- "dev/gpgap-equiv-out"
ids <- c("p022", "p023", "p094", "p109")
adj0 <- read.delim("dev/equiv-data/adjudication.tsv", stringsAsFactors = FALSE,
                   colClasses = "character", na.strings = character())
adj_row <- function(id, verdict, agree, level, reason) {
  data.frame(id = id, verdict = verdict, agree = agree, level = level,
             proposal = "", reason = reason, stringsAsFactors = FALSE)
}
drop4 <- adj0[!adj0$id %in% ids, ]
dbd <- "nugget line differs by decision; rest is one model through the edited line"
variants <- list(
  # the record: L1 unexplained, cannot measure
  baseline = list(out = base, adj = adj0),
  # the lane's numbers: rows dropped, verdict from the edited program
  lane_equivalent = list(out = src_lane, adj = drop4),
  # the standard's verdict for every program that differs: p002's treatment,
  # with p136 keeping only what brms's own program supports (L4, weak)
  divergence_4_p136_L4 = list(out = src_lane, adj = rbind(
    drop4,
    do.call(rbind, lapply(ids, adj_row, verdict = "divergence by decision",
                          agree = "FALSE", level = "L1", reason = dbd)),
    adj_row("p136", "equivalent", "TRUE", "L4", "L4 only; L1 is the edited line"))),
  # all five divergences, p136 included
  divergence_5 = list(out = src_lane, adj = rbind(
    drop4,
    do.call(rbind, lapply(c(ids, "p136"), adj_row, verdict = "divergence by decision",
                          agree = "FALSE", level = "L1", reason = dbd)))),
  # closed-form maxima at brms's own nugget stand as L4 on three pairs;
  # p109 and p136 have no closed form
  closed_form_3 = list(out = src_lane, adj = rbind(
    drop4,
    do.call(rbind, lapply(c("p022", "p023", "p094"), adj_row, verdict = "equivalent",
                          agree = "TRUE", level = "L4",
                          reason = "maxima at brms's nugget by closed form")),
    adj_row("p109", "divergence by decision", "FALSE", "L1", dbd),
    adj_row("p136", "equivalent", "TRUE", "L4", "L4 only; L1 is the edited line"))))
rs <- "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
wd <- getwd()
for (nm in names(variants)) {
  v <- variants[[nm]]
  o <- file.path(root, paste0("out-", nm))
  dir.create(o)
  fs <- Sys.glob(file.path(v$out, "p*.tsv"))
  file.copy(fs, o)
  utils::write.table(v$adj, file.path(root, "dev/equiv-data/adjudication.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  setwd(root)
  system2(rs, c("dev/equiv-summary.R", paste0("out-", nm)), stdout = FALSE, stderr = FALSE)
  setwd(wd)
  s <- readLines(file.path(o, "summary.md"), warn = FALSE)
  g <- function(p) {
    x <- grep(p, s, value = TRUE)[1]
    sub("^.*?: ", "", sub(";.*$", "", sub("^- \\**", "", x)))
  }
  r <- read.delim(file.path(o, "results.tsv"), stringsAsFactors = FALSE)
  w <- r[r$id %in% c(ids, "p136"), c("id", "verdict", "agree", "top_level", "weak_evidence")]
  cat(sprintf("RATE %-22s agree %s | no defect %s | weak out %s | measured only %s\n", nm,
              g("the headline"), g("^- no defect"), g("sensitivity, weak"), g("among measured")))
  cat("    ", paste(sprintf("%s=%s(%s%s)", w$id, substr(w$verdict, 1, 5), w$top_level,
                            ifelse(w$weak_evidence, ",weak", "")), collapse = " "), "\n")
}
