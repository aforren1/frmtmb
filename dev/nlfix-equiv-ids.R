# The pairs of dev/equiv-data/sample.tsv this lane's fix can reach: every spec
# with a nonlinear body, an ordinal family or a compat row the lane changes.
# Writes dev/nlfix-equiv-ids.txt.
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R")) sys.source(f, envir = E)
smp <- read.delim("dev/equiv-data/sample.tsv")
keep <- character(0)
pat <- "nl = TRUE|nlf[(]|cumulative|sratio|cratio|acat|autocor *=|beta_binomial"
for (i in seq_len(nrow(smp))) {
  key <- paste(smp$feature_a[i], smp$feature_b[i], sep = " x ")
  mk <- E$eq_specs[[key]]
  if (is.null(mk)) next
  sp <- tryCatch(mk(), error = function(e) NULL)
  if (is.null(sp)) next
  txt <- paste(c(sp$bf, sp$bf_frm, sp$bf_brms, sp$family, sp$family_frm),
               collapse = " ")
  if (grepl(pat, txt) || "nl" %in% c(smp$feature_a[i], smp$feature_b[i]) ||
        grepl("beta_binomial", key, fixed = TRUE)) {
    keep <- c(keep, smp$id[i])
  }
}
writeLines(keep, "dev/nlfix-equiv-ids.txt")
cat(length(keep), "pairs:", keep, "\n")
