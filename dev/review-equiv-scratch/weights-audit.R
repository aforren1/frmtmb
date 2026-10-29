# Reviewer: audit two reach patterns. Run as
#   Rscript dev/review-equiv-scratch/weights-audit.R dev/review-equiv-scratch/wa
# (the lane script writes its weights.tsv into that scratch directory).
dir.create(commandArgs(trailingOnly = TRUE)[1], showWarnings = FALSE)
sink(tempfile()); source("dev/equiv-weights.R"); sink()
hitdocs <- function(p) names(docs)[vapply(docs, function(t)
  grepl(p, t, perl = TRUE), NA)]
cat("bar_crossing pattern:", pat[["bar_crossing"]], "\n")
for (nm in hitdocs(pat[["bar_crossing"]])) {
  m <- regmatches(docs[[nm]], gregexpr(pat[["bar_crossing"]], docs[[nm]],
                                        perl = TRUE))[[1]]
  cat("  ", nm, ":", paste(unique(m), collapse = " ; "), "\n")
}
cat("docs with '* (' or ': (' crossing a bar:",
    length(hitdocs("[*:]\\s*\\([^()]*\\|")), "\n")
# prior: documents whose only match is a prior-reading call
pr <- hitdocs(pat[["prior"]])
setting <- "set_prior\\(|(?<![A-Za-z_])prior\\s*=|(?<![A-Za-z_.])prior\\("
pr_set <- hitdocs(setting)
cat("prior: pattern hits", length(pr), "; hits that SET a prior",
    length(intersect(pr, pr_set)), "\n")
cat("  hit only by get_prior/sample_prior/prior_summary/etc.:\n")
for (nm in setdiff(pr, pr_set)) cat("   ", nm, "\n")
# frm_curve: documents credited only through conditional_effects()
fc <- hitdocs(pat[["frm_curve"]])
cat("frm_curve hits", length(fc), "; via conditional_effects only",
    length(setdiff(fc, hitdocs("frm_curve\\("))), "\n")
