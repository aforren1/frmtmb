# Reviewer: arm B's third-party identity, read with the padded "NA"
# strings stripped. The lane's writer pads numeric fields, so read.delim
# sees "           NA" and hands back a character column: a denominator
# read straight off that file counts 15 reference fits where 6 ran.

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
B <- utils::read.delim(paste0(D, "latent-2p3-hmm-B.tsv"),
                       stringsAsFactors = FALSE, strip.white = TRUE)
num <- function(v) suppressWarnings(as.numeric(trimws(as.character(v))))
for (nm in c("ll_ref", "ll_rel_gap", "est_max_gap", "tpm_max_gap", "s_ref")) {
  B[[nm]] <- num(B[[nm]])
}
idx <- !is.na(B$ll_ref)
cat("replicates:", nrow(B), " with a hmmTMB fit:", sum(idx),
    " seeds:", paste(B$seed[idx], collapse = " "), "\n")
cat("ll_rel_gap  :", format(range(abs(B$ll_rel_gap[idx])), digits = 3),
    " claim 6.1e-12 to 6.9e-11\n")
cat("est_max_gap :", format(range(B$est_max_gap[idx]), digits = 3),
    " claim 5.2e-06 to 1.3e-04\n")
cat("tpm_max_gap :", format(range(B$tpm_max_gap[idx]), digits = 3),
    " claim 0.039 to 0.180\n")
cat("mean se on a state mean:",
    format(range(c(B$se_mu1[idx], B$se_mu2[idx], B$se_mu3[idx])),
           digits = 3), " claim 0.008 to 0.009\n")
cat("1.3e-04 as a fraction of the smallest of those:",
    format(max(B$est_max_gap[idx]) /
             min(c(B$se_mu1[idx], B$se_mu2[idx], B$se_mu3[idx])),
           digits = 3), " claim 1.5e-02\n")
cat("identity_ok is the LABEL match, not the reference match:",
    "sum =", sum(B$identity_ok), "of", nrow(B), "\n")
