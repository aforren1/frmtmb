# Driver: one Rscript process per vignette so a crash inside a fit costs
# one vignette, not the whole audit.  Rscript run-all.R [cap] [vignettes...]
args <- commandArgs(trailingOnly = TRUE)
CAP <- if (length(args)) args[1] else "120"
MODE <- if (length(args) > 1) args[2] else "raw"
VIGS <- if (length(args) > 2) args[-(1:2)] else c(
  "brms_overview", "brms_multilevel", "brms_distreg", "brms_nonlinear",
  "brms_phylogenetics", "brms_monotonic", "brms_multivariate",
  "brms_missings", "brms_customfamilies"
)
HERE <- "C:/Users/adf44/source/r/frmtmb-wt-audit/dev/brms-port"
RSCRIPT <- file.path(R.home("bin"), "Rscript")
for (v in VIGS) {
  cat("=====", v, "\n")
  st <- system2(RSCRIPT, c(file.path(HERE, "run-vignette.R"), v, CAP, MODE),
                stdout = FALSE, stderr = FALSE)
  cat("  exit:", st, "\n")
}
