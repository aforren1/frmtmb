# Driver: one Rscript process per vignette so a crash inside a fit costs
# one vignette, not the whole audit.
#
#   Rscript run-all.R [cap] [mode] [vignettes...]
#
# PORT_LIB and PORT_OUT (see env.R) pass through to each child. The
# vignettes run at the same time, because each process writes only its
# own files and the machine of 2026-09-28 has the memory for nine.
HERE <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
args <- commandArgs(trailingOnly = TRUE)
CAP <- if (length(args)) args[1] else "120"
MODE <- if (length(args) > 1) args[2] else "raw"
VIGS <- if (length(args) > 2) args[-(1:2)] else c(
  "brms_overview", "brms_multilevel", "brms_distreg", "brms_nonlinear",
  "brms_phylogenetics", "brms_monotonic", "brms_multivariate",
  "brms_missings", "brms_customfamilies"
)
RSCRIPT <- file.path(R.home("bin"), "Rscript")
cl <- parallel::makePSOCKcluster(length(VIGS))
st <- parallel::parLapply(cl, VIGS, function(v, rs, script, cap, mode) {
  system2(rs, c(script, v, cap, mode), stdout = FALSE, stderr = FALSE)
}, rs = RSCRIPT, script = file.path(HERE, "run-vignette.R"), cap = CAP,
mode = MODE)
parallel::stopCluster(cl)
for (i in seq_along(VIGS)) cat("=====", VIGS[i], " exit:", st[[i]], "\n")
