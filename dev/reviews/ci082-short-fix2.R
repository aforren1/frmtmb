# Short check of lane ci082's fix round, FIX 2: the two phrases the
# nonzero-code branch of test-se-limit.R now requires, against each of
# the three endings conv_headroom_clause() writes after a code-8 stop.
#   Rscript dev/reviews/ci082-short-fix2.R
.libPaths(c("C:/Users/adf44/source/r/wt-ci082-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
ns <- asNamespace("frmtmb")
for (room in c(9.9e-7, 0.0657, Inf)) {
  w <- paste0("Optimizer did not report convergence: false convergence (8)",
              ns$conv_headroom_clause(room, 1e-3))
  ok <- grepl("false convergence (8)", w, fixed = TRUE) &&
    grepl("so the estimates are at the maximum", w, fixed = TRUE)
  cat(sprintf("RESULT room %g phrases-pass %s | %s\n", room, ok,
              substr(w, 61, 200)))
}
