# Reviewer, lane wt-priorform: what brms and base frmtmb do with a
# group-level term written twice, which brms accepts and the lane
# refuses. brms: count of sd parameters in the Stan data; frmtmb base:
# count of random-effect blocks in the frame. Seed 20260916.
#   Rscript dev/priorform-rev-twins.R brms|ref
mode <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (mode == "ref") "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
set.seed(20260916)
n <- 120
d <- data.frame(g = factor(rep(1:12, 10)), h = factor(rep(1:10, each = 12)),
                x = rnorm(n))
d$y <- rnorm(n) + rnorm(12)[d$g]
fs <- c("y ~ x + (1 | g) + (1 | g)", "y ~ x + (1 | g:h) + (1 | h:g)",
        "y ~ x + (x | g) + (x | g)")
if (mode == "brms") {
  suppressMessages(library(brms))
  for (f in fs) {
    sd <- standata(as.formula(f), data = d)
    cat(f, ": M_ =", paste(names(sd)[grepl("^M_", names(sd))], collapse = " "),
        "; N_ =", paste(names(sd)[grepl("^N_[0-9]", names(sd))], collapse = " "), "\n")
  }
} else {
  suppressMessages(library(frmtmb))
  for (f in fs) {
    fit <- suppressWarnings(frm(bf(as.formula(f)) + gaussian(), data = d))
    cat(f, ": blocks =", length(fit$frame$re_blocks), "; VarCorr names =",
        paste(names(VarCorr(fit)), collapse = ", "),
        "; logLik", sprintf("%.6f", as.numeric(logLik(fit))), "\n")
  }
  for (f in c("y ~ x + (1 | g)", "y ~ x + (1 | g:h)", "y ~ x + (x | g)")) {
    fit <- frm(bf(as.formula(f)) + gaussian(), data = d)
    cat(f, ": blocks =", length(fit$frame$re_blocks),
        "; logLik", sprintf("%.6f", as.numeric(logLik(fit))), "\n")
  }
}
