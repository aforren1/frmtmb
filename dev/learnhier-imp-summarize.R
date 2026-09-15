# Lane `learnhier` (item 2.2): the importance correction's success rate
# at 100 learners by 200 trials.
#
#   Rscript dev/learnhier-imp-summarize.R <dir>
#
# THE COUNT IS `imp_stalled()`'s AND NOT AN EYE'S. Three outcomes:
# `settled` (the iteration reached its own fixed point inside the round
# cap), `stalled` (the cap bit and the rounds are taking the same step
# every time, which is what a collapsed variance component does to the
# fixed-point map) and `slow` (the cap bit and the steps are still
# shrinking). A run that the correction refused outright is `refused`.
# Only `settled` is a success.
#
# EVERY RUN'S OWN SPREAD IS PRINTED BESIDE ITS VERDICT. The threshold
# `imp_stall_tol = 1e-2` sits in a gap that `dev/reviews/2026-09-08-
# debts.md` measured at a factor of 3.2 over 184 fits, not the three
# orders of magnitude the roxygen claimed: the widest stalled spread it
# found was 6.51e-03 and the narrowest still-moving one 2.10e-02, which
# are margins of 1.54 and 2.10. A rate near either edge would be a
# finding about the THRESHOLD as much as about the model, so the table
# says where this design's spreads sit relative to both.
source("dev/learnhier-env.R")
DIR <- commandArgs(trailingOnly = TRUE)[[1L]]

fs <- list.files(DIR, pattern = "^imp-.*[.]rds$", full.names = TRUE)
R <- lapply(fs, readRDS)
cat("records ", length(R), " from ", DIR, "\n", sep = "")
if (!length(R)) quit(save = "no")

g <- function(f, d = NA_real_) {
  vapply(R, function(r) {
    v <- f(r)
    if (is.null(v) || !length(v)) d else as.numeric(v)[[1L]]
  }, numeric(1))
}
gs <- function(f, d = NA_character_) {
  vapply(R, function(r) {
    v <- f(r)
    if (is.null(v) || !length(v)) d else as.character(v)[[1L]]
  }, character(1))
}

seeds <- g(function(r) r$seed)
cat("distinct seeds ", length(unique(seeds)), " of ", length(seeds),
    " files; draws ", paste(sort(unique(g(function(r) r$draws))),
                            collapse = ","), "\n", sep = "")
out <- gs(function(r) r$outcome)
cat("\noutcome, counted with frmtmb:::imp_stalled():\n")
print(table(out))
n <- length(out)
p <- mean(out == "settled")
cat("settled ", sum(out == "settled"), " of ", n, " = ",
    round(100 * p, 1), " percent, binomial se ",
    round(100 * sqrt(p * (1 - p) / n), 1), " points\n", sep = "")
if (p == 1) {
  cat("a rate of 1 in ", n, " runs still admits a failure rate up to ",
      round(100 * (1 - 0.05^(1 / n)), 1),
      " percent at 95 percent confidence (the rule of three)\n", sep = "")
}

sp <- g(function(r) r$spread)
cat("\nround-to-round spread (max - min) / mean, which is the statistic",
    " imp_stalled() thresholds at 1e-2:\n", sep = "")
cat("  min ", format(min(sp, na.rm = TRUE), digits = 3), ", median ",
    format(stats::median(sp, na.rm = TRUE), digits = 3), ", max ",
    format(max(sp, na.rm = TRUE), digits = 3), "\n", sep = "")
cat("  smallest spread here is ",
    format(min(sp, na.rm = TRUE) / 1e-2, digits = 3),
    "x the threshold, and ",
    format(min(sp, na.rm = TRUE) / 2.10e-2, digits = 3),
    "x the narrowest still-moving spread the 184-fit review measured\n",
    sep = "")

cat("\nrounds used (cap is frmtmb_control(importance_rounds = 5)):\n")
print(table(g(function(r) r$rounds)))
cat("\ndiagnostics: min ESS per draw, over the worst group of each run\n")
es <- g(function(r) r$ess_min)
cat("  min ", round(min(es, na.rm = TRUE), 3), ", median ",
    round(stats::median(es, na.rm = TRUE), 3), ", max ",
    round(max(es, na.rm = TRUE), 3),
    "; below the 0.25 floor on ", sum(es < 0.25, na.rm = TRUE), " of ",
    n, "\n", sep = "")
mc <- g(function(r) r$mcse)
cat("  Monte Carlo se of the corrected log-likelihood: median ",
    round(stats::median(mc, na.rm = TRUE), 3), ", max ",
    round(max(mc, na.rm = TRUE), 3), "\n", sep = "")
cat("\nseconds: Laplace median ",
    round(stats::median(g(function(r) r$lap_s)), 1), ", correction median ",
    round(stats::median(g(function(r) r$imp_s)), 1), ", ratio ",
    round(stats::median(g(function(r) r$imp_s)) /
            stats::median(g(function(r) r$lap_s)), 1), "\n", sep = "")

cat("\nwhat the correction MOVED, per component, Laplace to corrected\n")
keys <- names(R[[1L]]$lap)
tab <- data.frame(component = keys, stringsAsFactors = FALSE)
truth <- c(stats::qlogis(0.35), log(3), 0.5, 0.3, 0.5)
for (k in seq_along(keys)) {
  la <- g(function(r) r$lap[[k]])
  im <- g(function(r) r$imp[[k]])
  tab$truth[k] <- truth[[k]]
  tab$laplace[k] <- mean(la, na.rm = TRUE)
  tab$corrected[k] <- mean(im, na.rm = TRUE)
  tab$shift[k] <- mean(im - la, na.rm = TRUE)
  tab$shift_mcse[k] <- stats::sd(im - la, na.rm = TRUE) /
    sqrt(sum(is.finite(im - la)))
  tab$lap_err[k] <- mean(la - truth[[k]], na.rm = TRUE)
  tab$cor_err[k] <- mean(im - truth[[k]], na.rm = TRUE)
}
num <- vapply(tab, is.numeric, TRUE)
tab[num] <- lapply(tab[num], function(x) round(x, 4))
print(tab)

bad <- Filter(function(r) identical(r$outcome, "refused"), R)
for (r in bad) {
  cat("\nREFUSED seed ", r$seed, ": ",
      substr(gsub("[[:space:]]+", " ", r$why), 1L, 300L), "\n", sep = "")
}
wn <- unlist(lapply(R, function(r) r$warns))
if (length(wn)) {
  cat("\nwarnings raised, ", length(wn), " over ", n, " runs:\n", sep = "")
  print(table(substr(gsub("[[:space:]]+", " ", wn), 1L, 70L)))
}
