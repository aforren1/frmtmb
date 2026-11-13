# Lane ordstart, fix round: summarise one set of the review's sweep
# (dev/reviews/ordstart-review-sweep.R), base against lane, case by case:
# logLik moves, convergence-code changes, the degenerate warning gained
# or lost (the "degen" column of the sweep matches "ended at a
# degenerate boundary" only when the warning text has it, so it is read
# from the warning text here), the step-start warning, and the largest
# |outer parameter| at each moved point.
#   Rscript dev/ordstart-fix-sweep-sum.R base_prefix lane_prefix [...]
args <- commandArgs(trailingOnly = TRUE)
pre <- matrix(args, 2)
rows <- NULL
for (j in seq_len(ncol(pre))) {
  B <- utils::read.delim(paste0(pre[1, j], ".tsv"), stringsAsFactors = FALSE)
  A <- utils::read.delim(paste0(pre[2, j], ".tsv"), stringsAsFactors = FALSE)
  pa <- readRDS(paste0(pre[2, j], "-pars.rds"))
  M <- merge(B, A, by = "case", suffixes = c("_b", "_l"))
  M$maxpar_l <- vapply(M$case, function(cs) {
    if (is.null(pa[[cs]])) NA_real_ else max(abs(pa[[cs]]))
  }, 0)
  rows <- rbind(rows, M)
}
dg <- function(w) grepl("ended at a degenerate boundary", w, fixed = TRUE)
rows$degen_b <- dg(rows$warn_b)
rows$degen_l <- dg(rows$warn_l)
rows$dll <- rows$ll_l - rows$ll_b
mv <- abs(rows$dll) > 1e-6
op <- options(width = 200)
print(rows[mv, c("case", "ll_b", "ll_l", "dll", "code_b", "code_l",
                 "degen_b", "degen_l", "stepw_l", "maxpar_l")],
      digits = 7, row.names = FALSE)
passed <- !rows$degen_b
cat("\nfits", nrow(rows), "; passed by 0.81.0", sum(passed),
    "; logLik moved", sum(mv), "(up", sum(rows$dll > 1e-6), ", down",
    sum(rows$dll < -1e-6), ")\n")
cat("passed fits: logLik down", sum(passed & rows$dll < -1e-6),
    "; gained the degenerate warning", sum(passed & rows$degen_l), "\n")
cat("degenerate warning: base", sum(rows$degen_b), "lane", sum(rows$degen_l),
    "; lost on", sum(rows$degen_b & !rows$degen_l), "; step-start warning",
    sum(rows$stepw_l), "(on passed fits", sum(passed & rows$stepw_l), ")\n")
cat("convergence code 0 -> not 0:", sum(rows$code_b == 0 & rows$code_l != 0),
    "; not 0 -> 0:", sum(rows$code_b != 0 & rows$code_l == 0), "\n")
cat("moved down to a point with max|par| >= 740:",
    sum(rows$dll < -1e-6 & rows$maxpar_l >= 740, na.rm = TRUE), "\n")
r <- rows$elapsed_l / rows$elapsed_b
cat("time ratio lane/base: median", round(stats::median(r), 2), "range",
    paste(round(range(r), 2), collapse = " to "), "\n")
