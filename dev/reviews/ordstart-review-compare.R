# Review of lane ordstart: compare a base and a lane run of
# ordstart-review-sweep.R case by case.
#   Rscript ordstart-review-compare.R base_prefix lane_prefix
args <- commandArgs(trailingOnly = TRUE)
B <- utils::read.delim(paste0(args[1], ".tsv"), stringsAsFactors = FALSE)
A <- utils::read.delim(paste0(args[2], ".tsv"), stringsAsFactors = FALSE)
pb <- readRDS(paste0(args[1], "-pars.rds"))
pa <- readRDS(paste0(args[2], "-pars.rds"))
M <- merge(B, A, by = "case", suffixes = c("_b", "_l"))
M$dll <- M$ll_l - M$ll_b
M$dpar <- vapply(M$case, function(cs) {
  a <- pa[[cs]]
  b <- pb[[cs]]
  if (is.null(a) || is.null(b) || length(a) != length(b)) return(NA_real_)
  max(abs(a - b))
}, 0)
norm <- function(s) {
  w <- strsplit(s, " || ", fixed = TRUE)[[1]]
  w <- w[!grepl("jittered starts reached a higher", w)]
  w <- gsub("-?[0-9][0-9.e+-]*", "#", w)
  sort(unique(w))
}
M$new_warn <- vapply(seq_len(nrow(M)), function(i) {
  paste(setdiff(norm(M$warn_l[i]), norm(M$warn_b[i])), collapse = " | ")
}, "")
M$lost_warn <- vapply(seq_len(nrow(M)), function(i) {
  paste(setdiff(norm(M$warn_b[i]), norm(M$warn_l[i])), collapse = " | ")
}, "")
M$time_ratio <- M$elapsed_l / M$elapsed_b
out <- M[, c("case", "ll_b", "ll_l", "dll", "dpar", "code_b", "code_l",
             "degen_b", "degen_l", "stepw_l", "time_ratio")]
op <- options(width = 220)
print(out, digits = 6, row.names = FALSE)
ok_b <- M$degen_b %in% FALSE
cat("\nfits", nrow(M), "; errors base", sum(is.na(M$ll_b)), "lane",
    sum(is.na(M$ll_l)), "\n")
cat("logLik moved > 1e-6:", sum(abs(M$dll) > 1e-6, na.rm = TRUE),
    "; up:", sum(M$dll > 1e-6, na.rm = TRUE), "; down:",
    sum(M$dll < -1e-6, na.rm = TRUE), "\n")
cat("base passed degenerate check:", sum(ok_b), "\n")
cat("  of those, logLik down > 1e-6:", sum(ok_b & M$dll < -1e-6, na.rm = TRUE),
    "; gained degenerate warning:", sum(ok_b & M$degen_l %in% TRUE),
    "; gained step-start warning:", sum(ok_b & M$stepw_l %in% TRUE),
    "; estimates moved (max|dpar| > 1e-4) at |dll| <= 1e-6:",
    sum(ok_b & abs(M$dll) <= 1e-6 & M$dpar > 1e-4, na.rm = TRUE), "\n")
cat("base degenerate:", sum(M$degen_b %in% TRUE), "; lane degenerate:",
    sum(M$degen_l %in% TRUE), "; base degenerate and lane lower:",
    sum(M$degen_b %in% TRUE & M$dll < -1e-6, na.rm = TRUE),
    "; base degenerate and lane still degenerate but lower:",
    sum(M$degen_b %in% TRUE & M$degen_l %in% TRUE & M$dll < -1e-6,
        na.rm = TRUE), "\n")
cat("convergence code 0 -> nonzero:", sum(M$code_b == 0 & M$code_l != 0,
                                          na.rm = TRUE),
    "; nonzero -> 0:", sum(M$code_b != 0 & M$code_l == 0, na.rm = TRUE), "\n")
cat("time ratio median", stats::median(M$time_ratio, na.rm = TRUE),
    "range", paste(signif(range(M$time_ratio, na.rm = TRUE), 3),
                   collapse = " .. "), "\n")
cat("\nNEW warnings (lane, not base; step-start excluded, numbers masked):\n")
for (i in which(nzchar(M$new_warn))) cat(" ", M$case[i], ":", M$new_warn[i], "\n")
cat("\nLOST warnings (base, not lane):\n")
for (i in which(nzchar(M$lost_warn))) cat(" ", M$case[i], ":", M$lost_warn[i], "\n")
