# Lane ordstarts (backlog 6136): the tables of the main measurement, from
# the logs of dev/ordstarts-run.R (dev/ordstarts-out/main-*.tsv). Writes
# the moved seeds to dev/ordstarts-out/moves.tsv for
# dev/ordstarts-refits.R.
#   Rscript dev/ordstarts-sum.R > dev/ordstarts-out/sum.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
fs <- Sys.glob("dev/ordstarts-out/main-*.tsv")
X <- do.call(rbind, lapply(fs, function(f) {
  utils::read.delim(f, stringsAsFactors = FALSE, quote = "",
                    comment.char = "")
}))
X$design <- sprintf("%-10s %-5s n%-4d G%d", X$family, X$re, X$n, X$G)
A <- X[X$arm == "default", ]
B <- X[X$arm == "starts10", ]
key <- function(D) paste(D$family, D$re, D$n, D$seed)
# a design still running has a default row without its pair
A <- A[key(A) %in% key(B), ]
B <- B[match(key(A), key(B)), ]
stopifnot(identical(key(A), key(B)))
tse <- function(D) suppressWarnings(pmax(D$tse1, D$tse2, na.rm = TRUE))
P <- data.frame(design = A$design, family = A$family, re = A$re, n = A$n,
  G = A$G, seed = A$seed,
  llA = A$ll, llB = B$ll, d = B$ll - A$ll, codeA = A$code, codeB = B$code,
  clsA = A$cls, clsB = B$cls, ran = !is.na(B$ms_starts),
  reached = B$ms_reached, above = B$ms_degenerate_above,
  warnA = A$nwarn, degA = A$degen_warn, degB = B$degen_warn,
  hintA = A$re_hint, convA = A$conv_warn, stepB = B$step_msg,
  tseA = tse(A), tseB = tse(B), maxseA = A$maxse, maxseB = B$maxse,
  eigA = A$mineig, eigB = B$mineig, defA = A$se_deferred,
  defB = B$se_deferred, elA = A$elapsed, elB = B$elapsed,
  evA = A$evals, evB = B$evals, errA = A$err, errB = B$err,
  stringsAsFactors = FALSE)
P$moved <- abs(P$d) > 1e-6
P$up <- P$d > 1e-3
P$down <- P$d < -1e-3
cat("fits per arm", nrow(P), "; errors default", sum(P$errA != "" &
    !is.na(P$errA)), "starts10", sum(P$errB != "" & !is.na(P$errB)), "\n\n")
cat("== Per design (40 seeds each) ==\n")
cat("code1A: default code not 0; ran: the starts arm ran starts;",
    "reached: converged starts (sum);\nmoved: |dll| > 1e-6; up/down: by",
    "more than 1e-3; silent: up and the default fit gave no warning;\n",
    "degA/degB: degenerate warning; stepB: step-start message\n\n")
by_design <- split(P, P$design)
for (nm in names(by_design)) {
  D <- by_design[[nm]]
  sil <- D$up & D$warnA == 0
  cat(sprintf(paste("%s | err %d | code1A %2d | clsA proper %2d step %2d",
                    "flat %2d | ran %2d reached %3d of %3d | moved %2d up",
                    "%2d down %2d | silent %2d | degA %2d degB %2d hintA",
                    "%2d stepB %2d\n"),
              nm, sum(D$errA != "" & !is.na(D$errA)),
              sum(D$codeA != 0, na.rm = TRUE),
              sum(D$clsA == "proper", na.rm = TRUE),
              sum(D$clsA == "step", na.rm = TRUE),
              sum(D$clsA == "flat", na.rm = TRUE), sum(D$ran),
              sum(D$reached, na.rm = TRUE), 10L * sum(D$ran),
              sum(D$moved, na.rm = TRUE), sum(D$up, na.rm = TRUE),
              sum(D$down, na.rm = TRUE), sum(sil, na.rm = TRUE),
              sum(D$degA), sum(D$degB), sum(D$hintA), sum(D$stepB)))
}
cat("\n== All designs ==\n")
sil <- P$up & P$warnA == 0
cat(sprintf(paste("fits %d | code1A %d | ran %d | moved %d up %d down %d |",
                  "silent %d | degA %d degB %d | stepB %d\n"),
            nrow(P), sum(P$codeA != 0, na.rm = TRUE), sum(P$ran),
            sum(P$moved, na.rm = TRUE), sum(P$up, na.rm = TRUE),
            sum(P$down, na.rm = TRUE), sum(sil, na.rm = TRUE), sum(P$degA),
            sum(P$degB), sum(P$stepB)))
cat("class default -> starts10:\n")
print(table(default = P$clsA, starts10 = P$clsB))
cat("code default -> starts10:\n")
print(table(default = P$codeA, starts10 = P$codeB, useNA = "ifany"))
cat("\n== Moves (|dll| > 1e-6) ==\n")
M <- P[P$moved %in% TRUE, ]
M <- M[order(M$design, M$seed), ]
fmt <- function(v, dg = 3) ifelse(is.na(v), "NA", format(signif(v, dg)))
for (j in seq_len(nrow(M))) {
  m <- M[j, ]
  cat(sprintf(paste("%s seed %2d | %.4f -> %.4f d %+.4f | code %d->%d |",
                    "%s->%s | thrSE %s->%s | maxSE %s->%s | mineig %s->%s",
                    "| warnA %d degA %s degB %s stepB %s\n"),
              m$design, m$seed, m$llA, m$llB, m$d, m$codeA, m$codeB, m$clsA,
              m$clsB, fmt(m$tseA), fmt(m$tseB), fmt(m$maxseA),
              fmt(m$maxseB), fmt(m$eigA), fmt(m$eigB), m$warnA, m$degA,
              m$degB, m$stepB))
}
if (nrow(M)) {
  cat("\nup moves: d min", fmt(min(M$d[M$up])), "median",
      fmt(stats::median(M$d[M$up])), "max", fmt(max(M$d[M$up])), "\n")
  if (any(M$down)) {
    cat("down moves: d min", fmt(min(M$d[M$down])), "max",
        fmt(max(M$d[M$down])), "\n")
  }
}
cat("\n== Cost (single runs, two processes at once) ==\n")
cat("per design: elapsed sum default, sum starts10, ratio of sums;",
    "per-fit ratio median and range\nover the fits where starts ran;",
    "evals ratio median over the same\n")
for (nm in names(by_design)) {
  D <- by_design[[nm]]
  r <- D[D$ran, ]
  q <- r$elB / r$elA
  e <- r$evB / r$evA
  cat(sprintf(paste("%s | default %7.1fs starts10 %7.1fs ratio %5.1f |",
                    "ran %2d: per-fit ratio median %5.1f range %5.1f to",
                    "%5.1f | evals ratio median %5.1f\n"),
              nm, sum(D$elA), sum(D$elB), sum(D$elB) / sum(D$elA), nrow(r),
              if (nrow(r)) stats::median(q) else NA,
              if (nrow(r)) min(q) else NA, if (nrow(r)) max(q) else NA,
              if (nrow(r)) stats::median(e, na.rm = TRUE) else NA))
}
r <- P[P$ran, ]
cat(sprintf(paste("ALL | default %.1fs starts10 %.1fs ratio %.1f | ran %d:",
                  "per-fit ratio median %.1f range %.1f to %.1f\n"),
            sum(P$elA), sum(P$elB), sum(P$elB) / sum(P$elA), nrow(r),
            stats::median(r$elB / r$elA), min(r$elB / r$elA),
            max(r$elB / r$elA)))
cat("\n== Standard-error check timing ==\n")
cat("deferred at the fit (no covariance built): default",
    sum(P$defA, na.rm = TRUE), "starts10", sum(P$defB, na.rm = TRUE), "of",
    sum(!is.na(P$defA)), "fitted\n")
cat("default fits deferred whose threshold SE (forced) is above 100:",
    sum(P$defA & P$tseA > 100, na.rm = TRUE), "\n")
cat("default fits whose threshold SE (forced covariance) is above 100:",
    sum(P$tseA > 100, na.rm = TRUE), "; of those with a degenerate",
    "warning:", sum(P$tseA > 100 & P$degA, na.rm = TRUE), "\n")
cat("starts10 fits whose threshold SE is above 100:",
    sum(P$tseB > 100, na.rm = TRUE), "; with a degenerate warning:",
    sum(P$tseB > 100 & P$degB, na.rm = TRUE), "\n")
utils::write.table(M[, c("family", "re", "n", "G", "seed", "d", "clsA",
                         "clsB", "codeA", "codeB")],
                   "dev/ordstarts-out/moves.tsv", sep = "\t", quote = FALSE,
                   row.names = FALSE)
cat("\nwrote", nrow(M), "moves to dev/ordstarts-out/moves.tsv\n")
