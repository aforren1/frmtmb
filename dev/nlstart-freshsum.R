# Lane nlstart fix round: the review's freshsum with the reference (best of
# 16 to 17 started fits) taken from the review's own after-arm run, a[3],
# since the fix-round run of dev/review-nlstart-fresh.R is made without it.
# Usage: Rscript dev/nlstart-freshsum.R <new after dirs> <before dirs> <old after dirs>
# Review of lane nlstart: summary of review-nlstart-fresh*.R logs.
# Usage: Rscript dev/review-nlstart-freshsum.R <after dir> <before dir>
a <- commandArgs(trailingOnly = TRUE)
rd <- function(dirs) {
  out <- NULL
  for (k in seq_along(dirs)) {
    dd <- rd1(dirs[k]); dd$set <- k; out <- rbind(out, dd)
  }
  out
}
rd1 <- function(dir) {
  fs <- list.files(dir, "[.]log$", full.names = TRUE)
  for (f in fs) cat("read", f, ":", readLines(f, 1), "\n")
  ln <- unlist(lapply(fs, function(f) grep("^ROW ", readLines(f), value = TRUE)))
  p <- strsplit(ln, " ")
  data.frame(fam = sapply(p, `[`, 2), body = sapply(p, `[`, 3),
             seed = sapply(p, `[`, 4),
             ll = suppressWarnings(as.numeric(sapply(p, `[`, 5))),
             conv = sapply(p, `[`, 6), t = as.numeric(sapply(p, `[`, 7)),
             nwarn = as.integer(sapply(p, `[`, 8)), wclass = sapply(p, `[`, 9),
             multi = as.integer(sapply(p, `[`, 10)),
             ll_off = suppressWarnings(as.numeric(sapply(p, `[`, 11))),
             t_off = suppressWarnings(as.numeric(sapply(p, `[`, 13))),
             ref = suppressWarnings(as.numeric(sapply(p, `[`, 14))),
             est = sapply(p, `[`, 16), stringsAsFactors = FALSE)
}
A <- rd(strsplit(a[1], ",")[[1]]); B <- rd(strsplit(a[2], ",")[[1]])
R <- rd(strsplit(a[3], ",")[[1]])
A$ref <- R$ref[match(paste(A$set, A$fam, A$body, A$seed), paste(R$set, R$fam, R$body, R$seed))]
best_old <- pmax(R$ll, R$ll_off, na.rm = TRUE)[match(paste(A$set, A$fam, A$body, A$seed), paste(R$set, R$fam, R$body, R$seed))]
A$ref <- pmax(A$ref, best_old, na.rm = TRUE)
key <- function(d) paste(d$set, d$fam, d$body, d$seed)
B <- B[match(key(A), key(B)), ]
stopifnot(identical(key(A), key(B)))
best <- pmax(A$ref, A$ll, A$ll_off, B$ll, na.rm = TRUE)
ok <- !is.na(A$ll)
cat(sprintf("rows %d, default error (after) %d, (before) %d\n", nrow(A),
            sum(is.na(A$ll)), sum(is.na(B$ll))))
sh_a <- ok & best - A$ll > 1e-4
sh_b <- !is.na(B$ll) & best - B$ll > 1e-4
cat(sprintf("short > 1e-4: after %d (silent %d), before %d (silent %d)\n",
            sum(sh_a), sum(sh_a & A$nwarn == 0), sum(sh_b),
            sum(sh_b & B$nwarn == 0)))
cat(sprintf("short > 1: after %d (silent %d), before %d (silent %d)\n",
            sum(ok & best - A$ll > 1), sum(ok & best - A$ll > 1 & A$nwarn == 0),
            sum(!is.na(B$ll) & best - B$ll > 1),
            sum(!is.na(B$ll) & best - B$ll > 1 & B$nwarn == 0)))
dlt <- A$ll - B$ll
cat(sprintf("after - before: higher > 1e-4 %d, lower < -1e-4 %d, within %d\n",
            sum(dlt > 1e-4, na.rm = TRUE), sum(dlt < -1e-4, na.rm = TRUE),
            sum(abs(dlt) <= 1e-4, na.rm = TRUE)))
offd <- A$ll_off - B$ll
cat(sprintf("probe-off - before: |diff| > 1e-8 %d, max |diff| %.3g\n",
            sum(abs(offd) > 1e-8, na.rm = TRUE), max(abs(offd), na.rm = TRUE)))
cat(sprintf("multistart message on %d rows; of them gain over probe-off < 1e-2: %d\n",
            sum(A$multi == 1, na.rm = TRUE),
            sum(A$multi == 1 & A$ll - A$ll_off < 1e-2, na.rm = TRUE)))
cat("\nshort rows after (best = max of 16-17 starts and every arm):\n")
s <- A[sh_a, c("fam", "body", "seed", "ll", "conv", "nwarn", "multi")]
s$before <- B$ll[sh_a]; s$best <- best[sh_a]; s$short <- best[sh_a] - A$ll[sh_a]
s$est <- A$est[sh_a]
print(s[order(-s$short), ], row.names = FALSE)
cat("\nrows lower than before by > 1e-4:\n")
w <- which(dlt < -1e-4)
print(data.frame(A[w, c("fam", "body", "seed", "ll", "conv", "nwarn")],
                 before = B$ll[w], diff = dlt[w], est = A$est[w]),
      row.names = FALSE)
cat("\nshort by body (after / before):\n")
print(table(A$body[sh_a]))
print(table(B$body[sh_b]))
tr <- A$t / A$t_off
cat("\ntime on / off (same process, on first) by body, median:\n")
print(round(tapply(tr, A$body, median, na.rm = TRUE), 2))
cat(sprintf("overall sum on %.1f s, off %.1f s, ratio %.2f\n",
            sum(A$t, na.rm = TRUE), sum(A$t_off, na.rm = TRUE),
            sum(A$t, na.rm = TRUE) / sum(A$t_off, na.rm = TRUE)))
gm <- A$ll - A$ll_off
cat("\nmultistart rows, gain of the returned fit over probe-off (quantiles):\n")
print(signif(quantile(gm[A$multi == 1], c(0, 0.05, 0.25, 0.5, 1), na.rm = TRUE), 3))
sm <- which(A$multi == 1 & gm < 1e-2)
print(data.frame(A[sm, c("fam", "body", "seed", "ll", "ll_off")], gain = gm[sm]),
      row.names = FALSE)
cat(sprintf("rows where the probe ran and found nothing (multi 0, curved body): %d\n",
            sum(A$multi == 0 & ok & !A$body %in% character(0))))
