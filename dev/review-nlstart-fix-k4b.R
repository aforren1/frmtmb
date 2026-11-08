# Review of lane nlstart, short check: rows with four curved coefficients.
# Usage: Rscript dev/review-nlstart-fix-k4b.R <fix dirs> <first-build dirs>
a <- commandArgs(trailingOnly = TRUE)
rd1 <- function(dir) {
  fs <- list.files(dir, "[.]log$", full.names = TRUE)
  for (f in fs) cat("read", f, ":", readLines(f, 1), "\n")
  ln <- unlist(lapply(fs, function(f) grep("^ROW ", readLines(f), value = TRUE)))
  p <- strsplit(ln, " ")
  data.frame(fam = sapply(p, `[`, 2), body = sapply(p, `[`, 3),
             seed = sapply(p, `[`, 4),
             ll = suppressWarnings(as.numeric(sapply(p, `[`, 5))),
             t = as.numeric(sapply(p, `[`, 7)),
             nwarn = as.integer(sapply(p, `[`, 8)),
             multi = as.integer(sapply(p, `[`, 10)),
             ll_off = suppressWarnings(as.numeric(sapply(p, `[`, 11))),
             t_off = suppressWarnings(as.numeric(sapply(p, `[`, 13))),
             ref = suppressWarnings(as.numeric(sapply(p, `[`, 14))))
}
rd <- function(dirs) {
  out <- NULL
  for (k in seq_along(dirs)) { d <- rd1(dirs[k]); d$set <- k; out <- rbind(out, d) }
  out
}
A <- rd(strsplit(a[1], ",")[[1]]); O <- rd(strsplit(a[2], ",")[[1]])
key <- function(d) paste(d$set, d$fam, d$body, d$seed)
O <- O[match(key(A), key(O)), ]
k4 <- A$body %in% c("sig4", "two_sat", "two_sat_rev", "five") & !is.na(A$ll)
best <- pmax(A$ref, O$ref, A$ll, O$ll, A$ll_off, na.rm = TRUE)
g_off <- A$ll - A$ll_off
g_old <- A$ll - O$ll
cat(sprintf("rows with four curved coefficients: %d\n", sum(k4)))
cat(sprintf("  fix build gains over probe off by > 0.01: %d; over the first build (4 patterns) by > 0.01: %d\n",
            sum(k4 & g_off > 0.01), sum(k4 & g_old > 0.01)))
cat(sprintf("  short of the best by > 1e-4: fix %d, first build %d, probe off %d\n",
            sum(k4 & best - A$ll > 1e-4), sum(k4 & best - O$ll > 1e-4),
            sum(k4 & best - A$ll_off > 1e-4)))
cat("  gain over first build > 0.01, by body:\n")
print(table(A$body[k4 & g_old > 0.01]))
cat("  gain over probe off > 0.01, by body:\n")
print(table(A$body[k4 & g_off > 0.01]))
r <- A$t / A$t_off
cat(sprintf("  time on/off median: fix %.2f, first build %.2f; sums fix %.1f/%.1f s\n",
            median(r[k4]), median((O$t / O$t_off)[k4]),
            sum(A$t[k4]), sum(A$t_off[k4])))
nk <- !k4 & !is.na(A$ll)
cat(sprintf("rows with 2 or 3 curved coefficients: %d, gain over off > 0.01: %d, time on/off median %.2f\n",
            sum(nk), sum(nk & g_off > 0.01), median(r[nk])))
# a cheaper trigger: run all 16 only when the first fit (0.79.1's default,
# which probe off equals) warned. Base logs as a third argument.
if (length(a) > 2) {
  Bs <- rd(strsplit(a[3], ",")[[1]])
  Bs <- Bs[match(key(A), key(Bs)), ]
  need <- k4 & g_old > 0.01
  cat(sprintf("rows that needed more than the first build's 4 patterns: %d; their first fit warned: %d, silent: %d (gain on silent: %s)\n",
              sum(need), sum(need & Bs$nwarn > 0), sum(need & Bs$nwarn == 0),
              paste(signif(sort(g_old[need & Bs$nwarn == 0]), 3), collapse = " ")))
  g2 <- !k4 & !is.na(A$ll) & g_off > 0.01
  cat(sprintf("2-3 coefficient rows that gained: %d; first fit warned: %d, silent: %d\n",
              sum(g2), sum(g2 & Bs$nwarn > 0), sum(g2 & Bs$nwarn == 0)))
}
