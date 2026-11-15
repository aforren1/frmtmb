# Lane nlstop2: the acceptance sets of lane nlstop on both arms.
#   fresh: the review's 425 fresh rows (dev/nlstop2-fresh.R, -fresh2.R);
#          the reference is the best of 16 or 17 started fits of lane
#          nlstop's before arm (frmtmb-wt-nlstop/dev/nlstop-out/fresh*-before,
#          same data) and of both arms' defaults.
#   scan:  the 450 rows (dev/nlstop2-scan.R); reference the best of 8
#          random starts of lane nlstop's before arm.
#   limit: the review's limit battery (dev/nlstop2-limit.R), warning counts.
# Usage: Rscript dev/nlstop2-accsum.R <after tag> <before tag>
#   (directories dev/nlstop2-out/<set>-<tag>)
a <- commandArgs(TRUE)
ta <- a[1]; tb <- a[2]
out <- "dev/nlstop2-out"
old <- "C:/Users/adf44/source/r/frmtmb-wt-nlstop/dev/nlstop-out"
rows <- function(dirs, pat) {
  r <- NULL
  for (k in seq_along(dirs)) {
    fs <- list.files(dirs[k], pat, full.names = TRUE)
    for (f in fs) {
      l <- readLines(f)
      cat("read", f, ":", l[1], "|", tail(l, 1), "\n")
      ln <- grep("^ROW ", l, value = TRUE)
      if (!length(ln)) next
      p <- strsplit(ln, " ")
      r <- rbind(r, data.frame(set = k, f = I(p), stringsAsFactors = FALSE))
    }
  }
  r
}
col <- function(r, i, num = TRUE) {
  v <- vapply(r$f, function(x) if (length(x) >= i) x[i] else NA_character_,
              "")
  if (num) suppressWarnings(as.numeric(v)) else v
}
cat("== fresh rows\n")
fa <- rows(file.path(out, c(paste0("fresh-", ta), paste0("fresh2-", ta))),
           "[.]log$")
fb <- rows(file.path(out, c(paste0("fresh-", tb), paste0("fresh2-", tb))),
           "[.]log$")
fr <- rows(file.path(old, c("fresh-before", "fresh2-before")), "[.]log$")
key <- function(r) paste(r$set, col(r, 2, FALSE), col(r, 3, FALSE),
                         col(r, 4, FALSE))
fb <- fb[match(key(fa), key(fb)), ]
fr <- fr[match(key(fa), key(fr)), ]
stopifnot(identical(key(fa), key(fb)), identical(key(fa), key(fr)))
A <- col(fa, 5); B <- col(fb, 5); R <- col(fr, 14)
best <- pmax(R, A, B, na.rm = TRUE)
wa <- col(fa, 8); wb <- col(fb, 8)
cat(sprintf("rows %d, default error after %d, before %d\n", nrow(fa),
            sum(is.na(A)), sum(is.na(B))))
for (arm in c("after", "before")) {
  L <- if (arm == "after") A else B
  W <- if (arm == "after") wa else wb
  s <- !is.na(L) & best - L > 1e-4
  s1 <- !is.na(L) & best - L > 1
  cat(sprintf("%s: short > 1e-4 %d (silent %d), short > 1 %d (silent %d)\n",
              arm, sum(s), sum(s & W == 0), sum(s1), sum(s1 & W == 0)))
}
d <- A - B
cat(sprintf("after - before: higher > 1e-4 %d, lower < -1e-4 %d, within %d\n",
            sum(d > 1e-4, na.rm = TRUE), sum(d < -1e-4, na.rm = TRUE),
            sum(abs(d) <= 1e-4, na.rm = TRUE)))
wca <- col(fa, 9, FALSE); wcb <- col(fb, 9, FALSE)
cat(sprintf("frmtmb_nl_limit after %d before %d; frmtmb_nl_pole after %d\n",
            sum(grepl("nl_limit", wca)), sum(grepl("nl_limit", wcb)),
            sum(grepl("nl_pole", wca))))
cat(sprintf("time of default fits: after %.1f s, before %.1f s\n",
            sum(col(fa, 7), na.rm = TRUE), sum(col(fb, 7), na.rm = TRUE)))
show <- function(i) {
  data.frame(fam = col(fa, 2, FALSE)[i], body = col(fa, 3, FALSE)[i],
             seed = col(fa, 4, FALSE)[i], after = A[i], before = B[i],
             best = best[i], short_after = best[i] - A[i],
             warn_after = wca[i], warn_before = wcb[i])
}
cat("rows whose after differs from before by > 1e-4:\n")
print(show(which(abs(d) > 1e-4)), row.names = FALSE)
cat("rows short after by > 1e-4:\n")
print(show(which(!is.na(A) & best - A > 1e-4)), row.names = FALSE)

cat("\n== scan rows\n")
sa <- rows(file.path(out, paste0("scan-", ta)), "[.]txt$")
sb <- rows(file.path(out, paste0("scan-", tb)), "[.]txt$")
sr <- rows(file.path(old, "scan-before"), "[.]txt$")
k2 <- function(r) paste(col(r, 2, FALSE), col(r, 3, FALSE), col(r, 4, FALSE),
                        col(r, 5, FALSE))
sb <- sb[match(k2(sa), k2(sb)), ]
sr <- sr[match(k2(sa), k2(sr)), ]
stopifnot(identical(k2(sa), k2(sb)), identical(k2(sa), k2(sr)))
A <- col(sa, 6); B <- col(sb, 6); R <- col(sr, 7)
best <- pmax(R, A, B, na.rm = TRUE)
cat(sprintf("rows %d, default NA after %d before %d\n", nrow(sa),
            sum(is.na(A)), sum(is.na(B))))
for (arm in c("after", "before")) {
  L <- if (arm == "after") A else B
  W <- col(if (arm == "after") sa else sb, 10)
  s <- !is.na(L) & best - L > 1e-4
  cat(sprintf("%s: short > 1e-4 %d (silent %d), warned rows %d\n", arm,
              sum(s), sum(s & W == 0), sum(W > 0, na.rm = TRUE)))
}
d <- A - B
cat(sprintf("after - before: higher > 1e-4 %d, lower < -1e-4 %d, within %d\n",
            sum(d > 1e-4, na.rm = TRUE), sum(d < -1e-4, na.rm = TRUE),
            sum(abs(d) <= 1e-4, na.rm = TRUE)))
wca <- col(sa, 11, FALSE); wcb <- col(sb, 11, FALSE)
cat(sprintf("warning classes differ on %d rows\n", sum(wca != wcb)))
if (any(wca != wcb)) {
  i <- which(wca != wcb)
  print(data.frame(fam = col(sa, 2, FALSE)[i], seed = col(sa, 3, FALSE)[i],
                   slope = col(sa, 4, FALSE)[i], body = col(sa, 5, FALSE)[i],
                   after = A[i], before = B[i], wa = wca[i], wb = wcb[i]),
        row.names = FALSE)
}
cat(sprintf("time of default fits: after %.1f s, before %.1f s\n",
            sum(col(sa, 9), na.rm = TRUE), sum(col(sb, 9), na.rm = TRUE)))

cat("\n== limit battery\n")
for (f in c("A", "B", "Cprior", "Cstart")) {
  for (tag in c(ta, tb)) {
    l <- readLines(file.path(out, paste0("limit-", tag), paste0(f, ".txt")))
    r <- strsplit(grep("^ROW ", l, value = TRUE), " ")
    code <- vapply(r, `[`, "", 7); wc <- vapply(r, `[`, "", 8)
    cat(sprintf("%s %s: rows %d err %d warned %d nl_limit %d nl_pole %d\n",
                f, tag, length(r), sum(code == "ERR"),
                sum(code != "ERR" & wc != "-"), sum(grepl("nl_limit", wc)),
                sum(grepl("nl_pole", wc))))
  }
}
