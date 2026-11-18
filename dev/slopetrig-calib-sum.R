# Lane slopetrig: summarize the calibration grid (dev/slopetrig-calib.R).
# Short is against the best log-likelihood of the four arms of a case,
# by more than 1e-3. Prints one line per cell, then, for candidate ratio
# thresholds R, what a rule "spread < 0.05 OR ratio > R" engages, what
# it misses and what it costs where it engages and is not needed.
# Usage: Rscript slopetrig-calib-sum.R <tsv> [<tsv> ...]
fs <- commandArgs(TRUE)
d <- do.call(rbind, lapply(fs, utils::read.delim))
best <- pmax(d$ll_ref, d$ll_F, d$ll_T, d$ll_D, na.rm = TRUE)
k <- 1e-3
d$sh_F <- best - d$ll_F
d$sh_T <- best - d$ll_T
d$sh_D <- best - d$ll_D
d$sh_ref <- best - d$ll_ref
d$short_F <- d$sh_F > k
d$tr <- d$t_T / d$t_F
cat("cases", nrow(d), "; NA ll: F", sum(is.na(d$ll_F)), "T", sum(is.na(d$ll_T)),
    "D", sum(is.na(d$ll_D)), "ref", sum(is.na(d$ll_ref)), "\n")
cat("short > 1e-3: F", sum(d$short_F, na.rm = TRUE),
    "T", sum(d$sh_T > k, na.rm = TRUE), "D(base)", sum(d$sh_D > k, na.rm = TRUE),
    "ref", sum(d$sh_ref > k, na.rm = TRUE), "; max short T",
    signif(max(d$sh_T, na.rm = TRUE), 3), "\n")
key <- paste(d$design, d$c0, d$r, d$fam, d$form)
cell <- do.call(rbind, lapply(split(seq_len(nrow(d)), factor(key, unique(key))),
  function(i) {
    x <- d[i, ]
    data.frame(design = x$design[1], c0 = x$c0[1], r = signif(x$r[1], 2),
               fam = x$fam[1], form = x$form[1], n = nrow(x),
               ratio = signif(stats::median(x$ratio), 3),
               spread_min = signif(min(x$spread), 2),
               shF = sum(x$short_F), silF = sum(x$short_F & x$w_F == 0),
               maxF = signif(max(x$sh_F), 3), shT = sum(x$sh_T > k),
               wT = sum(x$w_T), wF = sum(x$w_F),
               engB = sum(x$as_D), shB = sum(x$sh_D > k),
               tT_tF = signif(stats::median(x$tr), 3))
  }))
options(width = 220)
print(cell, row.names = FALSE)
sp <- d$spread < 0.05
cat("\nby spread rule alone (spread < 0.05): engaged", sum(sp),
    "; short F not engaged", sum(d$short_F & !sp), "\n")
nsp <- d[!sp, ]
cat("among spread >= 0.05: short F", sum(nsp$short_F), "of", nrow(nsp),
    "; min ratio of a short F", signif(min(nsp$ratio[nsp$short_F]), 4),
    "; max ratio of a not-short F", signif(max(nsp$ratio[!nsp$short_F]), 4),
    "\n")
cat("short F cases with spread >= 0.05, sorted by ratio:\n")
s <- nsp[nsp$short_F, c("design", "c0", "r", "fam", "form", "seed", "ratio",
                        "spread", "sh_F", "w_F")]
print(s[order(s$ratio), ][seq_len(min(15, nrow(s))), ], row.names = FALSE)
cat("not-short F cases with spread >= 0.05 and ratio > 10, by ratio:\n")
s2 <- nsp[!nsp$short_F & nsp$ratio > 10, c("design", "c0", "r", "fam", "form",
                                           "seed", "ratio", "spread", "sh_F")]
print(utils::tail(s2[order(s2$ratio), ], 15), row.names = FALSE)
for (R in c(5, 10, 15, 20, 25, 30, 40, 50, 100)) {
  eng <- sp | d$ratio > R
  byr <- !sp & d$ratio > R
  cat(sprintf(paste0("R %4g: engaged %d (by ratio %d) ; short F missed %d",
                     " (max %.3g) ; engaged not needed %d (by ratio %d)",
                     " ; median T/F time where engaged by ratio and not",
                     " needed %.3g ; short T where engaged %d\n"),
              R, sum(eng), sum(byr), sum(d$short_F & !eng),
              if (any(d$short_F & !eng)) max(d$sh_F[d$short_F & !eng]) else 0,
              sum(eng & !d$short_F), sum(byr & !d$short_F),
              if (any(byr & !d$short_F)) {
                stats::median(d$tr[byr & !d$short_F])
              } else NA,
              sum(eng & d$sh_T > k)))
}
