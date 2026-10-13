# Paired comparison: merged build (rellib-r8) against the lane's recorded
# cox() study rds (wt-moparity-lib), seed by seed.
r8 <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel070-log/studies/moparity/dev/moparity-log/cox-study-lane.rds")
ln <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-moparity/dev/moparity-log/cox-study-lane.rds")
stopifnot(identical(r8$seed, ln$seed))
d <- r8$ll - ln$ll
cat(sprintf("seeds %d | r8 - lane logLik: identical %d, |d| > 0 %d, > 1e-10 %d, > 1e-8 %d, max %.3g, min %.3g\n",
            nrow(r8), sum(r8$ll == ln$ll), sum(d != 0), sum(abs(d) > 1e-10),
            sum(abs(d) > 1e-8), max(d), min(d)))
cat(sprintf("gap (polished max - fit): r8 max %.3g, lane max %.3g; r8 > 1e-6 %d, lane %d\n",
            max(r8$gap), max(ln$gap), sum(r8$gap > 1e-6), sum(ln$gap > 1e-6)))
cat("codes r8:", paste(names(table(r8$code)), table(r8$code), sep = "="),
    "| lane:", paste(names(table(ln$code)), table(ln$code), sep = "="), "\n")
cat("warned r8:", sum(nzchar(r8$warn)), "| lane:", sum(nzchar(ln$warn)), "\n")
cat("evals r8 total", sum(r8$evals), "median", median(r8$evals),
    "| lane total", sum(ln$evals), "median", median(ln$evals),
    "| seeds with different evals", sum(r8$evals != ln$evals), "\n")
de <- r8$evals - ln$evals
if (any(de != 0)) {
  cat("evals r8 - lane, per differing seed:\n")
  print(data.frame(seed = r8$seed[de != 0], r8 = r8$evals[de != 0],
                   lane = ln$evals[de != 0], d_ll = d[de != 0]))
}
big <- order(-abs(d))[1:5]
cat("largest |d logLik| seeds:\n")
print(data.frame(seed = r8$seed[big], ll_r8 = sprintf("%.12f", r8$ll[big]),
                 ll_lane = sprintf("%.12f", ln$ll[big]), d = d[big],
                 gap_r8 = r8$gap[big], gap_lane = ln$gap[big]))
cat("minw < 1e-6 at fit r8", sum(r8$minw < 1e-6), "lane", sum(ln$minw < 1e-6), "\n")
