# Draw the survey sample from the core pair table.
#
# Design: stratified by compat status; within a stratum, n_h draws WITH
# replacement, each pair drawn with probability proportional to its
# reach weight (dev/equiv-weights.R). Every draw is then an independent
# Bernoulli trial whose mean is the reach-weighted equivalence rate of
# the stratum, so the plain share over draws estimates it without bias
# and a Wilson interval on n_h draws is valid. A pair drawn twice is
# measured once and counted twice: the cost is the unique pairs, the
# precision is the draws.
#
# The sample is frozen in sample.tsv so that a later build reruns the
# SAME pairs; redraw only when the pair table itself changes.
#
# Usage: Rscript dev/equiv-sample.R [outdir] [seed]
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "dev/equiv-data"
seed <- if (length(args) >= 2) as.integer(args[2]) else 20261028L

cp <- read.delim("dev/equiv-frame-out/compat.tsv")
cp <- cp[cp$in_core_only, ]
# The sample was drawn from the weights frozen in weights-drawn.tsv.
# weights.tsv is the current reach (the bar_crossing pattern was
# corrected after the draw, review 2026-10-28); redraw from it only at
# the next frame change, by passing it as the third argument.
wfile <- if (length(args) >= 3) args[3] else
  file.path(out, "weights-drawn.tsv")
w <- read.delim(wfile)
cp$w <- w$weight[match(cp$feature_a, w$feature)] *
  w$weight[match(cp$feature_b, w$feature)]
stopifnot(!anyNA(cp$w))

n_draws <- c(untested = 70L, conditional = 70L, works = 40L,
             refused = 20L)
set.seed(seed)
rows <- list()
for (s in names(n_draws)) {
  idx <- which(cp$status == s)
  d <- sample(idx, n_draws[[s]], replace = TRUE, prob = cp$w[idx])
  tab <- table(d)
  ii <- as.integer(names(tab))
  rows[[s]] <- data.frame(cp[ii, c("feature_a", "kind_a", "feature_b",
                                   "kind_b", "status")],
                          weight = cp$w[ii], draws = as.integer(tab))
}
smp <- do.call(rbind, rows)
smp <- smp[order(match(smp$status, names(n_draws)), -smp$draws,
                 smp$feature_a, smp$feature_b), ]
smp$id <- sprintf("p%03d", seq_len(nrow(smp)))
smp <- smp[, c("id", setdiff(names(smp), "id"))]
rownames(smp) <- NULL

# stratum frame: count of pairs and share of total reach weight
fr <- aggregate(w ~ status, cp, function(x) c(n = length(x), w = sum(x)))
fr <- data.frame(status = fr$status, N = fr$w[, "n"], W = fr$w[, "w"])
fr$W_share <- fr$W / sum(fr$W)
fr$N_share <- fr$N / sum(fr$N)
fr$draws <- n_draws[fr$status]
fr$unique <- as.integer(table(smp$status)[fr$status])
print(fr, row.names = FALSE)
cat("unique pairs:", nrow(smp), " draws:", sum(smp$draws), " seed:", seed,
    "\n")
utils::write.table(smp, file.path(out, "sample.tsv"), sep = "\t",
                   row.names = FALSE, quote = TRUE, qmethod = "double")
utils::write.table(fr, file.path(out, "strata.tsv"), sep = "\t",
                   row.names = FALSE, quote = TRUE, qmethod = "double")
feats <- sort(table(c(smp$feature_a, smp$feature_b)), decreasing = TRUE)
cat(length(feats), "features in the sample\n")
print(feats)
