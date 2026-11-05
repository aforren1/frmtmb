# Scores the rules on the logs of dev/review-refuse-f2-rule.R. Hit: warned at a
# largest shift >= 1 exact SE. False alarm: warned below 0.5. Between: either.
fs <- list.files("dev/review-refuse-out/f2", pattern = "^f2-rule-(table1|table2)", full.names = TRUE)
l <- unlist(lapply(fs, function(f) grep("^RULE t[12] .*truth_max", readLines(f, warn = FALSE), value = TRUE)))
g <- function(key) as.numeric(sub(paste0(".* ", key, " ([^ ]+).*"), "\\1", l))
d <- data.frame(table = substr(l, 6, 7), fam = sub("RULE t[12] ([a-z_A-Z.]+) .*", "\\1", l),
                truth = g("truth_max"), T0 = g("T0"), R1 = suppressWarnings(g("R1")), R1s = suppressWarnings(g("R1s")),
                T2 = g("T2"), R2 = g("R2"), gap = g("gap_logLik"), stringsAsFactors = FALSE)
cat("fits scored:", nrow(d), "; by table:", table(d$table), "\n")
rules <- list(
  "R0 |T| > 1" = abs(d$T0) > 1,
  "R1 gate, per parameter >= 0.5" = abs(d$T0) > 1 & !is.na(d$R1) & d$R1 >= 0.5,
  "R1 gate, per parameter >= 0.3" = abs(d$T0) > 1 & !is.na(d$R1) & d$R1 >= 0.3,
  "R1s gate, one sweep >= 0.5" = abs(d$T0) > 1 & !is.na(d$R1s) & d$R1s >= 0.5,
  "R1s gate, one sweep >= 0.3" = abs(d$T0) > 1 & !is.na(d$R1s) & d$R1s >= 0.3,
  "R2 TK V grad >= 0.5" = d$R2 >= 0.5,
  "R2 TK V grad >= 0.3" = d$R2 >= 0.3,
  "R2 TK V grad >= 0.2" = d$R2 >= 0.2,
  "TK level |T2| > 1" = abs(d$T2) > 1)
for (tb in c("t1", "t2", "all")) {
  k <- if (tb == "all") rep(TRUE, nrow(d)) else d$table == tb
  cat(sprintf("--- %s: %d fits; truth >= 1: %d; 0.5 to 1: %d; below 0.5: %d\n", tb, sum(k), sum(k & d$truth >= 1), sum(k & d$truth >= 0.5 & d$truth < 1), sum(k & d$truth < 0.5)))
  for (nm in names(rules)) {
    w <- rules[[nm]]; w[is.na(w)] <- FALSE
    cat(sprintf("%-32s hit %2d of %2d; false alarm %2d of %2d; warns on %2d of 0.5-1; missed largest shift %.2f\n", nm,
                sum(k & d$truth >= 1 & w), sum(k & d$truth >= 1), sum(k & d$truth < 0.5 & w), sum(k & d$truth < 0.5),
                sum(k & d$truth >= 0.5 & d$truth < 1 & w), max(c(0, d$truth[k & !w]))))
  }
}
cat(sprintf("TK total against integrate total: cor %.3f; median |diff| %.3f; R2 against truth: cor %.3f, median |R2 - truth| %.3f\n",
            cor(d$T2, d$T0, use = "complete.obs"), median(abs(d$T2 - d$T0), na.rm = TRUE), cor(d$R2, d$truth, use = "complete.obs"), median(abs(d$R2 - d$truth), na.rm = TRUE)))
print(round(d[d$table == "t1", c("truth", "T0", "R1", "T2", "R2")], 2)[seq_len(min(5, sum(d$table == "t1"))), ])
bad <- d[d$truth >= 0.75 & d$R2 < 0.5, ]
cat("largest misses of R2 >= 0.5:\n"); print(head(bad[order(-bad$truth), c("fam", "truth", "T0", "R1", "T2", "R2")], 8))
fa <- d[d$truth < 0.4 & d$R2 >= 0.5, ]
cat("false alarms of R2 >= 0.5:\n"); print(head(fa[order(-fa$R2), c("fam", "truth", "T0", "R1", "T2", "R2")], 8))
