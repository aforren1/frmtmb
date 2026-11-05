# The package's rule against the reviewer's exact shifts. A hit warns at a shift of
# 1 SE or more, a false alarm warns below 0.5 SE (dev/reviews/2026-11-02-refuse.md,
# RF1: 75 of 76 hits, 3 of 39 false alarms, largest silent 1.01).
#   Rscript dev/rel079-rf1-score.R
O <- "C:/Users/adf44/source/r/frmtmb-wt-refuse/dev/review-refuse-out/f2"
fs <- list.files(O, pattern = "^f2-rule-(table1|table2)", full.names = TRUE)
a <- unlist(lapply(fs, function(f) grep("^RULE t[12] .*truth_max", readLines(f, warn = FALSE), value = TRUE)))
g <- function(l, key) suppressWarnings(as.numeric(sub(paste0(".* ", key, " ([^ ]+).*"), "\\1", l)))
d <- data.frame(tag = sub("^RULE (t[12] .*) truth_max.*", "\\1", a), truth = g(a, "truth_max"))
L <- "dev/rel079-log/rr/rf1"
b <- unlist(lapply(list.files(L, pattern = "log$", full.names = TRUE), function(f) {
  grep("^PKG t[12] ", readLines(f, warn = FALSE), value = TRUE)
}))
e <- data.frame(tag = sub("^PKG (t[12] .*) tk .*", "\\1", b), warn = grepl("warn TRUE", b),
                mv = g(b, "mv"), tk = g(b, "tk"))
m <- merge(d, e, by = "tag")
cat("joined", nrow(m), "of", nrow(d), "truths and", nrow(e), "fits\n")
hit <- m$truth >= 1
fa <- m$truth < 0.5
cat(sprintf("hit %d of %d; false alarm %d of %d; largest silent %.2f\n", sum(m$warn[hit]), sum(hit),
            sum(m$warn[fa]), sum(fa), max(m$truth[!m$warn])))
tb1 <- startsWith(m$tag, "t1")
cat(sprintf("false alarm, table 1: %d of %d\n", sum(m$warn[fa & tb1]), sum(fa & tb1)))
print(m[(hit & !m$warn) | (fa & m$warn), ], row.names = FALSE)
