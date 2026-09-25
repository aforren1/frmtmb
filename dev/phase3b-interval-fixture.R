# Punch round 2, item 3: the test fixture of Rmpfr interval masses, the
# review's 511 (dev/phase3b-review2/interval.rds) and the 240 of
# dev/phase3b-interval-ref-small.R, written at 17 significant digits.
# Usage: Rscript dev/phase3b-interval-fixture.R
k <- c("t1", "t2", "v", "a", "w", "ref")
p <- rbind(readRDS("dev/phase3b-review2/interval.rds")[, k],
           readRDS("dev/phase3b-log/interval-small.rds")[, k])
p[] <- lapply(p, function(z) sprintf("%.17g", z))
fx <- "extensions/frmtmb.eam/tests/testthat/fixtures"
utils::write.csv(p, file.path(fx, "wiener-interval-ref.csv"),
                 row.names = FALSE, quote = FALSE)
cat(nrow(p), "rows\n")
