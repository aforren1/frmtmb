# Punch round 2, item 3: the interval mass against the review's 511
# Rmpfr intervals (dev/phase3b-review2/interval.rds, whose `ref` column
# is the 500-bit log mass and `refF1` the log F(t1)).
# Usage: Rscript dev/phase3b-interval-check.R <lib> [points.rds out.txt]
# Output: dev/phase3b-log/interval-check.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb.eam")
p <- readRDS(if (length(a) > 1) a[2] else "dev/phase3b-review2/interval.rds")
out_file <- if (length(a) > 2) a[3] else "dev/phase3b-log/interval-check.txt"
got <- ns$ddm_rt_linterval_b(p$t1, p$t2, p$v, p$a, p$w, 0)
err <- abs(got - p$ref)
b <- cut(log10(p$mass_over_F1), c(-Inf, -16, -12, -8, -4, 0.5))
out <- c(sprintf("package: %s", find.package("frmtmb.eam")),
         sprintf("points: %d; abs error of the log mass: max %.3g, median %.3g",
                 nrow(p), max(err), median(err)),
         "by mass / F(t1) band (max abs log error, count):",
         paste(capture.output(print(rbind(max = tapply(err, b, max),
                                          n = table(b)))), collapse = "\n"))
w <- order(-err)[1:6]
out <- c(out, "worst 6:",
         paste(capture.output(print(cbind(p[w, c("v", "a", "w", "u1", "r",
                                                  "ref")], got = got[w],
                                          err = err[w],
                                          mass_over_F1 = p$mass_over_F1[w]),
                                    digits = 5)), collapse = "\n"))
# the tape: gradient finite over the points
tp <- RTMB::MakeTape(function(q) {
  sum(ns$ddm_rt_linterval_b(p$t1, p$t2, q[1] + 0 * p$v + p$v, exp(q[2]) * p$a,
                            p$w, 0))
}, c(0, 0))
out <- c(out, sprintf("gradient finite at the points: %s",
                      all(is.finite(tp$jacobian(c(0, 0))))))
writeLines(out, out_file)
cat(out, sep = "\n")
