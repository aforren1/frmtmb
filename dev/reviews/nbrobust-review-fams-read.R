# Review of lane nbrobust: the 0.83.0 rows of the family check whose
# value is off by more than 1e-9 on or off the tape, beside the lane's
b <- read.delim("dev/reviews/nbrobust-review-fams-base.tsv")
l <- read.delim("dev/reviews/nbrobust-review-fams-lib.tsv")
stopifnot(identical(b$mu, l$mu), identical(b$lsh, l$lsh))
w <- which(!(b$absv <= 1e-9) | !(b$absvoff <= 1e-9) | is.na(b$absvoff) | is.na(b$absv))
out <- data.frame(fam = b$fam, link = b$link, rate = b$rate,
                  mu = signif(b$mu, 4), lsh = signif(b$lsh, 4),
                  vref = signif(b$vref, 8),
                  base_tape = signif(b$absv, 3), base_off = signif(b$absvoff, 3),
                  lane_tape = signif(l$absv, 3), lane_off = signif(l$absvoff, 3))
print(out[w, ], row.names = FALSE)
cat("lane: max abs value error on the tape", max(l$absv), " off", max(l$absvoff),
    "\n")
