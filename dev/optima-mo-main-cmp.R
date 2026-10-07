# Lane optima: the plain ls ~ mo(income) fits, lane against base, paired
# by seed (dev/optima-mo-main.R output).
a <- utils::read.delim("dev/optima-log/mo-main-base.tsv")
b <- utils::read.delim("dev/optima-log/mo-main-lane3.tsv")
m <- merge(a, b, by = "seed")
d <- m$ll_fit.y - m$ll_fit.x
cat("plain mo(): lane - base logLik over", nrow(m), "seeds: max",
    format(max(d), digits = 3), "min", format(min(d), digits = 3),
    "; |difference| > 1e-6:", sum(abs(d) > 1e-6), "\n")
