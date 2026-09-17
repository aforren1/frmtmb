## Reviewer recheck: the point evidence ratio on draws, base against lane,
## on the gauss_prior model of dev/brmsnames-rev2-natural.R (same data and
## per-model seed, so the same draws).
##   Rscript dev/brmsnames-rev2-er.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
src <- readLines("dev/brmsnames-rev2-natural.R")
a <- grep("^if [(]arm != \"compare\"[)]", src)
b <- grep("if [(]length[(]only[)][)]", src)
body <- c(src[1:(a - 1L)], src[(grep("q[(]library[(]frmtmb[)][)]", src)[1]):b])
only <- "gauss_prior"
eval(parse(text = body))
M <- ms$gauss_prior
fit <- q(frm(M$f, family = M$fam, data = d))
set.seed(sum(utf8ToInt("gauss_prior")))
ds <- q(frm_sample(fit, chains = 2, iter = 200, refresh = 0, seed = 9,
                   prior = M$prior))
h <- q(hypothesis(ds, c("x = 0", "x - 0.3 = 0", "x + Intercept = 0")))
print(if (is.data.frame(h)) h else h$hypothesis, digits = 10)
