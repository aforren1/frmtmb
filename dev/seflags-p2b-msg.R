# Lane seflags, punch round 2b: the ill-conditioning warning as a user
# sees it, degree 8 without random effects, seed 7.
#   Rscript dev/seflags-p2b-msg.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
set.seed(7)
d <- data.frame(x = stats::runif(200, 1, 2))
d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.2)
for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
f <- frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, data = d)
print(fixef(f))
