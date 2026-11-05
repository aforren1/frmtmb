# Is the bernoulli osa defect older than the lane? REV_ARM=base.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
set.seed(1)
d <- data.frame(x = rnorm(200), g = factor(rep(1:20, each = 10)))
d$bin <- rbinom(200, 1, plogis(.2 + .4 * d$x))
f <- frm(bin ~ x, family = bernoulli(), data = d)
r <- as.matrix(residuals(f, type = "osa"))[, 1]
cat(rev_arm, "bernoulli plain osa: mean", mean(r), " sd", sd(r), "\n")
f <- frm(bin ~ x + (1 | g), family = bernoulli(), data = d)
r <- as.matrix(residuals(f, type = "osa"))[, 1]
cat(rev_arm, "bernoulli RE osa: mean", mean(r), " sd", sd(r), "\n")
d$tr <- 1L
f <- frm(bin | trials(tr) ~ x, family = binomial(), data = d)
r <- as.matrix(residuals(f, type = "osa"))[, 1]
cat(rev_arm, "binomial 1 trial osa: mean", mean(r), " sd", sd(r), "\n")
print(packageVersion("frmtmb"))
