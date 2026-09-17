## Punch round 2, MINOR 2: the scale tier's switch to bernoulli() fits the
## model base fitted as binomial() without trials. Construction (a) of
## dev/famlink-rev2-ext-trials.R (seed 20260908, n = 200), whose base arm
## logged logLik -112.7997206422. Also the brms-coexistence line with brms
## attached, where `bernoulli` is brms's constructor.
## Usage: Rscript dev/famlink-p2-scale-bern.R
ARM <- "lane"
source("dev/famlink-rev-common.R")
set.seed(20260908L)
n <- 200L
d <- data.frame(x = rnorm(n), g1 = factor(rep_len(1:20, n)),
                g2 = factor(sample.int(10, n, TRUE)))
d$y <- rbinom(n, 1L, plogis(-0.5 + 0.8 * d$x +
                              rnorm(20, 0, 0.5)[as.integer(d$g1)]))
f <- frm(bf(y ~ x + (1 | g1) + (1 | g2)), family = bernoulli(), data = d)
cat(sprintf("(a) bernoulli(): logLik %.10f\n", logLik(f)))
suppressMessages(library(brms))
set.seed(1234)
dc <- data.frame(x = rnorm(80), g = gl(8, 10))
dc$y <- 1 + 2 * dc$x + rnorm(80)
dc$z <- rbinom(80, 1, plogis(0.5 * dc$x))
fb <- frm(z ~ x + (1 | g), data = dc, family = bernoulli())
cat(sprintf("coexistence: bernoulli is %s; fit logLik %.10f, family %s\n",
            environmentName(environment(bernoulli)), logLik(fb),
            family(fb)$family))
