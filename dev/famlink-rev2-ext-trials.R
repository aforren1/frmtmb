## Recheck, priority 3: extension constructions reaching the trials refusal.
## (a) frmtmb.sample's gated scale tier (test-scale.R line 54, FRMTMB_SCALE_TESTS)
##     at its small size: bf(y ~ x + (1 | g1) + (1 | g2)), stats::binomial(), 0/1 y.
## (b) frmtmb.latent hmm(K = 2, binomial()) on 0/1 data without trials().
## Seeds in-line. Usage: Rscript dev/famlink-rev2-ext-trials.R <lane|base>
ARM <- commandArgs(trailingOnly = TRUE)[1]
source("dev/famlink-rev-common.R")
suppressMessages(library(frmtmb.latent))
set.seed(20260908L)
n <- 200L
d <- data.frame(x = rnorm(n), g1 = factor(rep_len(1:20, n)), g2 = factor(sample.int(10, n, TRUE)))
d$y <- rbinom(n, 1L, plogis(-0.5 + 0.8 * d$x + rnorm(20, 0, 0.5)[as.integer(d$g1)]))
r <- tryCatch({ f <- frm(bf(y ~ x + (1 | g1) + (1 | g2)), family = stats::binomial(), data = d); sprintf("fit, logLik %.10f", logLik(f)) },
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("(a) scale-tier binomial:", r, "\n")
set.seed(20260916)
T <- 300; st <- numeric(T); st[1] <- 1
for (t in 2:T) st[t] <- if (runif(1) < 0.9) st[t - 1] else 3 - st[t - 1]
dh <- data.frame(t = seq_len(T), y = rbinom(T, 1, c(0.2, 0.8)[st]))
r2 <- tryCatch({ f <- frm(y ~ 1, family = hmm(K = 2, binomial(), time = t), data = dh); sprintf("fit, logLik %.10f", logLik(f)) },
               error = function(e) paste("ERROR:", conditionMessage(e)))
cat("(b) hmm(K = 2, binomial()) no trials:", r2, "\n")
