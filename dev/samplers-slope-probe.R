# Why stanli is slow on big_slope: gradient time and one short chain,
# against rstan's gradient time on the same brms program.
source("dev/samplers-designs.R")
samplers_lib()
suppressPackageStartupMessages(library(stanli))
spec <- samplers_big("big_slope")
d <- spec$data
code <- brms::make_stancode(y ~ x + (1 + x | g), data = d, family = brms::bernoulli())
sdat <- brms::make_standata(y ~ x + (1 + x | g), data = d, family = brms::bernoulli())
m <- stanli_model(code = code, data = unclass(sdat))
set.seed(1)
q <- rnorm(length(unconstrain(m, list(), partial = TRUE)), 0, 0.5)
t0 <- proc.time()[[3]]; for (i in 1:200) log_prob_grad(m, q); tg <- (proc.time()[[3]] - t0) / 200
cat("stanli gradient seconds:", tg, "npar", length(q), "\n")
t0 <- proc.time()[[3]]
f <- sample_model(m, chains = 1, warmup = 100, samples = 100, seed = 1, refresh = 20, parallel_chains = 1)
cat("stanli 1 chain 100/100 seconds:", proc.time()[[3]] - t0, "\n")
print(f$report)
print(table(f$sampler[, , "treedepth__"]))
