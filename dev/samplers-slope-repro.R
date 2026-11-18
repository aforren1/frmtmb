# Reproduce the stanli stall on big_slope (dev/samplers-q1.R ran over
# 44 CPU minutes in sample_model() with no end): the q1 call with
# progress printed, one chain at a time.
source("dev/samplers-designs.R")
samplers_lib()
suppressPackageStartupMessages(library(stanli))
spec <- samplers_big("big_slope")
bargs <- list(brms::bf(y ~ x + (1 + x | g)), data = spec$data,
              family = brms::bernoulli())
code <- do.call(brms::make_stancode, bargs)
sdat <- do.call(brms::make_standata, bargs)
m <- stanli_model(code = code, data = unclass(sdat))
for (ch in 1:4) {
  t0 <- proc.time()[[3]]
  cat("chain", ch, "\n")
  f <- sample_model(m, chains = 4, warmup = 1000, samples = 1000, seed = 1,
                    refresh = 100, parallel_chains = 1)
  cat("4 chains seconds:", proc.time()[[3]] - t0, "\n")
  print(f$report)
  break
}
