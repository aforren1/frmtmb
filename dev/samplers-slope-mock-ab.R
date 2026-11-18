# Is brm(backend = "mock") slow on big_slope for any stanfit, or only for
# stanli's? Short chains (4 x 100/100) from each sampler, then the mock
# backend on each, rstan's first. Each step timestamped.
source("dev/samplers-designs.R")
samplers_lib()
suppressPackageStartupMessages(library(stanli))
st <- function(what) cat(format(Sys.time(), "%H:%M:%S"), what, "\n")
spec <- samplers_big("big_slope")
bargs <- list(brms::bf(y ~ x + (1 + x | g)), data = spec$data,
              family = brms::bernoulli())
code <- do.call(brms::make_stancode, bargs)
sdat <- do.call(brms::make_standata, bargs)
mock <- function(sf) do.call(brms::brm, c(bargs, list(
  backend = "mock", mock_fit = sf, rename = TRUE)))
rm_ <- rstan::stan_model(model_code = code)
rf <- rstan::sampling(rm_, data = sdat, chains = 4, iter = 200, seed = 1,
                      refresh = 0, cores = 1)
st("rstan sampled")
b1 <- mock(rf)
st("rstan mock done")
m <- stanli_model(code = code, data = unclass(sdat))
f <- sample_model(m, chains = 4, warmup = 100, samples = 100, seed = 1,
                  refresh = 0, parallel_chains = 1)
st("stanli sampled")
b2 <- mock(as_stanfit(f))
st("stanli mock done")
