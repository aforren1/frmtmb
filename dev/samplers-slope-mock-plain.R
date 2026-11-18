# As dev/samplers-slope-mock.R, with as_stanfit(plain = TRUE): is the
# stall in the stanfit subclass? Where the q1 big_slope run stalls: as_stanfit() or the
# brms mock backend. Each step timed; stanli runs 4 chains x 1000/1000.
source("dev/samplers-designs.R")
samplers_lib()
suppressPackageStartupMessages(library(stanli))
spec <- samplers_big("big_slope")
bargs <- list(brms::bf(y ~ x + (1 + x | g)), data = spec$data,
              family = brms::bernoulli())
code <- do.call(brms::make_stancode, bargs)
sdat <- do.call(brms::make_standata, bargs)
m <- stanli_model(code = code, data = unclass(sdat))
f <- sample_model(m, chains = 4, warmup = 1000, samples = 1000, seed = 1,
                  refresh = 0, parallel_chains = 1)
st <- function(what) cat(format(Sys.time(), "%H:%M:%S"), what, "\n")
st("sampled")
cat("columns:", length(f$columns), "\n")
sf <- as_stanfit(f)
st("as_stanfit done")
cat("stanfit pars:", length(sf@sim$fnames_oi), "\n")
d <- as.matrix(sf)
st("as.matrix done")
sf <- as_stanfit(f, plain = TRUE)
st("plain as_stanfit done")
bf <- do.call(brms::brm, c(bargs, list(backend = "mock", mock_fit = sf,
                                        rename = TRUE)))
st("mock done")
