# Bisect the stanli stall on big_slope. Arg "fit": frm() fits the model
# in the same process before stanli samples; "load": frmtmb is loaded
# only; "sample": frmtmb.sample (and so rstan) is loaded; "eq": the
# survey instrument is loaded as dev/samplers-q1.R loads it; "none": none. stanli then samples 4 chains x 1000/1000 with
# the refresh given as the second argument (default 500).
mode <- commandArgs(trailingOnly = TRUE)[1]
refresh <- as.integer(commandArgs(trailingOnly = TRUE)[2] %||% "500")
source("dev/samplers-designs.R")
samplers_lib()
spec <- samplers_big("big_slope")
if (mode %in% c("fit", "load")) suppressPackageStartupMessages(library(frmtmb))
if (mode == "sample") loadNamespace("frmtmb.sample")
if (mode == "eq") eq <- samplers_eq_env()
if (mode == "fit") {
  f <- frm(y ~ x + (1 + x | g), data = spec$data, family = bernoulli())
  cat("frm done\n")
}
suppressPackageStartupMessages(library(stanli))
bargs <- list(brms::bf(y ~ x + (1 + x | g)), data = spec$data,
              family = brms::bernoulli())
code <- do.call(brms::make_stancode, bargs)
sdat <- do.call(brms::make_standata, bargs)
m <- stanli_model(code = code, data = unclass(sdat))
cat(format(Sys.time(), "%H:%M:%S"), "model ready\n")
t0 <- proc.time()[[3]]
f <- sample_model(m, chains = 4, warmup = 1000, samples = 1000, seed = 1,
                  refresh = refresh, parallel_chains = 1)
cat(mode, "4 chains seconds:", proc.time()[[3]] - t0, "\n")
