# Reviewer probe 4b: what the flat run's prior_summary() resolves to when
# passed back. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
dd <- re_data(1L)
q <- function(expr) suppressWarnings(suppressMessages(expr))
fit2 <- q(frm(bf(y ~ x + (1 + x + z | g)), data = dd))
ds <- q(frm_sample(fit2, prior = "flat", chains = 1, iter = 300, warmup = 150,
                   seed = 1, refresh = 0))
ps <- prior_summary(ds)
cat("class of prior_summary():", class(ps), "\n")
show_entries(q(frmtmb.sample:::sample_resolve_priors(fit2, "flat"))$ri)
cat("-- passed back\n")
show_entries(q(frmtmb.sample:::sample_resolve_priors(fit2, ps))$ri)
