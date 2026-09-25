# Is the sampler's prior the same function on both arms for mix2_dpar,
# or only the same at the estimate? Evaluate the negative log prior of
# the resolved sample-route entries at 500 random points; save for a
# cross-arm identical() and ulp comparison. Also print entry order.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
suppressPackageStartupMessages(library(frmtmb.sample))
source(file.path(R2_ROOT, "dev/mvprior-review2/r2-cases.R"))
d <- r2_data()
ns <- asNamespace("frmtmb")
pr <- set_prior("normal(0, 1)", class = "b", dpar = "mu1") +
  set_prior("normal(-2, 2)", class = "Intercept", dpar = "mu1") +
  set_prior("normal(3, 2)", class = "Intercept", dpar = "mu2") +
  set_prior("student_t(3, 0, 2.5)", class = "sigma1")
uf <- suppressMessages(frm(bf(ym ~ x), family = mixture(gaussian(), gaussian()),
                           data = d, dry_run = "objective"))
fit <- suppressWarnings(suppressMessages(frm(bf(ym ~ x), family = mixture(gaussian(), gaussian()), data = d)))
ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(uf, pr))$ri
for (e in ri$entries) cat(e$comp, paste(e$idx, collapse = ","),
                          paste(unlist(e$dist), collapse = " "), "\n")
f <- ns$neg_log_prior_fn(ri$entries)
est <- fit$estimates
set.seed(99)
vals <- vapply(1:500, function(i) {
  p <- lapply(est, function(v) v + rnorm(length(v), 0, 1))
  f(p)
}, 0)
cat("nlp at estimate", format(f(est), digits = 17), "\n")
saveRDS(vals, file.path(R2_ROOT, "dev/mvprior-review2", paste0("r2-nlp-", r2_arm, ".rds")))
