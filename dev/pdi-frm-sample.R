# dev/pdi-frm-sample.R -- the end-to-end symptom, at the frmtmb level.
#
# Fit the same model twice with frm_sample(), once under the default priors
# and once under priors = "flat", and compare the draws. They must differ.
# Run this against the stock container tmbstan and against a tmbstan whose
# reverse-mode log_prob was patched (see dev/prior-dropping-investigation.md).
#
#   Rscript dev/pdi-frm-sample.R

options(warn = 1, digits = 12)
suppressPackageStartupMessages(library(frmtmb))

cat("##### PDI frm_sample symptom\n")
cat("tmbstan lib: ", find.package("tmbstan"), "\n", sep = "")
cat("stan version: ", rstan::stan_version(), "\n", sep = "")

set.seed(1)
n <- 40L
d <- data.frame(x = rnorm(n), g = factor(rep(1:8, each = 5)))
d$y <- 1 + 0.8 * d$x + rep(rnorm(8, 0, 0.7), each = 5) + rnorm(n, 0, 0.5)

run <- function(pr) {
  set.seed(99)
  s <- frm_sample(y ~ x + (1 | g), data = d, family = gaussian(),
                  priors = pr, chains = 1, iter = 800, warmup = 400,
                  seed = 42, refresh = 0)
  as.matrix(s)
}

a <- run(NULL)      # brms-style defaults
b <- run("flat")

keep <- intersect(colnames(a), colnames(b))
keep <- keep[!grepl("^b\\[", keep)]
cat("\nparameter        default_sd        flat_sd     max|default-flat|\n")
for (k in keep) {
  cat(sprintf("%-14s %14.8f %14.8f %18.8f\n", k, sd(a[, k]), sd(b[, k]),
              max(abs(a[, k] - b[, k]))))
}
cat("\nidentical(default, flat) = ", identical(a[, keep], b[, keep]), "\n",
    sep = "")
cat("max |default - flat| over all shared columns = ",
    max(abs(a[, keep] - b[, keep])), "\n", sep = "")
cat("\n##### done\n")
