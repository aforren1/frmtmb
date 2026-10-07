# Reviewer of lane surface: brms 2.23.0 fits for claims 1 and 4.
#   Rscript dev/surface-rev-brms.R
# Writes dev/surface-rev-out/brms-ppmix.rds and brms-multiple.rds.
# Compiles two Stan programs (rstan backend).
source("dev/surface-rev-env.R")
rev_env("lane")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")

## claim 1: pp_mixture on a two-gaussian mixture, the lane's MC data
set.seed(4)
dm <- data.frame(y = c(rnorm(150, -1.5), rnorm(150, 1.5)))
t0 <- proc.time()
if (!file.exists("dev/surface-rev-out/brms-ppmix.rds")) {
bm <- brm(y ~ 1, family = mixture(gaussian, gaussian), data = dm,
          chains = 4, iter = 3000, seed = 1, refresh = 0,
          backend = "rstan")
cat("mixture fit", (proc.time() - t0)[3], "s\n")
print(summary(bm))
pm <- pp_mixture(bm)
cat("brms pp_mixture: class", class(pm), "dim", dim(pm), "\n")
str(dimnames(pm))
saveRDS(list(pm = pm, data = dm, summ = posterior_summary(bm)),
        "dev/surface-rev-out/brms-ppmix.rds")
}

## claim 4: brm_multiple on nhanes, m = 5, mice seed 1
data("nhanes", package = "mice")
imp <- mice::mice(nhanes, m = 5, print = FALSE, seed = 1)
t0 <- proc.time()
bmu <- brm_multiple(bmi ~ age * chl, data = imp, chains = 2, iter = 4000,
                    seed = 1, refresh = 0, backend = "rstan")
cat("brm_multiple", (proc.time() - t0)[3], "s\n")
print(summary(bmu))
cat("rhats:\n"); print(bmu$rhats)
cat("nrow(bmu$data)", nrow(bmu$data), "| data equals complete(imp, 1):",
    isTRUE(all.equal(bmu$data[, c("bmi", "age", "chl")],
                     mice::complete(imp, 1)[, c("bmi", "age", "chl")],
                     check.attributes = FALSE)), "\n")
ce <- conditional_effects(bmu, "chl")[[1]]
saveRDS(list(fixef = fixef(bmu), summ = posterior_summary(bmu),
             ce = ce, data1 = bmu$data,
             draws = as.data.frame(bmu)),
        "dev/surface-rev-out/brms-multiple.rds")
cat("DONE\n")
