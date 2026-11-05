# mi() on a latent response beside a latent residual, once through frm_sample(), on
# the merged build: frm() refuses the pair, and the sampler holds every value.
# Seed 12301. brms 2.23.0 reads the same model.
.libPaths(c("C:/Users/adf44/source/r/rellib-r17",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-rel079/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
set.seed(12301)
d <- expand.grid(t = 1:8, g = factor(1:10))
d <- d[order(d$g, d$t), ]
d$x <- rnorm(nrow(d))
d$y <- rgamma(nrow(d), 5, 5 / exp(0.3 * d$x))
d$y[c(3L, 20L, 41L)] <- NA
r <- tryCatch(frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)), family = Gamma(link = "log"), data = d),
              error = function(e) conditionMessage(e))
cat("frm:", substr(as.character(r)[1], 1, 180), "\n")
t0 <- proc.time()[[3]]
ds <- frm_sample(bf(y | mi() ~ x + ar(t, g, cov = TRUE)), family = Gamma(link = "log"), data = d,
                 chains = 1, iter = 400, refresh = 0, seed = 19)
cat("frm_sample ok in", round(proc.time()[[3]] - t0, 1), "s\n")
m <- as.matrix(ds$draws)
print(colnames(m)[grepl("Ymi|sderr|ar\\[", colnames(m))])
print(round(colMeans(m[, grepl("Ymi|^sderr|^ar\\[", colnames(m)), drop = FALSE]), 3))
ep <- posterior_epred(ds)
cat("epred", dim(ep), "any NA:", anyNA(ep), "\n")
