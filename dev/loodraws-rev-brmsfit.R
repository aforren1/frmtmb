# Reviewer of lane loodraws: brms 2.23.0 fits for the same-draws
# comparison (dev/loodraws-rev-same.R). Data seed 3, sampler seed 1.
#   Rscript dev/loodraws-rev-brmsfit.R
# Writes dev/loodraws-rev-out/brmsfits.rds; log
# dev/loodraws-rev-out/brmsfit.txt.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/loodraws-rev-out"
dir.create(out_dir, showWarnings = FALSE)
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "loo",
    format(packageVersion("loo")), "\n")
set.seed(3)
n <- 80
dd <- data.frame(x = rnorm(n), g = factor(rep(1:8, 10)))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(8, 0, 0.7)[dd$g])
dd$y2 <- rnorm(n, -1 + 0.3 * dd$x)
dd$y01 <- rbinom(n, 1, plogis(-0.3 + 0.8 * dd$x))
dd$nt <- sample(5:15, n, TRUE)
dd$s <- rbinom(n, dd$nt, plogis(0.2 + 0.5 * dd$x))
dd$w <- runif(n, 0.5, 2)
q <- function(e) suppressWarnings(suppressMessages(e))
fit1 <- function(form, fam) {
  t0 <- proc.time()[["elapsed"]]
  f <- q(brm(form, data = dd, family = fam, chains = 4, iter = 2000,
             seed = 1, refresh = 0, cores = 1))
  cat(deparse(form)[1], "s:", round(proc.time()[["elapsed"]] - t0, 1), "\n")
  f
}
fits <- list(
  gauss = fit1(bf(y ~ x), gaussian()),
  wt = fit1(bf(y | weights(w) ~ x), gaussian()),
  bern = fit1(bf(y01 ~ x), bernoulli()),
  binom = fit1(bf(s | trials(nt) ~ x), binomial()),
  ri = fit1(bf(y ~ x + (1 | g)), gaussian()),
  mv = fit1(bf(y ~ x) + bf(y2 ~ x) + set_rescor(FALSE), gaussian()))
fits$thin <- q(update(fits$gauss, thin = 2, recompile = FALSE, seed = 1,
                      refresh = 0))
cat("thin ndraws", ndraws(fits$thin), "chains", nchains(fits$thin), "\n")
saveRDS(list(dd = dd, fits = fits), file.path(out_dir, "brmsfits.rds"))
cat("done\n")
