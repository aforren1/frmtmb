# Reviewer of lane optima, task 3: brms 2.23.0 with its default priors
# (simo ~ dirichlet(1)) on the data of optima-rev-sample2-arm.R.
#   Rscript dev/optima-rev-sample3-brms.R   (after sample2 wrote data.rds)
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "rstan",
    format(packageVersion("rstan")), "\n")
out <- "C:/Users/adf44/source/r/frmtmb-wt-optima/dev/optima-rev-out"
# same construction as sample2, so this need not wait for data.rds
set.seed(1); n <- 100
x <- sample(0:3, n, TRUE)
weak <- data.frame(x = x, y = 0.05 * x + rnorm(n))
set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
strong <- data.frame(income, ls)

res <- list()
for (nm in c("weak", "strong")) {
  f <- if (nm == "weak") y ~ mo(x) else ls ~ mo(income)
  d <- if (nm == "weak") weak else strong
  print(get_prior(f, data = d))
  t0 <- Sys.time()
  b <- brm(f, data = d, chains = 2, iter = 2000, seed = 1, cores = 1,
           refresh = 0)
  cat("brm time", format(Sys.time() - t0), "\n")
  print(prior_summary(b))
  dm <- as.matrix(b)
  sm <- rstan::summary(b$fit)$summary
  print(round(sm[, c("mean", "sd", "n_eff", "Rhat")], 4))
  sp <- rstan::get_sampler_params(b$fit, inc_warmup = FALSE)
  cat("divergent per chain:", sapply(sp, function(x) sum(x[, "divergent__"])),
      "\n")
  res[[nm]] <- list(draws = dm, summary = sm,
                    div = sapply(sp, function(x) sum(x[, "divergent__"])))
}
saveRDS(res, file.path(out, "draws-brms.rds"))
