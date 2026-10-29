# Lane meme: brms 2.23.0's prior-draw names for me() hyperparameters
# with coefficient-specific priors, after renaming (one small fit).
# Output: dev/meme-out/brms-names.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
set.seed(11)
n <- 40
d <- data.frame(x = rnorm(n), z = rnorm(n), sx = 0.3,
                sz = runif(n, 0.2, 0.4), g = factor(rep(1:8, each = 5)))
d$xg <- rnorm(8)[d$g]
d$sg <- 0.25
d$y <- 1 + 0.5 * d$x + rnorm(n)
pr <- c(prior(normal(0, 2), class = meanme),
        prior(normal(1, 0.5), class = meanme, coef = mez),
        prior(exponential(1), class = sdme),
        prior(lkj(2), class = corme))
f1 <- suppressWarnings(suppressMessages(brm(
  bf(y ~ me(x, sx) + me(z, sz) + me(xg, sg, gr = g)), data = d,
  prior = pr, sample_prior = "yes", chains = 1, iter = 200, seed = 1,
  refresh = 0)))
cat("variables:", grep("^Xme", variables(f1), value = TRUE, invert = TRUE),
    "\n", fill = 78)
cat("prior_draws:", names(prior_draws(f1)), "\n")
print(as.data.frame(prior_summary(f1))[, c("prior", "class", "coef",
                                           "group", "lb", "ub",
                                           "source")])
cat("\nwith gr = g correlated and a group corme prior:\n")
d$zg <- rnorm(8)[d$g]
f2 <- suppressWarnings(suppressMessages(brm(
  bf(y ~ me(xg, sg, gr = g) + me(zg, sg, gr = g)), data = d,
  prior = c(prior(normal(0, 2), class = meanme),
            prior(exponential(1), class = sdme),
            prior(lkj(3), class = corme, group = g)),
  sample_prior = "yes", chains = 1, iter = 200, seed = 1, refresh = 0)))
cat("variables:", grep("^Xme", variables(f2), value = TRUE, invert = TRUE),
    "\n", fill = 78)
cat("prior_draws:", names(prior_draws(f2)), "\n")
print(as.data.frame(prior_summary(f2))[, c("prior", "class", "coef",
                                           "group", "lb", "ub",
                                           "source")])
