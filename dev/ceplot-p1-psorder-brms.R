# Lane ceplot punch 1: brms 2.23.0's variables() order of the population
# coefficients on nonlinear models, to pin the order
# posterior_samples(pars = "^b_") must follow (review m5).
# brms at algorithm = "fixed_param", one chain of one iteration.
#   Rscript dev/ceplot-p1-psorder-brms.R > dev/ceplot-log/p1-psorder-brms.txt
# Data seed 11.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
set.seed(11)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n))
d$yp <- 2 * exp(0.3 * d$x) + 0.2 * d$z + rnorm(n, 0, exp(0.1 * d$w))
qb <- function(e) suppressMessages(suppressWarnings(e))
fits <- list(
  nl = bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
  nl_sigma = bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1 + w, sigma ~ w,
                nl = TRUE),
  nl_three = bf(yp ~ a * exp(b * x) + c, a ~ 1 + z, b ~ 1, c ~ 0 + w,
                nl = TRUE),
  lin_sigma = bf(yp ~ x + z, sigma ~ w)
)
for (nm in names(fits)) {
  b <- qb(brm(fits[[nm]], data = d, algorithm = "fixed_param", chains = 1,
              iter = 1, warmup = 0, refresh = 0, seed = 1, silent = 2,
              init = 0))
  cat(nm, ":", grep("^b_", variables(b), value = TRUE), "\n")
}
cat("done\n")
