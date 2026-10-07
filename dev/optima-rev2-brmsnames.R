# Reviewer of lane optima, re-check (d): brms 2.23.0's simo names for
# mo() in sigma and in an nlpar (dev/optima-rev2-readers.R's data).
#   Rscript dev/optima-rev2-brmsnames.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
set.seed(5)
n <- 300
d <- data.frame(x1 = sample(0:3, n, TRUE), x2 = sample(0:4, n, TRUE),
                z = rnorm(n), g = factor(sample(1:15, n, TRUE)))
u <- rnorm(15, sd = 0.5)
d$y <- c(0, 1, 1, 2)[d$x1 + 1] + c(0, 0.3, 0.6, 0.6, 1)[d$x2 + 1] +
  0.3 * d$z + u[d$g] + rnorm(n, sd = exp(c(0, 0.2, 0.2, 0.4)[d$x1 + 1]))
d$yn <- 2 * exp(-0.2 * c(0, 1, 1, 3)[d$x1 + 1]) * (1 + 0.3 * d$z) +
  rnorm(n, sd = 0.3)
for (f in list(bf(y ~ z, sigma ~ mo(x1)),
               bf(yn ~ a * (1 + bz * z), a ~ mo(x1), bz ~ 1, nl = TRUE))) {
  b <- suppressWarnings(brm(f, data = d, chains = 1, iter = 200, seed = 1,
                            refresh = 0,
                            prior = if (!is.null(f$pforms$a))
                              c(prior(normal(2, 1), nlpar = a),
                                prior(normal(0, 1), nlpar = bz))))
  cat(grep("simo|^bsp", variables(b), value = TRUE), "\n")
}
