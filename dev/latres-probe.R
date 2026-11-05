# Probe: brms 2.23.0's Stan code for autocor terms on a poisson response.
# Prints the program for a few terms so the latent-residual model can be
# read off brms's own code.
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
d <- expand.grid(t = 1:5, g = factor(1:8))
d$x <- rnorm(nrow(d))
d$y <- rpois(nrow(d), exp(0.3 + 0.2 * d$x))
M <- diag(nrow(d))
for (f in c("y ~ x + ar(t, g, cov = TRUE)", "y ~ x + ar(t, g)",
            "y ~ x + ar(t, g, p = 2)",
            "y ~ x + arma(t, g, cov = TRUE)", "y ~ x + cosy(t, g)",
            "y ~ x + unstr(t, g)", "y ~ x + fcor(M)")) {
  cat("\n\n=====", f, "\n")
  r <- try(stancode(as.formula(f), data = d, family = poisson(),
                    data2 = list(M = M)))
  if (!inherits(r, "try-error")) print(r)
}
cat("\n\n===== default_prior ar(cov = TRUE)\n")
print(default_prior(y ~ x + ar(t, g, cov = TRUE), data = d,
                    family = poisson()))
cat("\n\n===== default_prior ar(cov = FALSE)\n")
print(default_prior(y ~ x + ar(t, g), data = d, family = poisson()))
