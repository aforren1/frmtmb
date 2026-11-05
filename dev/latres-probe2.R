# Probe: where brms 2.23.0 adds the latent residual and its AR part on
# families whose program applies an inverse link in the model block, and
# on a term of an ordinal family. Prints the model block only.
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
d <- expand.grid(t = 1:5, g = factor(1:8))
d$x <- rnorm(nrow(d))
d$y <- rpois(nrow(d), exp(0.3 + 0.2 * d$x)) + 1
d$o <- factor(pmin(d$y, 4), ordered = TRUE)
d$b <- as.integer(d$y > 2)
mb <- function(code) {
  l <- strsplit(code, "\n")[[1]]
  i <- grep("^model \\{", l)
  j <- grep("^generated quantities", l)
  cat(l[i:(j - 1L)], sep = "\n")
}
show <- function(f, fam, ...) {
  cat("\n=====", deparse(f), fam$family, fam$link, "\n")
  r <- try(stancode(f, data = d, family = fam, ...))
  if (!inherits(r, "try-error")) mb(r)
}
show(y ~ x + ar(t, g), Gamma(link = "log"))
show(y ~ x + ar(t, g), Gamma())
show(y ~ x + ar(t, g), poisson("sqrt"))
show(b ~ x + ar(t, g), bernoulli("probit"))
show(b ~ x + ar(t, g, cov = TRUE), bernoulli("probit"))
show(y ~ x + ar(t, g, cov = TRUE), Gamma())
show(y ~ x + ar(t, g), lognormal())
show(y ~ x + ar(t, g), zero_inflated_poisson())
show(o ~ x + ar(t, g), cumulative())
show(o ~ x + ar(t, g, cov = TRUE), cumulative())
show(o ~ x + cs(x) + ar(t, g, cov = TRUE), sratio())
show(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x), lognormal())
show(bf(y ~ x + ar(t, g), sigma ~ x), lognormal())
show(y | weights(x + 3) ~ x + ar(t, g, cov = TRUE), poisson())
show(y | cens(b) ~ x + ar(t, g, cov = TRUE), poisson())
show(y | trunc(lb = 0) ~ x + ar(t, g, cov = TRUE), poisson())
show(bf(y ~ a + x, a ~ 1, nl = TRUE) + acformula(~ar(t, g, cov = TRUE)),
     poisson(), prior = set_prior("normal(0,1)", nlpar = "a"))
show(y ~ x + (1 | g) + ar(t, g, cov = TRUE), poisson())
show(y ~ x + ar(g = g), poisson())
