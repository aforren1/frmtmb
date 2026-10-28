# What brms 2.23.0 generates for sample_prior = "yes" / "only": the
# generated quantities blocks of a set of models chosen to cover each
# prior class. No compile; stancode() only.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)),
                m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(60)
d$yb <- rbinom(60, 1, 0.5)
d$yp <- rbeta(60, 2, 2)
d$yc <- rpois(60, 3)
d$t <- rep(1:10, 6)
show <- function(lab, ...) {
  cat("\n=====", lab, "=====\n")
  sc <- tryCatch(stancode(...), error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n"); NULL
  }, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n"); NULL
  })
  if (is.null(sc)) return(invisible())
  sc <- strsplit(sc, "\n")[[1]]
  i <- grep("^generated quantities", sc)
  pr <- grep("prior_|lprior|_rng", sc, value = TRUE)
  cat(if (length(i)) sc[i:length(sc)] else "(no generated quantities)",
      sep = "\n")
  cat("--- lprior / prior lines:\n")
  cat(grep("lprior", sc, value = TRUE), sep = "\n")
}
show("A default gaussian RE", y ~ x + (1 | g), data = d,
     sample_prior = "yes")
show("B proper b, cor", y ~ x + z + (1 + x | g), data = d,
     prior = c(prior(normal(0, 1), class = b),
               prior(normal(0, 2), class = b, coef = z),
               prior(normal(1, 3), class = Intercept),
               prior(exponential(2), class = sigma),
               prior(lkj(2), class = cor),
               prior(cauchy(0, 1), class = sd, group = g, coef = x)),
     sample_prior = "yes")
show("C bounded b", y ~ x, data = d,
     prior = c(prior(normal(0, 1), class = b, lb = 0),
               prior(normal(0, 1), class = b, coef = x, ub = 2)),
     sample_prior = "yes")
show("C2 bounded b class lb ub", y ~ x + z, data = d,
     prior = c(set_prior("normal(0, 1)", class = "b", lb = 0, ub = 1)),
     sample_prior = "yes")
show("D mo simo", y ~ mo(m), data = d, sample_prior = "yes")
show("E beta phi", yp ~ x, data = d, family = Beta(),
     sample_prior = "yes")
show("F negbin shape", yc ~ x, data = d, family = negbinomial(),
     sample_prior = "yes")
show("G only default (flat b)", y ~ x, data = d, sample_prior = "only")
show("G2 only proper", y ~ x, data = d,
     prior = prior(normal(0, 1), class = b), sample_prior = "only")
show("H no intercept center", bf(y ~ 0 + Intercept + x), data = d,
     prior = prior(normal(0, 1), class = b), sample_prior = "yes")
show("I dpar sigma formula", bf(y ~ x, sigma ~ z), data = d,
     prior = c(prior(normal(0, 1), class = b),
               prior(normal(0, 1), class = b, dpar = sigma)),
     sample_prior = "yes")
show("J ar", y ~ x + ar(time = t, gr = g), data = d,
     prior = prior(normal(0, 1), class = b), sample_prior = "yes")
show("K student nu", y ~ x, data = d, family = student(),
     sample_prior = "yes")
show("L bernoulli", yb ~ x + (1 | g), data = d, family = bernoulli(),
     prior = prior(normal(0, 1), class = b), sample_prior = "yes")
show("M smooth", y ~ s(x), data = d, sample_prior = "yes")
show("N cumulative", yc ~ x, data = transform(d, yc = pmin(yc, 5) + 1),
     family = cumulative(), sample_prior = "yes")
show("O mvbf rescor", bf(mvbind(y, z) ~ x) + set_rescor(TRUE), data = d,
     sample_prior = "yes")
show("P nl", bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d,
     prior = c(prior(normal(1, 2), nlpar = a),
               prior(normal(0, 1), nlpar = b)),
     sample_prior = "yes")
show("Q sd uniform bounded", y ~ x + (1 | g), data = d,
     prior = c(prior(uniform(0, 5), class = sd, lb = 0, ub = 5)),
     sample_prior = "yes")
show("R mixture", y ~ x, data = d, family = mixture(gaussian, gaussian),
     prior = c(prior(normal(0, 3), Intercept, dpar = mu1),
               prior(normal(1, 3), Intercept, dpar = mu2)),
     sample_prior = "yes")
show("S hurdle_poisson hu", yc ~ x, data = d, family = hurdle_poisson(),
     sample_prior = "yes")
show("T sd for two groups", y ~ x + (1 | g) + (1 | m), data = d,
     prior = prior(normal(0, 1), class = sd, group = m),
     sample_prior = "yes")
