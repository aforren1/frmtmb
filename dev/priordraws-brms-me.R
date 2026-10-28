# brms 2.23.0 on me() priors under sample_prior = "yes" and "only", and
# on coefficient-specific ordinal threshold priors (fix round 1, items
# 5 and 6). stancode() and default_prior(), then one small fit each for
# the prior_ names.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
set.seed(1)
d <- data.frame(xo = rnorm(60), sdx = 0.3, z = rnorm(60))
d$x <- d$xo + rnorm(60, 0, 0.3)
d$y <- 1 + 0.5 * d$xo + rnorm(60)
d$yo <- cut(d$y, c(-Inf, 0, 1, 2, Inf), labels = FALSE)
show_block <- function(sc, pat) {
  sc <- strsplit(sc, "\n")[[1]]
  cat(grep(pat, sc, value = TRUE), sep = "\n")
}
q <- function(lab, expr) {
  cat("\n=====", lab, "=====\n")
  tryCatch(expr, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}
q("me default_prior", print(default_prior(y ~ me(x, sdx), data = d)))
q("me yes stancode", show_block(stancode(
  y ~ me(x, sdx), data = d, prior = prior(normal(0, 1), class = b),
  sample_prior = "yes"), "prior_|lprior|meanme|sdme|Xme"))
q("me only default", stancode(y ~ me(x, sdx), data = d,
                              prior = prior(normal(0, 1), class = b),
                              sample_prior = "only"))
pr_me <- c(prior(normal(0, 1), class = b),
           prior(normal(0, 2), class = meanme),
           prior(exponential(1), class = sdme))
q("me only with meanme prior", show_block(stancode(
  y ~ me(x, sdx), data = d, prior = pr_me, sample_prior = "only"),
  "prior_only|target|lprior"))
q("me yes fit", {
  f <- suppressWarnings(suppressMessages(brm(
    y ~ me(x, sdx), data = d, prior = pr_me, sample_prior = "yes",
    chains = 1, iter = 200, seed = 1, refresh = 0)))
  cat("variables:", variables(f), "\n")
  cat("prior_draws:", names(prior_draws(f)), "\n")
})
q("ordinal coef stancode", show_block(stancode(
  yo ~ z, data = d, family = cumulative(),
  prior = c(prior(normal(0, 2), class = Intercept, coef = 1),
            prior(normal(0, 1), class = b)),
  sample_prior = "yes"), "prior_|ordered|Intercept"))
q("ordinal coef fit", {
  f <- suppressWarnings(suppressMessages(brm(
    yo ~ z, data = d, family = cumulative(),
    prior = c(prior(normal(0, 2), class = Intercept, coef = 1),
              prior(normal(0, 1), class = b)),
    sample_prior = "yes", chains = 1, iter = 200, seed = 1, refresh = 0)))
  cat("prior_draws:", names(prior_draws(f)), "\n")
  print(round(sapply(prior_draws(f), sd), 2))
})
q("ordinal all coef, flat class", show_block(stancode(
  yo ~ z, data = d, family = cumulative(),
  prior = c(prior(normal(0, 2), class = Intercept, coef = 1),
            set_prior("", class = "Intercept")),
  sample_prior = "yes"), "prior_"))
