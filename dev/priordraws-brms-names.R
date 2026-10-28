# brms 2.23.0's prior draw columns after fitting (sample_prior = "yes"):
# their names, their position among variables(), and what
# prior_draws() returns. One case set per process:
#   Rscript dev/priordraws-brms-names.R <set>
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
set <- commandArgs(TRUE)[1]
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)),
                h = factor(rep(1:5, each = 12)), m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(60)
d$yb <- rbinom(60, 1, 0.5)
d$yp <- rbeta(60, 2, 2)
d$yc <- rpois(60, 3)
d$yo <- pmin(d$yc, 4) + 1
d$t <- rep(1:10, 6)
cases <- list(
  "1" = list(
    A = list(y ~ x + (1 | g)),
    B = list(y ~ x + z + (1 + x | g),
             prior = c(prior(normal(0, 1), class = b),
                       prior(normal(0, 2), class = b, coef = z),
                       prior(normal(1, 3), class = Intercept),
                       prior(exponential(2), class = sigma),
                       prior(lkj(2), class = cor),
                       prior(cauchy(0, 1), class = sd, group = g,
                             coef = x))),
    B2 = list(y ~ x + (1 + x | g) + (1 | h),
              prior = c(prior(normal(0, 1), class = sd))),
    B3 = list(y ~ x + z,
              prior = c(prior(normal(0, 2), class = b, coef = z))),
    C2 = list(y ~ x + z,
              prior = set_prior("normal(0, 1)", class = "b", lb = 0,
                                ub = 1)),
    H = list(bf(y ~ 0 + Intercept + x),
             prior = prior(normal(0, 1), class = b))
  ),
  "2" = list(
    D = list(y ~ mo(m)),
    D2 = list(y ~ mo(m) + x, prior = prior(normal(0, 1), class = b)),
    I = list(bf(y ~ x, sigma ~ z),
             prior = c(prior(normal(0, 1), class = b),
                       prior(normal(0, 1), class = b, dpar = sigma))),
    O = list(bf(mvbind(y, z) ~ x) + set_rescor(TRUE),
             prior = c(prior(normal(0, 1), class = b, resp = y))),
    R = list(y ~ x, family = mixture(gaussian, gaussian),
             prior = c(prior(normal(0, 3), Intercept, dpar = mu1),
                       prior(normal(1, 3), Intercept, dpar = mu2))),
    P = list(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
             prior = c(prior(normal(1, 2), nlpar = a),
                       prior(normal(0, 1), nlpar = b)))
  ),
  "3" = list(
    N = list(yo ~ x, family = cumulative(),
             prior = prior(normal(0, 1), class = b)),
    M = list(y ~ s(x)),
    E = list(yp ~ x, family = Beta()),
    K = list(y ~ x, family = student()),
    S = list(yc ~ x, family = hurdle_poisson()),
    AR = list(y ~ x + ar(time = t, gr = h),
              prior = prior(normal(0, 0.5), class = ar)),
    GP = list(y ~ gp(x), prior = prior(normal(0, 1), class = lscale))
  )
)[[set]]
for (nm in names(cases)) {
  cs <- cases[[nm]]
  args <- c(list(formula = cs[[1]], data = d, sample_prior = "yes",
                 chains = 1, iter = 200, seed = 1, refresh = 0),
            cs[-1])
  fit <- tryCatch(suppressWarnings(suppressMessages(do.call(brm, args))),
                  error = function(e) e)
  cat("\n=====", nm, "=====\n")
  if (inherits(fit, "error")) {
    cat("ERROR:", conditionMessage(fit), "\n")
    next
  }
  v <- variables(fit)
  cat("variables:", paste(v, collapse = " "), "\n")
  pd <- prior_draws(fit)
  cat("prior_draws names:", paste(names(pd), collapse = " "), "\n")
  cat("prior_draws dim:", dim(pd), "\n")
  print(utils::head(pd, 2))
  cat("ndraws", ndraws(fit), "nrow as_draws_df inc_warmup",
      nrow(as_draws_df(fit, inc_warmup = TRUE)), "\n")
  print(fit$prior)
}
