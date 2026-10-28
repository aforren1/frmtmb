# Reviewer: sample_prior = "yes" on models the lane did not test.
#   Rscript dev/review-priordraws-yes.R <brms|frm> <set>
# Saves the prior draws of each case to dev/review-priordraws-log/yes-<side>-<case>.rds
# Data seed 1; sampler seed 5; 2 chains of 3000 (1000 warmup) = 4000 draws.
side <- commandArgs(TRUE)[1]
set <- commandArgs(TRUE)[2]
.libPaths(c(if (side == "frm") "C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-stancache")
LOG <- "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-log"
if (side == "brms") {
  suppressPackageStartupMessages(library(brms))
  P <- function(...) brms::set_prior(...)
} else {
  suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
  P <- function(...) frmtmb::set_prior(...)
  cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
}
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)),
                h = factor(rep(1:5, each = 12)), m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(60)
d$yc <- rpois(60, 3)
d$yo <- pmin(d$yc, 4) + 1
d$t <- rep(1:12, 5)
d$ynl <- 2 * exp(0.3 * d$x) + rnorm(60, 0, 0.3)

cases <- list(
  "1" = list(
    T1 = list(bf(y ~ x + (1 | g)) + bf(z ~ x) + set_rescor(FALSE), gaussian(),
              P("normal(0, 1)", class = "b", resp = "y") +
                P("student_t(3, 0, 1)", class = "sd", resp = "y") +
                P("exponential(1)", class = "sigma", resp = "z") +
                P("normal(0, 2)", class = "b", coef = "x", resp = "z")),
    T2 = list(bf(y ~ x + (1 | g) + (1 | h), sigma ~ x), gaussian(),
              P("normal(0, 1)", class = "Intercept", dpar = "sigma") +
                P("student_t(3, 0, 1)", class = "sd", group = "h") +
                P("normal(0, 0.5)", class = "b", dpar = "sigma")),
    T3 = list(bf(y ~ x + (1 + x | g) + (1 + x + z | h)), gaussian(),
              P("lkj(1.5)", class = "cor") +
                P("lkj(3)", class = "cor", group = "h") +
                P("gamma(2, 2)", class = "sd"))
  ),
  "2" = list(
    T4 = list(bf(y ~ x + z), gaussian(),
              P("student_t(3, 0, 1)", class = "b", lb = 0) +
                P("normal(1, 1)", class = "b", coef = "z") +
                P("inv_gamma(3, 2)", class = "sigma")),
    T5 = list(bf(ynl ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
              gaussian(),
              P("normal(1, 2)", nlpar = "a") +
                P("normal(0, 1)", nlpar = "b", lb = 0) +
                P("exponential(1)", class = "sd", nlpar = "a")),
    T6 = list(bf(yc ~ x + (1 | g)), negbinomial(),
              P("gamma(2, 0.5)", class = "shape") +
                P("logistic(0, 1)", class = "b"))
  ),
  "3" = list(
    T7 = list(bf(yc ~ x), zero_inflated_poisson(),
              P("beta(2, 2)", class = "zi") +
                P("normal(1, 1)", class = "Intercept")),
    T8 = list(bf(yo ~ x), cumulative(),
              P("normal(0, 2)", class = "Intercept", coef = "1") +
                P("normal(0, 1)", class = "b")),
    T9 = list(bf(y ~ x + ma(time = t, gr = h)), gaussian(),
              P("normal(0, 0.5)", class = "ma")),
    T10 = list(bf(y ~ x + ar(time = t, gr = h, cov = TRUE)), gaussian(),
               P("normal(0, 0.5)", class = "ar")),
    T11 = list(bf(y ~ x + cosy(time = t, gr = h)), gaussian(),
               P("beta(2, 2)", class = "cosy"))
  )
)[[set]]

for (nm in names(cases)) {
  cs <- cases[[nm]]
  cat("\n=====", nm, "=====\n")
  t0 <- Sys.time()
  fit <- tryCatch(suppressWarnings(suppressMessages(
    if (side == "brms") {
      brm(cs[[1]], data = d, family = cs[[2]], prior = cs[[3]],
          sample_prior = "yes", chains = 2, iter = 3000, warmup = 1000,
          seed = 5, refresh = 0)
    } else {
      frm_sample(cs[[1]], data = d, family = cs[[2]], prior = cs[[3]],
                 sample_prior = "yes", chains = 2, iter = 3000,
                 warmup = 1000, seed = 5, refresh = 0)
    })), error = function(e) e)
  if (inherits(fit, "error")) {
    cat("ERROR:", conditionMessage(fit), "\n")
    next
  }
  pd <- prior_draws(fit)
  cat("names:", paste(names(pd), collapse = " "), "\n")
  cat("dim:", dim(pd), " elapsed", format(Sys.time() - t0), "\n")
  saveRDS(pd, file.path(LOG, sprintf("yes-%s-%s.rds", side, nm)))
}
