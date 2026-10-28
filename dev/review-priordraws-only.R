# Reviewer: sample_prior = "only" on models the lane did not test, and
# posterior_predict of the result (the prior predictive).
#   Rscript dev/review-priordraws-only.R <brms|frm> <set>
# Saves parameter draws and yrep rows 1, 2 of each case to
# dev/review-priordraws-log/only-<side>-<case>.rds
# Data seed 1; sampler seed 9; 4 chains of 3000 (1000 warmup).
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
                h = factor(rep(1:5, each = 12)))
d$y <- 1 + 0.5 * d$x + rnorm(60)
d$w <- runif(60, 0.5, 2)
d$o <- rnorm(60, 0, 0.5)
d$cc <- rep(c(0, 1, 0, -1), 15)
d$t <- rep(1:12, 5)
d$xo <- d$x + rnorm(60, 0, 0.3)
d$sdx <- 0.3
d$ynl <- 2 * exp(0.3 * d$x) + rnorm(60, 0, 0.3)
base <- P("normal(0, 1)", class = "b") + P("normal(1, 2)", class = "Intercept") +
  P("exponential(1)", class = "sigma")
cases <- list(
  "1" = list(
    O1 = list(bf(y | weights(w) ~ x + offset(o) + (1 | g)),
              base + P("exponential(1)", class = "sd")),
    O2 = list(bf(y | cens(cc) ~ x), base)
  ),
  "2" = list(
    O3 = list(bf(y ~ x + ar(time = t, gr = h, cov = TRUE)),
              base + P("normal(0, 0.5)", class = "ar")),
    O4 = list(bf(ynl ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
              P("normal(1, 2)", nlpar = "a") + P("normal(0, 1)", nlpar = "b") +
                P("exponential(1)", class = "sigma"))
  ),
  "3" = list(
    O5 = list(bf(y ~ me(xo, sdx)),
              if (side == "frm") base else base + P("normal(0, 1)", class = "meanme") +
                P("exponential(1)", class = "sdme")),
    O6 = list(bf(y | trunc(lb = -3) ~ x), base)
  )
)[[set]]

for (nm in names(cases)) {
  cs <- cases[[nm]]
  cat("\n=====", nm, "=====\n")
  t0 <- Sys.time()
  fit <- tryCatch(suppressWarnings(suppressMessages(
    if (side == "brms") {
      brm(cs[[1]], data = d, family = gaussian(), prior = cs[[2]],
          sample_prior = "only", chains = 4, iter = 3000, warmup = 1000,
          seed = 9, refresh = 0)
    } else {
      frm_sample(cs[[1]], data = d, family = gaussian(), prior = cs[[2]],
                 sample_prior = "only", chains = 4, iter = 3000,
                 warmup = 1000, seed = 9, refresh = 0)
    })), error = function(e) e)
  if (inherits(fit, "error")) {
    cat("ERROR:", conditionMessage(fit), "\n")
    next
  }
  dr <- as.data.frame(posterior::as_draws_df(fit))
  cat("variables:", paste(head(names(dr), 30), collapse = " "), "\n")
  set.seed(3)
  yr <- tryCatch(suppressWarnings(posterior_predict(fit)),
                 error = function(e) e)
  if (inherits(yr, "error")) {
    cat("PP ERROR:", conditionMessage(yr), "\n")
    yr <- NULL
  } else {
    yr <- yr[, 1:3, drop = FALSE]
  }
  rh <- tryCatch(max(posterior::summarise_draws(
    posterior::as_draws_array(fit), "rhat")$rhat, na.rm = TRUE),
    error = function(e) NA)
  cat("max rhat", rh, " elapsed", format(Sys.time() - t0), "\n")
  saveRDS(list(draws = dr, yrep = yr, rhat = rh),
          file.path(LOG, sprintf("only-%s-%s.rds", side, nm)))
}
