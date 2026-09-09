# lane tmbstan, probe 08: the false-alarm rate and the cost of the
# grad_log_prob assertion, across model shapes rather than one.
#
# probe 07 measured the separation on one gaussian model. What decides
# whether this becomes a test is (a) does the clean gap stay at zero on
# other shapes, since a check that fires on a correct model is worse
# than no check, and (b) what does it cost.
#
# Only FIXED-EFFECTS models. With a random block, fit$obj marginalizes
# by Laplace while tmbstan samples the random effects too, so obj$par
# and the unconstrained vector are different objects and the identity
# does not apply. The defect is in the generated model, not in any one
# model's tape, so a fixed-effects model tests it fully.
#
# SEEDS: 4021, 4023, 4024, 4025 for the four data sets; sampler seed 11.
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-tmbstan/dev/stan-cache")
library(frmtmb)
library(frmtmb.sample)

cases <- list()
cases$gaussian <- local({
  set.seed(4021)
  dd <- data.frame(x = stats::rnorm(80))
  dd$y <- stats::rnorm(80, 1 + 0.5 * dd$x, 1)
  frm(bf(y ~ x) + gaussian(), data = dd)
})
cases$poisson <- local({
  set.seed(4023)
  dd <- data.frame(x = stats::rnorm(80))
  dd$y <- stats::rpois(80, exp(0.4 + 0.3 * dd$x))
  frm(bf(y ~ x) + poisson(), data = dd)
})
cases$bernoulli <- local({
  set.seed(4024)
  dd <- data.frame(x = stats::rnorm(120))
  dd$y <- stats::rbinom(120, 1, stats::plogis(0.2 + 0.8 * dd$x))
  frm(bf(y ~ x) + bernoulli(), data = dd)
})
cases$sigma_pred <- local({
  set.seed(4025)
  dd <- data.frame(x = stats::rnorm(120))
  dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x, exp(-0.2 + 0.3 * dd$x))
  frm(bf(y ~ x, sigma ~ x) + gaussian(), data = dd)
})

out <- list()
for (nm in names(cases)) {
  fit <- cases[[nm]]
  t0 <- proc.time()[["elapsed"]]
  sf <- suppressWarnings(suppressMessages(
    as_tmbstan(fit, chains = 1L, iter = 20L, warmup = 10L, refresh = 0,
               seed = 11)))
  t_fit <- proc.time()[["elapsed"]] - t0
  np <- rstan::get_num_upars(sf)
  set.seed(4030)
  pts <- list(rep(0, np), fit$obj$par + 0.3, fit$obj$par - 0.7,
              stats::rnorm(np), fit$obj$par * 3)
  worst_clean <- 0
  best_broken <- Inf
  ulp <- 0
  bitwise <- TRUE
  t1 <- proc.time()[["elapsed"]]
  for (p in pts) {
    u <- as.numeric(p)
    g_stan <- as.numeric(rstan::grad_log_prob(sf, u))
    g_obj <- -as.numeric(fit$obj$gr(u))
    g_bad <- -u                     # the defect's closed form
    s <- max(abs(g_obj))
    worst_clean <- max(worst_clean, max(abs(g_stan - g_obj)) / s)
    best_broken <- min(best_broken, max(abs(g_bad - g_obj)) / s)
    ulp <- max(ulp, max(abs(g_stan - g_obj)) /
                 max(.Machine$double.eps * abs(g_obj)))
    bitwise <- bitwise && identical(g_stan, g_obj)
  }
  t_assert <- proc.time()[["elapsed"]] - t1
  out[[nm]] <- data.frame(
    model = nm, npars = np, points = length(pts),
    clean_rel = worst_clean, broken_rel = best_broken,
    max_ulp = ulp, bitwise = bitwise,
    fit_s = round(t_fit, 2), assert_s = round(t_assert, 3),
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, out)
print(res, row.names = FALSE, digits = 6)

cat("\n== verdict ==\n")
cat("worst clean relative gap over all models and points:",
    format(max(res$clean_rel), digits = 4), "\n")
cat("best broken relative gap over all models and points:",
    format(min(res$broken_rel), digits = 4), "\n")
cat("bitwise identical everywhere:", all(res$bitwise), "\n")
cat("total as_tmbstan seconds:", sum(res$fit_s), "\n")
cat("total assertion seconds :", sum(res$assert_s), "\n")
