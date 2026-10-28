# Reviewer, short check of fix round 1.
#   Rscript dev/review-priordraws-fr1.R <brms|frm>
# brms: y | mi(sdy) under "only" (fit), me() stancode with flat and
# proper meanme/sdme. frm: the same mi(sdy) fit, the mi refusals, the
# warmup mismatch refusal, update(ds, 5).
side <- commandArgs(TRUE)[1]
.libPaths(c(if (side == "frm") "C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
LOG <- "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-log"
q <- function(expr) suppressWarnings(suppressMessages(expr))
err <- function(expr) tryCatch({expr; "OK"}, error = function(e) paste("ERROR:", conditionMessage(e)))
set.seed(1)
d <- data.frame(x = rnorm(40))
d$y <- 1 + d$x + rnorm(40)
d$sdy <- 0.4
d$xo <- d$x + rnorm(40, 0, 0.3)
d$sdx <- 0.3
dn <- d
dn$y[c(2, 9)] <- NA
if (side == "brms") {
  suppressPackageStartupMessages(library(brms))
  pr <- prior(normal(0, 1), class = b) + prior(normal(1, 2), class = Intercept) +
    prior(exponential(1), class = sigma)
  fit <- q(brm(bf(y | mi(sdy) ~ x), data = d, prior = pr, sample_prior = "only",
               chains = 4, iter = 3000, warmup = 1000, seed = 9, refresh = 0))
  dr <- as.data.frame(as_draws_df(fit))
  cat("brms mi(sdy) only variables:", head(names(dr), 8), "\n")
  set.seed(3)
  yr <- posterior_predict(fit)[, 1:3]
  saveRDS(list(draws = dr, yrep = yr), file.path(LOG, "fr1-mi-brms.rds"))
  # brms with missing y under mi(sdy): samples (flat) or refuses?
  cat("brms mi(sdy) with 2 NA, only, stancode:",
      err(stancode(bf(y | mi(sdy) ~ x), data = dn, prior = pr,
                   sample_prior = "only")), "\n")
  # me(): flat defaults and proper priors
  cat("brms me() default priors only:",
      err(stancode(bf(y ~ me(xo, sdx)), data = d, prior = pr,
                   sample_prior = "only")), "\n")
  print(subset(default_prior(bf(y ~ me(xo, sdx)), data = d),
               class %in% c("meanme", "sdme")))
  pme <- pr + prior(normal(0, 2), class = meanme) + prior(exponential(1), class = sdme)
  cat("brms me() proper meanme/sdme only:",
      err(stancode(bf(y ~ me(xo, sdx)), data = d, prior = pme,
                   sample_prior = "only")), "\n")
  sc <- stancode(bf(y ~ me(xo, sdx)), data = d, prior = pme, sample_prior = "yes")
  cat("brms me() yes, prior_ lines:\n",
      paste(" ", trimws(grep("prior_.*_rng", strsplit(sc, "\n")[[1]], value = TRUE)),
            collapse = "\n"), "\n")
} else {
  suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
  cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
  options(mc.cores = 1)
  pr <- set_prior("normal(0, 1)", class = "b") +
    set_prior("normal(1, 2)", class = "Intercept") +
    set_prior("exponential(1)", class = "sigma")
  fit <- q(frm_sample(bf(y | mi(sdy) ~ x), data = d, family = gaussian(),
                      prior = pr, sample_prior = "only", chains = 4,
                      iter = 3000, warmup = 1000, seed = 9, refresh = 0))
  dr <- as.data.frame(as_draws_df(fit))
  cat("frm mi(sdy) only variables:", head(names(dr), 8), "\n")
  set.seed(3)
  yr <- tryCatch(posterior_predict(fit)[, 1:3], error = function(e) conditionMessage(e))
  if (is.character(yr)) cat("frm posterior_predict:", yr, "\n")
  j <- frmtmb.sample:::draws_par_index(fit$fit)[["miss"]][1L]
  saveRDS(list(draws = dr, yrep = if (is.matrix(yr)) yr, yl1 = fit$draws[, j]),
          file.path(LOG, "fr1-mi-frm.rds"))
  cat("frm mi(sdy) with 2 NA, only:",
      err(q(frm_sample(bf(y | mi(sdy) ~ x), data = dn, family = gaussian(),
                       prior = pr, sample_prior = "only", chains = 1,
                       iter = 100, seed = 1, refresh = 0))), "\n")
  cat("frm plain mi() 2 NA, only:",
      err(q(frm_sample(bf(y | mi() ~ x), data = dn, family = gaussian(),
                       prior = pr, sample_prior = "only", chains = 1,
                       iter = 100, seed = 1, refresh = 0))), "\n")
  cat("frm me() only:",
      err(q(frm_sample(bf(y ~ me(xo, sdx)), data = d, family = gaussian(),
                       prior = pr, sample_prior = "only", chains = 1,
                       iter = 100, seed = 1, refresh = 0))), "\n")
  ds <- q(frm_sample(bf(y ~ x), data = d, family = gaussian(), prior = pr,
                     sample_prior = "yes", chains = 2, iter = 300, seed = 1,
                     refresh = 0))
  cat("frm me() yes names: ",
      paste(names(prior_draws(q(frm_sample(bf(y ~ me(xo, sdx)), data = d,
        family = gaussian(), prior = pr, sample_prior = "yes", chains = 1,
        iter = 200, seed = 1, refresh = 0)))), collapse = " "), "\n")
  bad1 <- ds
  bad1$draws <- bad1$draws[, -2]
  cat("mismatch (column dropped):", err(as_draws_array(bad1, inc_warmup = TRUE)), "\n")
  cat("no mismatch, inc_warmup = FALSE on that object:", err(as_draws_array(bad1)), "\n")
  cat("update(ds, 5):", err(update(ds, 5)), "\n")
  cat("update(ds, y ~ x):", err(update(ds, y ~ x)), "\n")
}
