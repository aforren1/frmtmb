# frm_sample() on the latent-residual model against brms 2.23.0, same data,
# brms's default priors on both sides: does the route sample the same model?
# Seed 167 for the data, 11 for the samplers. Usage: Rscript dev/latres-sample.R
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
  library(rstan)
})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- function(a, b) if (is.null(a)) b else a
sys.source("dev/equiv-lib.R", envir = hel)

set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:30))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 5, sd = 0.7)))
eta <- 0.8 + 0.4 * d$x + e
d$count <- rpois(n, exp(eta))

iter <- 1500L
args <- commandArgs(TRUE)
which <- if (length(args)) strsplit(args[[1]], ",")[[1]] else
  c("ar(t, g, cov = TRUE)", "ar(t, g)")
for (tm in which) {
  cat("\n=====", tm, "\n")
  f <- as.formula(paste("count ~ x +", tm))
  fb <- brms::bf(f)
  cat("--- default_prior, brms then frm_sample route\n")
  print(as.data.frame(brms::default_prior(fb, data = d, family = brms::brmsfamily("poisson")))[,
        c("prior", "class", "coef", "lb", "ub")])
  fit <- frm(frmtmb::bf(f), family = poisson(), data = d)
  print(default_prior(frmtmb::bf(f), family = poisson(), data = d,
                      route = "sample"))
  t0 <- proc.time()[[3]]
  dr <- frm_sample(frmtmb::bf(f), family = poisson(), data = d, chains = 2,
                   iter = iter, warmup = 700, seed = 11, refresh = 0)
  cat("frm_sample seconds", proc.time()[[3]] - t0, "\n")
  print(prior_summary(dr))
  print(utils::head(variables(dr), 12))
  sm <- as.data.frame(posterior::summarise_draws(as_draws_df(dr)))
  print(sm[sm$variable %in% c("b_Intercept", "b_x", "ar[1]", "sderr"),
           c("variable", "mean", "sd", "rhat", "ess_bulk")])
  code <- brms::stancode(fb, data = d, family = brms::brmsfamily("poisson"))
  sdat <- brms::standata(fb, data = d, family = brms::brmsfamily("poisson"))
  mod <- hel$brms_stan_model(code)
  sf <- suppressWarnings(rstan::sampling(mod, data = sdat, chains = 2,
                                         iter = iter, warmup = 700, seed = 11,
                                         refresh = 0, pars = "zerr",
                                         include = FALSE))
  sb <- rstan::summary(sf, pars = c("b_Intercept", "b", "ar", "sderr"))$summary
  print(round(sb[, c("mean", "sd", "Rhat", "n_eff")], 3))
  ep <- tryCatch(posterior_epred(dr), error = function(e) conditionMessage(e))
  if (is.matrix(ep)) {
    cat("epred dim", dim(ep), " mean row1", mean(ep[, 1]), "\n")
    print(summary(colMeans(ep)))
  } else print(ep)
  er <- rstan::extract(sf, pars = "err")$err
  mu_b <- colMeans(exp(
    as.vector(rstan::extract(sf, "b_Intercept")$b_Intercept) +
      outer(rstan::extract(sf, "b")$b[, 1], d$x) + er))
  if (is.matrix(ep) && grepl("cov = TRUE", tm)) {
    cat("posterior mean of the expected counts, frm_sample vs brms: mean ",
        "abs diff", mean(abs(colMeans(ep) - mu_b)), "\n")
  }
  for (mth in c("log_lik", "posterior_predict")) {
    v <- tryCatch(do.call(mth, list(dr)), error = function(e) conditionMessage(e))
    cat(mth, ":", if (is.character(v)) substr(v, 1, 300) else dim(v), "\n")
  }
  nd <- d[d$g == "1", ]
  v <- tryCatch(posterior_epred(dr, newdata = nd),
                error = function(e) conditionMessage(e))
  cat("epred newdata:", if (is.character(v)) substr(v, 1, 300) else dim(v), "\n")
  v <- tryCatch(posterior_predict(dr, newdata = nd),
                error = function(e) conditionMessage(e))
  cat("predict newdata:", if (is.character(v)) substr(v, 1, 300) else dim(v), "\n")
}
