# frm_sample() on an exact gp() field, centered (rellib-r7) against
# non-centered (the lane), and brms 2.23.0 on the same data. Fixtures:
# lane surface's y ~ gp(x) (60 points, data seed 5,
# dev/surface-gp-sweep.R) and a gp(x, by = f) on the same x. 4 chains
# of 1000 (500 warmup), seed 1. Per quantity: mean, sd, R-hat, bulk
# ESS; per run: divergences, transitions at the maximum tree depth,
# and the sampler's wall time.
#   Rscript dev/moparity-gp-ncp.R base|lane|brms
args <- commandArgs(TRUE)
arm <- args[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
set.seed(5)
n <- 60
d9 <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d9$y <- 0.5 + sin(d9$x) + stats::rnorm(n, 0, 0.3)
d9$f <- factor(rep(c("a", "b"), length.out = n))
d9$y2 <- d9$y + ifelse(d9$f == "b", 0.5 * cos(d9$x), 0)
cases <- list(gp = y ~ gp(x), gpby = y2 ~ gp(x, by = f))
summ <- function(M, nch, cols) {
  arr <- array(M[, cols], c(nrow(M) / nch, nch, length(cols)),
               dimnames = list(NULL, NULL, cols))
  s <- posterior::summarise_draws(posterior::as_draws_array(arr),
                                  "mean", "sd", "rhat", "ess_bulk",
                                  "mcse_mean", "mcse_sd")
  as.data.frame(s)
}
res <- list()
for (nm in names(cases)) {
  fo <- cases[[nm]]
  routes <- if (arm == "brms") "brms" else c("formula", "fit")
  for (rt in routes) {
    t0 <- proc.time()[["elapsed"]]
    if (rt == "brms") {
      b <- brms::brm(fo, data = d9, chains = 4, iter = 1000, seed = 1,
                     cores = 1, refresh = 0)
      t1 <- proc.time()[["elapsed"]]
      M <- as.matrix(b)
      sp <- rstan::get_sampler_params(b$fit, inc_warmup = FALSE)
      td <- 10
    } else {
      s <- if (rt == "formula") {
        suppressWarnings(suppressMessages(frm_sample(
          bf(fo), data = d9, family = gaussian(), chains = 4, iter = 1000,
          seed = 1, cores = 1, refresh = 0)))
      } else {
        fit <- frm(bf(fo), family = gaussian(), data = d9)
        suppressWarnings(suppressMessages(frm_sample(
          fit, chains = 4, iter = 1000, seed = 1, cores = 1, refresh = 0)))
      }
      t1 <- proc.time()[["elapsed"]]
      M <- s$draws
      sp <- rstan::get_sampler_params(s$stanfit, inc_warmup = FALSE)
      td <- s$stanfit@stan_args[[1]]$control$max_treedepth %||% 10
    }
    cols <- grep("^(b_Intercept|sigma|sdgp_|lscale_)", colnames(M),
                 value = TRUE)
    S <- summ(M, 4, cols)
    div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
    tdh <- sum(sapply(sp, function(x) sum(x[, "treedepth__"] >= td)))
    cat(sprintf("\n== %s %s %s: divergences %d, at max tree depth %d, sampler seconds %.1f\n",
                nm, arm, rt, div, tdh, t1 - t0))
    print(S, digits = 4, row.names = FALSE)
    res[[paste(nm, rt)]] <- list(S = S, div = div, tdh = tdh,
                                 secs = t1 - t0)
  }
}
saveRDS(res, file.path("dev/moparity-log", paste0("gp-ncp-", arm, ".rds")))
