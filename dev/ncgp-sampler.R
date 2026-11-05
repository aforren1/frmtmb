# frm_sample() on the same tape as frm(): the density at frm()'s maximum, and its
# draws on dev/ncgp-post.R's low-noise design (seed 20261105), brms's default
# priors, same chains and control as that script's brms run. Output NCGPS lines.
#   Rscript dev/ncgp-sampler.R <true sigma> <n> [iter] [delta]
args <- commandArgs(TRUE)
s <- as.numeric(args[1])
n <- as.integer(args[2])
iter <- if (length(args) > 2) as.integer(args[3]) else 1500L
delta <- if (length(args) > 3) as.numeric(args[4]) else 0.99
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20261105)
x <- round(stats::runif(n, 0, 6), 3)
d <- data.frame(x = x, y = sin(2 * x) + s * stats::rnorm(n))
fit <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian()))
ours <- -as.numeric(fit$obj$env$f(fit$obj$env$last.par.best))
nll <- frmtmb::build_objective(fit$frame)
cat(sprintf("NCGPS s=%.0e the fit's joint %.10f; build_objective at its estimates %.10f; diff %.2e\n",
            s, ours, -nll(fit$estimates), ours + nll(fit$estimates)))
plan <- frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, NULL)
cat(sprintf("NCGPS non-centered blocks to transform: %d; left centered: %d\n",
            length(plan$idx), length(plan$centered)))
t0 <- Sys.time()
fs <- suppressWarnings(frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(),
                                  chains = 4L, cores = 4L, iter = iter,
                                  warmup = iter %/% 2L, seed = 20261105,
                                  control = list(adapt_delta = delta, max_treedepth = 14L),
                                  refresh = 0))
tb <- as.numeric(Sys.time() - t0, units = "secs")
dr <- posterior::as_draws_df(fs)
vars <- grep("^(sigma|sdgp_|lscale_)", names(dr), value = TRUE)
print(as.data.frame(posterior::summarise_draws(posterior::subset_draws(dr, variable = vars),
  ~quantile(.x, c(.05, .5, .95)), "rhat", "ess_bulk")))
cat(sprintf("NCGPS s=%.0e frm_sample %.0fs\n", s, tb))
