# frm_sample() against brms at low noise, and on y in other units.
#
#   Rscript dev/gpgap-sampler.R <arm> <route> <n> <iter> <chains> <delta> <seed> 
#     [<max_treedepth> [<metric>]]
#     arm    low: sin(2 x) + 1e-3 e;  units: 1000 sin(2 x) + 100 e (the
#            review's M7, where brms's absolute line breaks)
#     route  frm_sample, laplace (frm_sample(laplace = TRUE)), base (rellib-r16's
#            frm_sample), baselap (rellib-r16's laplace), or brms
#
# Prints sigma's quantiles, Rhat and ESS, sdgp's and lscale's, divergences, wall
# seconds. No thinning; a run is one process.
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
route <- a[2]
n <- as.integer(a[3])
iter <- as.integer(a[4])
chains <- as.integer(a[5])
delta <- as.numeric(a[6])
seed <- as.integer(a[7])
td <- as.integer(if (length(a) >= 8) a[8] else 12L)
metric <- if (length(a) >= 9) a[9] else "diag_e"
libs <- if (route %in% c("base", "baselap")) {
  "C:/Users/adf44/source/r/rellib-r16"
} else {
  c("C:/Users/adf44/source/r/wt-gpgap-lib",
    "C:/Users/adf44/source/r/rellib-r16")
}
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(find.package("frmtmb")), "sample:",
    dirname(find.package("frmtmb.sample")), "\n")
set.seed(seed)
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
y <- if (arm == "low") sin(2 * x) + 1e-3 * e0 else 1000 * sin(2 * x) + 100 * e0
d <- data.frame(x = x, y = y)
warned <- character()
t0 <- proc.time()[["elapsed"]]
ctl <- list(adapt_delta = delta, max_treedepth = td, metric = metric)
keep <- function(w) {
  warned <<- c(warned, paste0(class(w)[1], ": ",
                              substr(conditionMessage(w), 1, 100)))
  invokeRestart("muffleWarning")
}
fit <- withCallingHandlers(if (route == "brms") {
  brms::brm(brms::bf(y ~ gp(x)), data = d, family = gaussian(), chains = chains,
            iter = iter, warmup = iter %/% 2, seed = seed, refresh = 0,
            cores = chains, silent = 2, control = ctl)
} else if (route %in% c("fit", "fitlap")) {
  # frm_sample(fit), where frm() collapses at low noise
  f0 <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian()))
  cat(sprintf("frm() sigma %.3g\n", as.numeric(sigma(f0))))
  frm_sample(f0, chains = chains, iter = iter, warmup = iter %/% 2, seed = seed,
             refresh = 0, cores = chains, control = ctl,
             laplace = route == "fitlap")
} else {
  frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(), chains = chains,
             iter = iter, warmup = iter %/% 2, seed = seed, refresh = 0,
             cores = chains, control = ctl,
             laplace = route %in% c("laplace", "baselap"))
}, warning = keep)
wall <- proc.time()[["elapsed"]] - t0
dd <- posterior::as_draws_df(fit)
vars <- posterior::variables(dd)
sf <- if (route == "brms") fit$fit else fit$stanfit
sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
divs <- sum(vapply(sp, function(m) sum(m[, "divergent__"]), 0))
tree <- sum(vapply(sp, function(m) sum(m[, "treedepth__"] >= td), 0))
for (pat in c("^sigma$", "^sdgp_", "^lscale_")) {
  for (nm in grep(pat, vars, value = TRUE)) {
    m <- posterior::extract_variable_matrix(dd, nm)
    q <- stats::quantile(as.numeric(m), c(0.05, 0.5, 0.95))
    cat(sprintf(paste0("GPSAMP %s %s n %d seed %d %-10s q05 %.4g q50 %.4g q95 ",
                       "%.4g rhat %.3f ess %.0f\n"),
                arm, route, n, seed, nm, q[1], q[2], q[3],
                posterior::rhat(m), posterior::ess_bulk(m)))
  }
}
cat(sprintf("GPSAMP %s %s n %d seed %d iter %d chains %d delta %.3f td %d metric %s divergent %d treedepth %d wall %.0f\n",
            arm, route, n, seed, iter, chains, delta, td, metric, divs, tree,
            wall))
cat("warnings:", paste(unique(warned), collapse = " | "), "\n")
