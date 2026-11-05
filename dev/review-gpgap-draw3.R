# Short check of fix round 1, item 1: why does the drawn field at true sigma
# 1e-3 mix worse than brms's? Same data as dev/gpgap-sampler.R (seed 6601,
# n 15), same settings on every arm, Stan seed as given.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-draw.R <log> <route> <stan seed> <iter> <td>
#     frm        frm_sample(), formula route, default init ("random")
#     frm_warm   frm_sample(), every chain from the target's typical set: brms's
#                sigma, sdgp and lscale, field through the fit's factor
#     frm_fit    frm_sample(fit), the mode start at frm()'s collapsed sigma
#     brms       brms, rows as drawn
#     brms_sort  brms, rows sorted on x, so its factor is over sorted positions
a <- commandArgs(trailingOnly = TRUE)
route <- a[1]
sseed <- as.integer(a[2])
iter <- as.integer(a[3])
td <- as.integer(a[4])
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(6601)
n <- 15
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
d <- data.frame(x = x, y = sin(2 * x) + as.numeric(Sys.getenv("REV_S")) * e0)
if (route == "brms_sort") d <- d[order(d$x), ]
chains <- 4L
ctl <- list(adapt_delta = 0.95, max_treedepth = td)
warned <- character()
keep <- function(w) {
  warned <<- c(warned, gsub("[\r\n]+", " ", substr(conditionMessage(w), 1, 90)))
  invokeRestart("muffleWarning")
}
t0 <- proc.time()[["elapsed"]]
arg <- list(chains = chains, iter = iter, warmup = iter %/% 2, seed = sseed,
            refresh = 0, cores = chains, control = ctl)
fit <- withCallingHandlers(switch(route,
  brms = , brms_sort = do.call(brms::brm, c(list(brms::bf(y ~ gp(x)), data = d,
                                                 family = gaussian(), silent = 2), arg)),
  frm = do.call(frm_sample, c(list(bf(y ~ gp(x)), data = d, family = gaussian()), arg)),
  frm_fit = {
    f0 <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
    cat(sprintf("fit: sigma %.3g\n", as.numeric(sigma(f0))))
    do.call(frm_sample, c(list(f0), arg))
  },
  frm_warm = {
    # from the marginal's target (gpgap-findings table): sigma 1.03e-3, sdgp 2.17,
    # lscale 0.2539 on brms's scale
    f0 <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
    dmax <- diff(range(d$x))
    u <- sort(unique(d$x))
    sdg <- 2.17; rho <- 0.2539 * dmax; sg <- 1.03e-3
    K <- sdg^2 * exp(-outer(u, u, "-")^2 / (2 * rho^2)) + diag(1e-12 * sdg^2, length(u))
    J <- match(d$x, u)
    V <- K[J, J] + diag(sg^2, n)
    fm <- as.numeric(K[, J] %*% solve(V, d$y - mean(d$y)))
    z <- forwardsolve(t(chol(K)), fm)
    lp <- f0$obj$env$last.par.best
    nm <- names(lp)
    v <- as.numeric(lp)
    v[nm == "beta"] <- mean(d$y)
    v[nm == "betad"] <- log(sg)
    v[nm == "theta"] <- c(log(sdg), log(rho))
    v[nm == "b"] <- z
    set.seed(sseed)
    inits <- lapply(seq_len(chains), function(i) v + rnorm(length(v), 0, 0.05))
    do.call(frm_sample, c(list(f0, init = inits), arg))
  }), warning = keep)
wall <- proc.time()[["elapsed"]] - t0
dd <- posterior::as_draws_df(fit)
vars <- posterior::variables(dd)
sf <- if (startsWith(route, "brms")) fit$fit else fit$stanfit
sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
for (pat in c("^sigma$", "^sdgp_", "^lscale_")) {
  for (nm in grep(pat, vars, value = TRUE)) {
    m <- posterior::extract_variable_matrix(dd, nm)
    q <- stats::quantile(as.numeric(m), c(0.05, 0.5, 0.95))
    cat(sprintf("DRAW %s seed %d %-11s q05 %.4g q50 %.4g q95 %.4g rhat %.3f ess %.0f\n",
                route, sseed, nm, q[1], q[2], q[3], posterior::rhat(m),
                posterior::ess_bulk(m)))
  }
}
ss <- vapply(sp, function(m) mean(m[, "stepsize__"]), 0)
cat(sprintf("DRAW %s seed %d iter %d td %d divergent %d at-depth %d of %d stepsize %s wall %.0f\n",
            route, sseed, iter, td,
            sum(vapply(sp, function(m) sum(m[, "divergent__"]), 0)),
            sum(vapply(sp, function(m) sum(m[, "treedepth__"] >= td), 0)),
            sum(vapply(sp, nrow, 0L)), paste(signif(ss, 2), collapse = ","), wall))
cat("TOLD:", paste(unique(warned), collapse = " | "), "\n")
