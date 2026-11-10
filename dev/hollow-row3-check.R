# Lane hollow, Row 3: on the base build, frmtmb's log_lik() and
# posterior_predict() of a custom family come from the family's own
# lpdf and sim slots, not from a brms-style function found by name.
#   Rscript dev/hollow-row3-check.R > dev/hollow-out/row3-check.txt 2>&1
# Library: rellib-r19 then the user library (the base build).
.libPaths(c("C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-hollow/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
cat("lib:", find.package("frmtmb.sample"),
    format(packageVersion("frmtmb.sample")), "\n")

# A marker in each slot: the log density carries +0.125 and the
# simulator returns mu + 1000 exactly, so a value that came from either
# slot is recognizable, and neither changes the posterior of mu.
marker_ll <- 0.125
calls <- new.env()
calls$sim <- 0L
mycust <- custom_family(
  "mycust", dpars = "mu", links = list(mu = "identity"),
  lpdf = function(y, dpars, aterms) {
    -0.5 * (y - dpars$mu)^2 - 0.5 * log(2 * pi) + 0.125
  },
  init_dpars = list(mu = function(y, aterms) mean(y)),
  post = list(mean_fn = function(dpars, aterms) dpars$mu),
  sim = function(dpars, aterms, n) {
    calls$sim <- calls$sim + 1L
    dpars$mu + 1000
  })
# brms's dispatch looks these names up; frmtmb must not reach them
log_lik_mycust <- function(i, prep) stop("brms-style log_lik reached")
posterior_predict_mycust <- function(i, prep) {
  stop("brms-style posterior_predict reached")
}

set.seed(603)
dd <- data.frame(x = rnorm(30))
dd$y <- rnorm(30, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x), family = mycust, data = dd)
ds <- frm_sample(fit, chains = 1, iter = 150, warmup = 50, seed = 440,
                 refresh = 0)
ll <- log_lik(ds)
mu <- posterior_epred(ds)
pp <- posterior_predict(ds)
cat("dim log_lik:", dim(ll), " dim epred:", dim(mu), " dim predict:",
    dim(pp), " ndraws:", ndraws(ds), " nobs:", nrow(dd), "\n")
ref_ll <- -0.5 * sweep(mu, 2, dd$y)^2 - 0.5 * log(2 * pi) + marker_ll
cat("log_lik vs the family's lpdf at each draw's mu, max abs diff:",
    format(max(abs(ll - ref_ll)), digits = 3), "\n")
cat("log_lik vs the same density WITHOUT the marker, min abs diff:",
    format(min(abs(ll - (ref_ll - marker_ll))), digits = 6), "\n")
cat("posterior_predict - epred, range:",
    format(range(pp - mu), digits = 10), "\n")
cat("sim slot calls:", calls$sim, "\n")
cat("ROW3 CHECK END\n")
