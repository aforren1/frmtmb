# Probe 6 (review finding 7): does frmtmb's draws-level
# posterior_predict() of a wiener fit keep the draws x observations
# matrix for ONE row, or drop it to a vector? And what does the
# rowwise simulator the adapter calls return? Log:
# dev/bin2port-out/probe6.txt. Seed 420.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.sample)
})
cat("frmtmb", format(packageVersion("frmtmb")), "frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "frmtmb.sample",
    format(packageVersion("frmtmb.sample")), "\n")
set.seed(420)
dat <- ddm_simulate(200, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
fit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5), family = wiener(),
           data = dat)
ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 300,
                                  warmup = 150, refresh = 0, seed = 420))
p_all <- posterior_predict(ds)
cat("posterior_predict(draws), all rows: dim", dim(p_all), "\n")
try_dim <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    cat(label, ": ERROR", conditionMessage(r), "\n")
  } else {
    cat(label, ": class", paste(class(r), collapse = "/"), " dim",
        if (is.null(dim(r))) "NULL" else dim(r), "\n")
  }
}
cat("newdata columns:", names(dat), "\n")
try_dim("posterior_predict(draws, newdata = 1 row)",
        posterior_predict(ds, newdata = dat[1, , drop = FALSE]))
try_dim("posterior_predict(draws, newdata = 2 rows)",
        posterior_predict(ds, newdata = dat[1:2, , drop = FALSE]))
try_dim("posterior_predict(draws, newdata = all rows)",
        posterior_predict(ds, newdata = dat))
try_dim("predict(fit, newdata = 1 row, summary = FALSE)",
        predict(fit, newdata = dat[1, , drop = FALSE], summary = FALSE))
try_dim("predict(fit, summary = FALSE), all rows",
        predict(fit, summary = FALSE))
try_dim("simulate(fit, nsim = 3)", simulate(fit, nsim = 3))
# a family without a boundary term, the same newdata route
dg <- data.frame(y = rnorm(30), x = rnorm(30))
fg <- frm(bf(y ~ x), family = gaussian(), data = dg)
dsg <- suppressWarnings(frm_sample(fg, chains = 1, iter = 300,
                                   warmup = 150, refresh = 0, seed = 420))
try_dim("gaussian posterior_predict(draws, newdata = 1 row)",
        posterior_predict(dsg, newdata = dg[1, , drop = FALSE]))
# the rowwise contract the adapter calls: one value per row, a vector
w <- brmsfamily("wiener")
s <- sim_draw(list(family = w, dpars = list(mu = rep(0.8, 5),
                                            bs = rep(1.4, 5),
                                            ndt = rep(0.3, 5),
                                            bias = rep(0.5, 5)),
                   aterms = list(dec = rep(1, 5)), n = 5, autocor = NULL,
                   extra = NULL, max_iter = NULL))
cat("sim_draw(): class", class(s), " dim", if (is.null(dim(s))) "NULL" else
      dim(s), " length", length(s), "\n")
# the same contract for a family brms's helper returns as a vector
g <- sim_draw(list(family = brmsfamily("gaussian"),
                   dpars = list(mu = rep(0, 5), sigma = rep(1, 5)),
                   aterms = list(), n = 5, autocor = NULL, extra = NULL,
                   max_iter = NULL))
cat("sim_draw(gaussian): class", class(g), " dim",
    if (is.null(dim(g))) "NULL" else dim(g), "\n")
